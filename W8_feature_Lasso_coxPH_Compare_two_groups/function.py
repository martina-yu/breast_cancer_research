import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from sklearn.utils import resample

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sksurv.compare import compare_survival
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored, concordance_index_ipcw
from sksurv.metrics import brier_score
from sklearn.model_selection import StratifiedKFold
from lifelines.statistics import logrank_test
from sksurv.svm import FastSurvivalSVM
from sklearn.model_selection import GridSearchCV, ShuffleSplit
from sklearn.metrics import make_scorer
from sksurv.ensemble import GradientBoostingSurvivalAnalysis

class EarlyStoppingMonitor:
    def __init__(self, window_size, max_iter_without_improvement):
        self.window_size = window_size
        self.max_iter_without_improvement = max_iter_without_improvement
        self._best_step = -1

    def __call__(self, iteration, estimator, args):
        # continue training for first self.window_size iterations
        if iteration < self.window_size:
            return False

        # compute average improvement in last self.window_size iterations.
        # oob_improvement_ is the different in negative log partial likelihood
        # between the previous and current iteration.
        start = iteration - self.window_size + 1
        end = iteration + 1
        improvement = np.mean(estimator.oob_improvement_[start:end])

        if improvement > 1e-6:
            self._best_step = iteration
            return False  # continue fitting

        # stop fitting if there was no improvement
        # in last max_iter_without_improvement iterations
        diff = iteration - self._best_step
        return diff >= self.max_iter_without_improvement


def read_file(survival_data_path, lasso_file_path, comp_file_path):
    train_df_control = pd.read_csv(survival_data_path, index_col=0)
    gene_exp = train_df_control.iloc[:, :-10]
    y = train_df_control.iloc[:, -10:]

    print(gene_exp.shape, y.shape)

    clinical_val = train_df_control[['Menospausal Status', 'Age', 'Risk Group']]
    survival_val = train_df_control[['USI','Time to event if any (days)','IDFS Event']]
    print(clinical_val.shape, survival_val.shape)
    print(y.columns)

    file_path_lasso = lasso_file_path
    df_lasso = pd.read_excel(file_path_lasso)
    print(df_lasso.columns)

    file_path_comp = comp_file_path
    try:
        all_sheets = pd.read_excel(file_path_comp, sheet_name=None)
        
        print("Available sheets:", list(all_sheets.keys()))

        for sheet_name, comp_result in all_sheets.items():
            print(f"\n--- Reading sheet: '{sheet_name}' ---")
            print(comp_result.head())

    except FileNotFoundError:
        print(f"Error: The file '{file_path_comp}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
    gene_lasso = df_lasso['gene_name'][:6].to_list()
    gene_comp = comp_result['Gene'].tolist()
    print(gene_lasso, gene_comp)

    selected_genes_lasso = gene_exp[gene_lasso]
    selected_genes_comp = gene_exp[gene_comp]

    print(selected_genes_lasso.shape,selected_genes_comp.shape)
    
    return y, clinical_val, survival_val, gene_exp, df_lasso, comp_result, selected_genes_lasso, selected_genes_comp

def quantile_survival_plot(y_test, rs_df):
    quantiles = np.quantile(rs_df, [0.25, 0.5, 0.75])

    group_labels = [
        'Low Risk (Q1)', 
        'Medium-Low Risk (Q2)', 
        'Medium-High Risk (Q3)', 
        'High Risk (Q4)'
    ]
    masks = [
        rs_df <= quantiles[0],
        (rs_df > quantiles[0]) & (rs_df <= quantiles[1]),
        (rs_df > quantiles[1]) & (rs_df <= quantiles[2]),
        rs_df > quantiles[2]
    ]

    plt.figure(figsize=(10, 6))
    kmf = KaplanMeierFitter()

    for i, mask in enumerate(masks):
        if np.any(mask):
            time_data = y_test[mask]['time'].flatten()
            event_data = y_test[mask]['event'].flatten()
            
            kmf.fit(time_data, event_data, label=group_labels[i])
            kmf.plot_survival_function(ci_show=False)
        else:
            print(f"Skipping empty group: {group_labels[i]}")

    plt.title("Survival Curves based on Model's Risk Score Quartiles")
    plt.xlabel("Time to Event")
    plt.ylabel("Survival Probability")
    plt.legend()
    plt.grid(True)
    plt.show()
    
def quantile_survival_plot_for_folds(y_test, risk_scores_test, fold_num):
    quantiles = np.quantile(risk_scores_test, [0.25, 0.5, 0.75])

    group_labels = [
        'Low Risk (Q1)', 
        'Medium-Low Risk (Q2)', 
        'Medium-High Risk (Q3)', 
        'High Risk (Q4)'
    ]
    masks = [
        risk_scores_test <= quantiles[0],
        (risk_scores_test > quantiles[0]) & (risk_scores_test <= quantiles[1]),
        (risk_scores_test > quantiles[1]) & (risk_scores_test <= quantiles[2]),
        risk_scores_test > quantiles[2]
    ]

    plt.figure(figsize=(10, 6))
    kmf = KaplanMeierFitter()

    for i, mask in enumerate(masks):
        if np.any(mask):
            time_data = y_test[mask]['time'].flatten()
            event_data = y_test[mask]['event'].flatten()
            
            kmf.fit(time_data, event_data, label=group_labels[i])
            kmf.plot_survival_function(ci_show=False)
        else:
            print(f"Skipping empty group: {group_labels[i]}")

    plt.title(f"Survival Curves for Risk Quartiles (Fold {fold_num+1})")
    plt.xlabel("Time to Event")
    plt.ylabel("Survival Probability")
    plt.legend()
    plt.grid(True)
    plt.show()
    
def plot_km_curves(y_train, y_test, fold_num):
    """
    Plots Kaplan-Meier curves for the training and testing sets of a fold.

    Args:
        y_train (np.ndarray): The structured training data array.
        y_test (np.ndarray): The structured testing data array.
        fold_num (int): The current fold number.
    """
    kmf = KaplanMeierFitter()
    
    plt.figure(figsize=(8, 6))

    kmf.fit(y_train['time'].flatten(), event_observed=y_train['event'].flatten(), label='Training Set')
    kmf.plot_survival_function(ci_show=False)

    kmf.fit(y_test['time'].flatten(), event_observed=y_test['event'].flatten(), label='Testing Set')
    kmf.plot_survival_function(ci_show=False)
    
    plt.title(f'Kaplan-Meier Curves (Fold {fold_num + 1})')
    plt.xlabel('Time to Event')
    plt.ylabel('Survival Probability')
    plt.legend()
    plt.grid(True)
    plt.show()
    
def gene_exp_scale(selected_genes, clinical_val, gene_exp):
    for_scale = pd.concat([selected_genes,clinical_val['Age']],axis=1)
    print(for_scale.shape)

    df = gene_exp.copy()
    total_zeroes = (df == 0).sum().sum()
    total_el = df.size
    sparsity = total_zeroes/total_el * 100
    print(f'sparsity is totally {sparsity}% among matrix')

    selected_df = selected_genes.copy()
    total_zeroes = (selected_df == 0).sum().sum()
    total_el = selected_df.size
    sparsity = total_zeroes/total_el * 100
    print(f'sparsity is totally {sparsity}% among selected genes matrix')
    
    scaler = StandardScaler()
    # selected_genes_scaled = pd.DataFrame(scaler.fit_transform(selected_genes), columns=selected_genes.columns, index = selected_genes.index)
    age_gene_scaled = pd.DataFrame(scaler.fit_transform(for_scale), columns=for_scale.columns, index=for_scale.index)
    
    # selected_genes_scaled = pd.concat([selected_genes_scaled, age_gene_scaled], axis=1)
    
    X = pd.concat([age_gene_scaled, clinical_val[['Risk Group','Menospausal Status']]],axis=1)
    X = pd.get_dummies(X, columns=['Menospausal Status','Risk Group'], drop_first=True)

    print(X.columns)
    
    return X


def split_dataset(X, y_structured):
    X_train, X_test, y_train, y_test = train_test_split(X, y_structured, test_size=0.3, random_state=100, stratify=y_structured['event'])
    
    print('log-rank test for survival curve')
    y_all = np.concatenate([y_train, y_test])
    group_labels = np.array(['training'] * len(y_train) + ['testing'] * len(y_test))
    results = compare_survival(y_all, group_labels)
    print(f'log-rank test statistic, p-value: {results}')

    print("--- training set distribution ---")
    event_counts_train = np.bincount(y_train['event'])
    print(event_counts_train)
    print(f"False: {event_counts_train[0]}, True: {event_counts_train[1]}")
    print(f"event(happen) rate: {event_counts_train[1] / len(y_train) * 100:.2f}%")

    print("\n--- testing set distribution ---")
    event_counts_test = np.bincount(y_test['event'])
    print(event_counts_test)
    print(f"False: {event_counts_test[0]}, True: {event_counts_test[1]}")
    print(f"event(happen) rate: {event_counts_train[1] / len(y_train) * 100:.2f}%")
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6))

    time_train, survival_prob_train = kaplan_meier_estimator(y_train['event'], y_train['time'])
    ax1.step(time_train, survival_prob_train, where="post")
    ax1.set_title("training survival curve plot")
    ax1.set_xlabel("time(day)")
    ax1.set_ylabel("survival proportion")
    ax1.grid(True)

    time_test, survival_prob_test = kaplan_meier_estimator(y_test['event'], y_test['time'])
    ax2.step(time_test, survival_prob_test, where="post")
    ax2.set_title("testing survival curve plot")
    ax2.set_xlabel("time(day)")
    ax2.set_ylabel("survival proportion")
    ax2.grid(True)

    df_y = pd.DataFrame(y_structured)
    time_df, survival_prob_df = kaplan_meier_estimator(df_y['event'], df_y['time'])
    ax3.step(time_df, survival_prob_df, where="post")
    ax3.set_title("All dataset survival curve plot")
    ax3.set_xlabel("time(day)")
    ax3.set_ylabel("survival proportion")
    ax3.grid(True)

    plt.tight_layout()
    plt.show()
    
    return X_train, X_test, y_train, y_test

def random_survival_forest_with_sd(X_train, X_test, y_train, y_test, n, n_bootstraps=100):
    rsf_model = RandomSurvivalForest(n_estimators=100, random_state=42)
    rsf_model.fit(X_train, y_train)

    time_point = 365 * n

    c_indices = []
    auc_scores = []
    brier_scores = []

    for i in range(n_bootstraps):
        indices = resample(np.arange(len(X_test)), replace=True, random_state=i)
        X_test_boot = X_test.iloc[indices]
        y_test_boot = y_test[indices]
        
        # Use a try-except block to gracefully handle potential errors
        try:
            risk_scores_test = rsf_model.predict(X_test_boot)
            surv_funcs_test = rsf_model.predict_survival_function(X_test_boot)
            
            c_index = concordance_index_censored(y_test_boot['event'], y_test_boot['time'], risk_scores_test)[0]
            c_indices.append(c_index)
            
            auc_test = concordance_index_ipcw(y_train, y_test_boot, risk_scores_test, tau=time_point)[0]
            auc_scores.append(auc_test)

            surv_prob_test = np.array([f(time_point) for f in surv_funcs_test])
            score_test = brier_score(y_train, y_test_boot, surv_prob_test, time_point)[1]
            brier_scores.append(score_test)
            
        except ValueError as e:
            # Skip this iteration if a ValueError occurs
            print(f"Skipping bootstrap sample {i} due to an error: {e}")
            continue

    mean_c_index = np.mean(c_indices)
    std_c_index = np.std(c_indices)
    print(f"C-index (mean ± std): {mean_c_index:.4f} ± {std_c_index:.4f}")

    mean_auc = np.mean(auc_scores)
    std_auc = np.std(auc_scores)
    print(f"Time-Dependent AUC in {time_point} days (mean ± std): {mean_auc:.4f} ± {std_auc:.4f}")

    mean_brier_score = np.mean(brier_scores)
    std_brier_score = np.std(brier_scores)
    print(f"Brier Score in {time_point} days (mean ± std): {mean_brier_score:.4f} ± {std_brier_score:.4f}")

    risk_scores_original = rsf_model.predict(X_test)
    quantile_survival_plot(y_test, risk_scores_original)
    
    return

# def random_survival_forest(X_train, X_test, y_train, y_test, n):
    # rsf_model = RandomSurvivalForest(n_estimators=100, random_state=42)
    # rsf_model.fit(X_train, y_train)
    
    # risk_scores_test = rsf_model.predict(X_test)
    # surv_funcs_train = rsf_model.predict_survival_function(X_train)
    # surv_funcs_test = rsf_model.predict_survival_function(X_test)

    # print(f"Number of survival functions for the training set: {len(surv_funcs_train)}")
    # print(f"Number of survival functions for the testing set: {len(surv_funcs_test)}")

    # c_index_test = concordance_index_censored(y_test['event'], y_test['time'], risk_scores_test)[0]
    # print(f"C-index for the testing set: {c_index_test:.4f}")

    # time_point = 365 * n

    # print('------ Time-Dependent AUC ------')
    # auc_test = concordance_index_ipcw(
    #     y_train, y_test, risk_scores_test, tau=time_point
    # )[0]

    # print(f"testing set in {time_point} days has Time-Dependent AUC: {auc_test:.4f}")

    # print('------ Brier Score ------')
    # surv_funcs_test = rsf_model.predict_survival_function(X_test)
    # surv_prob_test = np.array([f(time_point) for f in surv_funcs_test])
    # score_test = brier_score(y_train, y_test, surv_prob_test, time_point)[1]
    # print(f"testing set in {time_point} days has Brier Score: {score_test}")

    # quantile_survival_plot(y_test, risk_scores_test)
    
    # return

def random_survival_forest_for_all(X, y_structured, n, n_splits):

    rsf_model = RandomSurvivalForest(n_estimators=100, random_state=42)
    skf = StratifiedKFold(n_splits, shuffle=True, random_state=42)

    c_index_scores = []
    auc_scores = []
    brier_scores = []

    time_point = 365 * n

    for fold_num, (train_index, test_index) in enumerate(skf.split(X, y_structured['event'])):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y_structured[train_index], y_structured[test_index]
        
        plot_km_curves(y_train, y_test, fold_num)
        
        print(f"--- Fold {fold_num + 1} ---")
        print("Performing log-rank test for survival curve similarity...")

        time_train = y_train['time'].flatten()
        event_train = y_train['event'].flatten()
        time_test = y_test['time'].flatten()
        event_test = y_test['event'].flatten()

        results = logrank_test(time_train, time_test, event_observed_A=event_train, event_observed_B=event_test)
        
        print(f'Log-rank test statistic: {results.test_statistic:.4f}, p-value: {results.p_value:.4f}')

        if results.p_value < 0.05:
            print("Warning: The survival curves of the training and testing sets are significantly different (p < 0.05).")
        else:
            print("The survival curves of the training and testing sets are not significantly different (p > 0.05).")
        print("-" * 30)
        
        train_event_rate = y_train['event'].mean()
        test_event_rate = y_test['event'].mean()
        print(f"  Training Set: Event Rate = {train_event_rate:.4f} ({len(y_train)} samples)")
        print(f"  Testing Set:  Event Rate = {test_event_rate:.4f} ({len(y_test)} samples)")
        print("-" * 30)
        
        rsf_model.fit(X_train, y_train)
        risk_scores_test = rsf_model.predict(X_test)
        surv_funcs_test = rsf_model.predict_survival_function(X_test)

        c_index_fold = concordance_index_censored(
            y_test['event'], y_test['time'], risk_scores_test
        )[0]
        print(f'c-index in {fold_num+1} is {c_index_fold}')
        c_index_scores.append(c_index_fold)
        
        auc_fold = concordance_index_ipcw(y_train, y_test, risk_scores_test, tau=time_point)[0]
        auc_scores.append(auc_fold)
        
        surv_prob_test = np.array([f(time_point) for f in surv_funcs_test])
        score_fold = brier_score(y_train, y_test, surv_prob_test, time_point)[1]
        brier_scores.append(score_fold)

        quantile_survival_plot_for_folds(y_test, risk_scores_test, fold_num)
        
    average_c_index = np.mean(c_index_scores)
    std_c_index = np.std(c_index_scores)

    print(f"C-index for each fold: {c_index_scores}")
    print(f"Average C-index (4-fold cross-validation): {average_c_index:.4f}")
    print(f"Standard deviation of C-index: {std_c_index:.4f}")
    print(f"Time-Dependent AUC ({time_point} days): {np.mean(auc_scores):.4f} +/- {np.std(auc_scores):.4f}")
    print(f"Brier Score ({time_point} days): {np.mean(brier_scores):.4f} +/- {np.std(brier_scores):.4f}")
    return

def c_index_scorer(estimator, X, y):
    """A custom scorer function for concordance_index_censored."""
    risk_scores = estimator.predict(X)
    return concordance_index_censored(y['event'], y['time'], risk_scores)[0]


def Survival_SVM_for_split(X_train, X_test, y_train, y_test, n, n_bootstraps=100):
    """
    Fits a Survival SVM model on pre-split data and evaluates it with
    C-index, Time-Dependent AUC, and Brier Score using bootstrapping.
    
    Args:
        X_train, X_test, y_train, y_test: Training and testing data.
        n (int): Number of years for time-dependent metrics (e.g., n=1 for 365 days).
        n_bootstraps (int): Number of bootstrap samples to use for standard deviation.
    """
    print('------ Running Survival SVM on a single split ------')
    
    time_point = 365 * n
    
    model = FastSurvivalSVM(max_iter=1000, tol=1e-5, random_state=0)
    param_grid = {"alpha": [2.0**v for v in range(-12, 13, 2)]}
    cv = ShuffleSplit(n_splits=10, test_size=0.3, random_state=0)
    
    gcv = GridSearchCV(model, param_grid, scoring=c_index_scorer, n_jobs=1, refit=True, cv=cv)
    
    gcv.fit(X_train, y_train)
    
    best_alpha = gcv.best_params_['alpha']
    print(f"Best alpha found: {best_alpha}")

    final_model = gcv.best_estimator_

    c_indices = []
    auc_scores = []
    brier_scores = []

    for i in range(n_bootstraps):
        indices = resample(np.arange(len(X_test)), replace=True, random_state=i)
        X_test_boot = X_test.iloc[indices]
        y_test_boot = y_test[indices]
        
        try:
            risk_scores_test = final_model.predict(X_test_boot)
            
            c_index = concordance_index_censored(y_test_boot['event'], y_test_boot['time'], risk_scores_test)[0]
            c_indices.append(c_index)
            
            auc_test = concordance_index_ipcw(y_train, y_test_boot, risk_scores_test, tau=time_point)[0]
            auc_scores.append(auc_test)
            
        except ValueError as e:
            print(f"Skipping bootstrap sample {i} due to an error: {e}")
            continue

    # Calculate and print the mean and standard deviation
    print('--- Metrics for a single split ---')
    print(f"C-index (mean ± std): {np.mean(c_indices):.4f} ± {np.std(c_indices):.4f}")
    print(f"Time-Dependent AUC in {time_point} days (mean ± std): {np.mean(auc_scores):.4f} ± {np.std(auc_scores):.4f}")
    
    return


def Survival_SVM_for_all(X, y_structured, n, n_splits, n_bootstraps=100):
    """
    Performs stratified k-fold cross-validation with a Survival SVM model,
    calculating metrics with bootstrapping for standard deviation.

    Args:
        X (pd.DataFrame): All feature data.
        y_structured (np.ndarray): All structured target data (event, time).
        n (int): Number of years for time-dependent metrics.
        n_splits (int): Number of folds for cross-validation.
        n_bootstraps (int): Number of bootstrap samples for standard deviation.
    """
    print(f'------ Running Survival SVM with {n_splits}-fold Stratified CV ------')
    
    time_point = 365 * n
    
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    c_index_means = []
    c_index_stds = []
    auc_means = []
    auc_stds = []

    for fold_num, (train_index, test_index) in enumerate(skf.split(X, y_structured['event'])):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y_structured[train_index], y_structured[test_index]
        
        print(f"\n--- Fold {fold_num + 1} ---")
        
        results = logrank_test(y_train['time'], y_test['time'], y_train['event'], y_test['event'])
        print(f'Log-rank test p-value: {results.p_value:.4f}')
        print(f"  Training Set: Event Rate = {y_train['event'].mean():.4f} ({len(y_train)} samples)")
        print(f"  Testing Set:  Event Rate = {y_test['event'].mean():.4f} ({len(y_test)} samples)")
        print("-" * 30)

        # Fit model on this fold's training data
        model = FastSurvivalSVM(alpha=0.1, max_iter=1000, tol=1e-5, random_state=0)
        model.fit(X_train, y_train)
        
        c_indices_fold = []
        auc_scores_fold = []
        
        for j in range(n_bootstraps):
            indices = resample(np.arange(len(X_test)), replace=True, random_state=j)
            X_test_boot = X_test.iloc[indices]
            y_test_boot = y_test[indices]
            
            try:
                risk_scores_test = model.predict(X_test_boot)
                
                c_index = concordance_index_censored(
                    y_test_boot['event'], y_test_boot['time'], risk_scores_test
                )[0]
                c_indices_fold.append(c_index)
                
                auc_test = concordance_index_ipcw(y_train, y_test_boot, risk_scores_test, tau=time_point)[0]
                auc_scores_fold.append(auc_test)
                
            except ValueError as e:
                print(f"  Skipping bootstrap sample {j} for Fold {fold_num+1} due to error: {e}")
                continue
        
        c_index_means.append(np.mean(c_indices_fold))
        c_index_stds.append(np.std(c_indices_fold))
        auc_means.append(np.mean(auc_scores_fold))
        auc_stds.append(np.std(auc_scores_fold))
        
        print(f"C-index for Fold {fold_num+1}: {c_index_means[-1]:.4f} ± {c_index_stds[-1]:.4f}")
        print(f"AUC ({time_point} days) for Fold {fold_num+1}: {auc_means[-1]:.4f} ± {auc_stds[-1]:.4f}")
        
    print("\n" + "=" * 50)
    print("Final Cross-Validation Results")
    print("=" * 50)
    print(f"Overall C-index (Mean ± STD): {np.mean(c_index_means):.4f} ± {np.std(c_index_means):.4f}")
    print(f"Overall AUC ({time_point} days) (Mean ± STD): {np.mean(auc_means):.4f} ± {np.std(auc_means):.4f}")
    
    

def survival_xgboost_for_split(X_train, X_test, y_train, y_test, n, n_bootstraps=100):
    """
    Fits a Gradient Boosting Survival Analysis model on pre-split data and evaluates it with
    C-index, Time-Dependent AUC, and Brier Score using bootstrapping for standard deviation.
    
    Args:
        X_train (pd.DataFrame): Training feature data.
        X_test (pd.DataFrame): Testing feature data.
        y_train (np.ndarray): Training structured target data (event, time).
        y_test (np.ndarray): Testing structured target data (event, time).
        n (int): Number of years for time-dependent metrics.
        n_bootstraps (int): Number of bootstrap samples to use for standard deviation.
    """
    print('------ Running Gradient Boosting Survival Analysis on a single split ------')
    
    time_point = 365 * n
    
    # Fit the model on the full training set
    model = GradientBoostingSurvivalAnalysis(n_estimators=100, learning_rate=0.1, max_depth=1, random_state=0)
    model.fit(X_train, y_train)
    
    # Arrays to store bootstrapped results
    c_indices = []
    auc_scores = []
    brier_scores = []

    for i in range(n_bootstraps):
        # Resample the test set with replacement
        indices = resample(np.arange(len(X_test)), replace=True, random_state=i)
        X_test_boot = X_test.iloc[indices]
        y_test_boot = y_test[indices]
        
        # Check if the bootstrapped sample has any events
        if not np.any(y_test_boot['event']):
            print(f"Skipping bootstrap sample {i}: All samples are censored.")
            continue
            
        try:
            risk_scores_test = model.predict(X_test_boot)
            
            c_index = concordance_index_censored(
                y_test_boot['event'], y_test_boot['time'], risk_scores_test
            )[0]
            c_indices.append(c_index)
            
            auc_test = concordance_index_ipcw(y_train, y_test_boot, risk_scores_test, tau=time_point)[0]
            auc_scores.append(auc_test)

            surv_funcs_test = model.predict_survival_function(X_test_boot)
            surv_prob_test = np.array([f(time_point) for f in surv_funcs_test])
            score_test = brier_score(y_train, y_test_boot, surv_prob_test, time_point)[1]
            brier_scores.append(score_test)
            
        except ValueError as e:
            # Catch other potential ValueErrors from sksurv metrics
            print(f"Skipping bootstrap sample {i} due to an error: {e}")
            continue

    # Calculate and print the mean and standard deviation
    print('\n--- Metrics for a single split ---')
    print(f"C-index (mean ± std): {np.mean(c_indices):.4f} ± {np.std(c_indices):.4f}")
    print(f"Time-Dependent AUC in {time_point} days (mean ± std): {np.mean(auc_scores):.4f} ± {np.std(auc_scores):.4f}")
    print(f"Brier Score in {time_point} days (mean ± std): {np.mean(brier_scores):.4f} ± {np.std(brier_scores):.4f}")
    
    return


def survival_xgboost_for_all_sample(X, y_structured, n, n_splits, n_bootstraps=100):
    """
    Performs stratified k-fold cross-validation with a Gradient Boosting Survival Analysis model,
    calculating metrics with bootstrapping for standard deviation.

    Args:
        X (pd.DataFrame): All feature data.
        y_structured (np.ndarray): All structured target data (event, time).
        n (int): Number of years for time-dependent metrics.
        n_splits (int): Number of folds for cross-validation.
        n_bootstraps (int): Number of bootstrap samples for standard deviation.
    """
    print(f'------ Running Gradient Boosting Survival Analysis with {n_splits}-fold Stratified CV ------')
    
    time_point = 365 * n
    
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    # Lists to store the mean and std of metrics for each fold
    c_index_means = []
    c_index_stds = []
    auc_means = []
    auc_stds = []
    brier_means = []
    brier_stds = []

    for fold_num, (train_index, test_index) in enumerate(skf.split(X, y_structured['event'])):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y_structured[train_index], y_structured[test_index]
        
        print(f"\n--- Fold {fold_num + 1} ---")
        
        # Log-rank test to compare survival curves
        results = logrank_test(y_train['time'], y_test['time'], y_train['event'], y_test['event'])
        print(f'Log-rank test p-value: {results.p_value:.4f}')
        print(f"  Training Set: Event Rate = {y_train['event'].mean():.4f} ({len(y_train)} samples)")
        print(f"  Testing Set:  Event Rate = {y_test['event'].mean():.4f} ({len(y_test)} samples)")
        print("-" * 30)

        # Fit model on this fold's training data
        model = GradientBoostingSurvivalAnalysis(n_estimators=100, learning_rate=0.1, max_depth=1, random_state=0)
        model.fit(X_train, y_train)
        
        c_indices_fold = []
        auc_scores_fold = []
        brier_scores_fold = []
        
        for j in range(n_bootstraps):
            indices = resample(np.arange(len(X_test)), replace=True, random_state=j)
            X_test_boot = X_test.iloc[indices]
            y_test_boot = y_test[indices]
            
            # Check if the bootstrapped sample has any events
            if not np.any(y_test_boot['event']):
                # In this case, just skip and print a message
                print(f"  Skipping bootstrap sample {j} for Fold {fold_num+1}: All samples are censored.")
                continue

            try:
                risk_scores_test = model.predict(X_test_boot)
                
                c_index = concordance_index_censored(
                    y_test_boot['event'], y_test_boot['time'], risk_scores_test
                )[0]
                c_indices_fold.append(c_index)
                
                auc_test = concordance_index_ipcw(y_train, y_test_boot, risk_scores_test, tau=time_point)[0]
                auc_scores_fold.append(auc_test)
                
                surv_funcs_test = model.predict_survival_function(X_test_boot)
                surv_prob_test = np.array([f(time_point) for f in surv_funcs_test])
                score_test = brier_score(y_train, y_test_boot, surv_prob_test, time_point)[1]
                brier_scores_fold.append(score_test)
                
            except ValueError as e:
                # Catch other potential ValueErrors from sksurv metrics
                print(f"  Skipping bootstrap sample {j} for Fold {fold_num+1} due to error: {e}")
                continue
        
        # Calculate mean and std for the current fold and store them
        c_index_means.append(np.mean(c_indices_fold))
        c_index_stds.append(np.std(c_indices_fold))
        auc_means.append(np.mean(auc_scores_fold))
        auc_stds.append(np.std(auc_scores_fold))
        brier_means.append(np.mean(brier_scores_fold))
        brier_stds.append(np.std(brier_scores_fold))
        
        print(f"C-index for Fold {fold_num+1}: {c_index_means[-1]:.4f} ± {c_index_stds[-1]:.4f}")
        print(f"AUC ({time_point} days) for Fold {fold_num+1}: {auc_means[-1]:.4f} ± {auc_stds[-1]:.4f}")
        print(f"Brier Score ({time_point} days) for Fold {fold_num+1}: {brier_means[-1]:.4f} ± {brier_stds[-1]:.4f}")
        
    # Final output: average of means and stds across all folds
    print("\n" + "=" * 50)
    print("Final Cross-Validation Results")
    print("=" * 50)
    print(f"Overall C-index (Mean ± STD): {np.mean(c_index_means):.4f} ± {np.std(c_index_means):.4f}")
    print(f"Overall AUC ({time_point} days) (Mean ± STD): {np.mean(auc_means):.4f} ± {np.std(auc_means):.4f}")
    print(f"Overall Brier Score ({time_point} days) (Mean ± STD): {np.mean(brier_means):.4f} ± {np.std(brier_means):.4f}")
    
    return