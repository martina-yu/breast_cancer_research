import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from sklearn.preprocessing import StandardScaler

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
    selected_genes_scaled = pd.DataFrame(scaler.fit_transform(selected_genes), columns=selected_genes.columns, index = selected_genes.index)
    age_gene_scaled = pd.DataFrame(scaler.fit_transform(for_scale), columns=for_scale.columns, index=for_scale.index)
    
    X = pd.concat([selected_genes_scaled, clinical_val[['Risk Group','Menospausal Status']]],axis=1)
    X = pd.get_dummies(X, columns=['Menospausal Status','Risk Group'], drop_first=True)

    print(X.columns)
    
    return X
