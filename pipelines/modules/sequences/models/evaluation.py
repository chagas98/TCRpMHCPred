#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import argparse
import json
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_curve, roc_auc_score, precision_recall_curve, average_precision_score, accuracy_score


def plot_histogram(y_pred, bins=30, save_dir=""):
    plt.hist(y_pred, bins=bins, edgecolor='black')
    plt.xlabel('Predicted Probability')
    plt.ylabel('Frequency')
    plt.xlim(0, 1)  # Fix x-axis range to [0, 1]
    plt.savefig(f'{save_dir}/hist_probs.png')  # Save the figure to a file
    plt.close()

def plot_loss_curve(fold_losses_train, fold_losses_val, save_dir):
    plt.figure()
    for i, (losses_train, losses_val) in enumerate(zip(fold_losses_train, fold_losses_val)):
        plt.plot(losses_train, label=f'Fold {i+1} Train', linestyle='dashed')
        plt.plot(losses_val, label=f'Fold {i+1} Val')
        
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    plt.ylim(0, 1) 
    plt.savefig(f'{save_dir}/loss_curve.png')
    plt.close()

def plot_confusion_matrix(y_true, y_pred_binary, save_dir):
    cm = confusion_matrix(y_true, y_pred_binary)
    disp = ConfusionMatrixDisplay(cm)
    disp.plot()
    plt.savefig(f'{save_dir}/confusion_matrix.png')
    plt.close()

    return cm

def plot_precision_recall_curve(y_true, y_pred, save_dir):
    """
    Plots the Precision-Recall curve.

    Args:
        y_true (array-like): True labels.
        y_pred (array-like): Predicted probabilities.
        output_path (str): Path to save the plot.
    """
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    avg_precision = average_precision_score(y_true, y_pred)

    plt.figure()
    plt.plot(recall, precision, label=f'AP = {avg_precision:.3f}')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend()
    plt.grid()
    plt.savefig(f'{save_dir}/precision_recall_curve.png')
    plt.close()

def plot_roc_curve(y_true, y_scores, save_dir):
    # Calcula FPR, TPR e thresholds
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    df = pd.DataFrame({'FPR': fpr, 'TPR': tpr, 'Thresholds': thresholds})
    df.to_csv('roc_curve.csv', index=False)

    # Plota a curva ROC
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve')
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig(f'{save_dir}/roc_curve.png')
    plt.close()

def calculate_auc_per_peptide(y_true, y_pred, peptides):
    """
    Calculate the AUC (FPR ≤ 0.1) for each unique peptide.

    Args:
        y_true (array-like): True labels.
        y_pred (array-like): Predicted probabilities.
        peptides (array-like): Peptide identifiers for each sample.

    Returns:
        dict: A dictionary with peptides as keys and their AUC scores as values.
    """
    auc_scores = {}
    seen=set()
    for peptide in peptides:

        if peptide not in seen:
            seen.add(peptide)
        else:
            continue

        mask = (peptides == peptide)

        y_true_peptide = y_true[mask]
        y_pred_peptide = y_pred[mask]

        if len(np.unique(y_true_peptide)) > 1:  # Ensure both classes are present
            auc = roc_auc_score(y_true_peptide, y_pred_peptide, max_fpr=0.1)
            auc_scores[peptide] = auc
        else:
            auc_scores[peptide] = None  # Not enough data to calculate AUC
    
    macro_auc = np.mean([auc for auc in auc_scores.values() if auc is not None])
    
    #parse auc_scores to a DataFrame
    auc_scores_df = pd.DataFrame(list(auc_scores.items()), columns=['Peptide', 'AUC'])
    
    return auc_scores_df, macro_auc

def results_analysis(df_predictions, experiment_name: str = '', test_set_name: str = '', threshold: int = 0.5,
                     description: str = '', model_type: str = '', test_readable: pd.DataFrame = pd.DataFrame(), 
                     train_readable: pd.DataFrame = pd.DataFrame(), save_dir='.'):

    print(f"Results Analysis")
    os.makedirs(save_dir, exist_ok=True)

    #merge the predictions with the readable data by TCR_ID
    print("Merging predictions with readable data...")
    print('Length of df_predictions:', len(df_predictions))
    print('Length of df_readable:', len(test_readable))

    if test_readable.empty or len(test_readable) == 0:
        df_merged = df_predictions
        print('No readable data provided, using predictions only.')
    else:
        df_merged = df_predictions.merge(test_readable, on='TCR_ID', how='left')    
        print('Length of merged df_predictions:', len(df_merged))
    
    if train_readable is not None and not train_readable.empty:
        train_length = len(train_readable)
        n_peptides_train = len(np.unique(train_readable['peptide']))
    else:
        train_length = 'NA'
        n_peptides_train = 'NA'

    # Metrics
    df_predictions['binary_preds'] = (df_predictions['prediction'] >= threshold).astype(int)
    binary_preds = df_predictions['binary_preds'].values
    y_true = df_predictions['Target'].values
    y_pred = df_predictions['prediction'].values
    peptides = df_merged['peptide'].values

    accuracy = accuracy_score(y_true, binary_preds)  
    auc_scores_df, macro_auc = calculate_auc_per_peptide(y_true, y_pred, peptides)
    pep_auc_scores_path = os.path.join(save_dir, 'peptide_auc_scores.csv')
    auc_scores_df.to_csv(pep_auc_scores_path, index=False)

    print(f'Number of unique peptides: {len(np.unique(peptides))}')
    print(f'Test Accuracy: {accuracy:.3f}')        
    print(f'Test Macro AUC (FPR < 0.1): {macro_auc:.3f}')

    #Save results
    plot_histogram(y_pred, bins=30, save_dir=save_dir)
    plot_precision_recall_curve(y_true, y_pred, save_dir)
    plot_roc_curve(y_true, y_pred, save_dir)    
    cm = plot_confusion_matrix(y_true, binary_preds, save_dir)
    tn, fp, fn, tp = cm.ravel()

    result = pd.DataFrame({
        'experiment_name': experiment_name,
        'test_set_name': test_set_name,
        'description': description,
        'model': model_type,
        'threshold': round(threshold, 3),
        'accuracy': round(accuracy, 3),
        'macro_auc': round(macro_auc, 3),
        'tn': int(tn),
        'fn': int(fn),
        'tp': int(tp),
        'fp': int(fp),
        'precision': round(tp / (tp + fp), 3) if (tp + fp) > 0 else 0,
        'recall': round(tp / (tp + fn), 3) if (tp + fn) > 0 else 0,
        'test_length': int(len(df_merged)),
        'npeptides_test': int(len(np.unique(peptides))),
        'train_length': train_length,
        'npeptides_train': n_peptides_train
        }, index=[0])
    
    result.to_csv(f'{save_dir}/results_{model_type}_{experiment_name}.csv', index=False)

    return result

def results_analysis_tulip(df_predictions, experiment_name: str = '', test_set_name: str = '',
                            description: str = '', model_type: str = '', save_dir='.',
                            df_readable: pd.DataFrame = pd.DataFrame(), df_length: int = 0):
    """
    Analyze the results for TULIP model predictions.
    This function assumes that the predictions are already in the correct format.
    """

    print(f"Results Analysis for TULIP")
    os.makedirs(save_dir, exist_ok=True)

    auc_score_df = df_predictions[['peptide', 'AUC0_1']].rename(columns={'AUC0_1': 'AUC', 'peptide': 'Peptide'})
    auc_score_df.to_csv(os.path.join(save_dir, 'peptide_auc_scores.csv'), index=False)
    macro_auc = auc_score_df['AUC'].mean()
    peptides = auc_score_df['Peptide'].values

    result = pd.DataFrame({
        'experiment_name': experiment_name,
        'test_set_name': test_set_name,
        'description': description,
        'model': model_type,
        'threshold': 'NA',
        'accuracy': 'NA',
        'macro_auc': round(macro_auc, 3),
        'tn': 'NA',
        'fn': 'NA',
        'tp': 'NA',
        'fp': 'NA',
        'precision': 'NA',
        'recall': 'NA',
        'test_length': df_length,
        'npeptides_test': int(len(peptides)),
        'train_length': 'NA',
        'npeptides_train': 'NA'
        }, index=[0])
    result.to_csv(f'{save_dir}/results_{model_type}_{experiment_name}.csv', index=False)

    return result


