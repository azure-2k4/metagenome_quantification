#!/usr/bin/env python3
"""
Generate visualization plots for method comparison.
"""
import json
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def plot_abundance_comparison(ground_truth, predictions, output_file):
    """
    Create bar plot comparing true vs predicted abundances.
    
    Args:
        ground_truth (pd.Series): Ground truth abundances
        predictions (pd.DataFrame): Predicted abundances from all methods
        output_file (str): Path to save plot
    """
    # Prepare data
    plot_data = pd.DataFrame({
        'Genome': ground_truth.index,
        'Ground Truth': ground_truth.values
    })
    
    for method in predictions.columns:
        if method != 'genome':
            plot_data[method] = predictions[method].values
    
    # Melt data for plotting
    plot_data_melted = pd.melt(
        plot_data,
        id_vars=['Genome'],
        var_name='Method',
        value_name='Abundance'
    )
    
    # Create plot
    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=plot_data_melted,
        x='Genome',
        y='Abundance',
        hue='Method'
    )
    
    plt.xticks(rotation=45)
    plt.title('True vs Predicted Abundances by Method')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_performance_comparison(runtime_data, output_file):
    """
    Create performance comparison plot.
    
    Args:
        runtime_data (pd.DataFrame): Runtime and memory usage data
        output_file (str): Path to save plot
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Runtime plot
    sns.barplot(
        data=runtime_data,
        x='method',
        y='runtime_minutes',
        ax=ax1
    )
    ax1.set_title('Runtime Comparison')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
    ax1.set_ylabel('Runtime (minutes)')
    
    # Memory plot
    sns.barplot(
        data=runtime_data,
        x='method',
        y='memory_gb',
        ax=ax2
    )
    ax2.set_title('Memory Usage Comparison')
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
    ax2.set_ylabel('Memory Usage (GB)')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_accuracy_metrics(metrics_data, output_file):
    """
    Create summary plot of accuracy metrics.
    
    Args:
        metrics_data (pd.DataFrame): Accuracy metrics for each method
        output_file (str): Path to save plot
    """
    # Prepare data
    metrics_melted = pd.melt(
        metrics_data,
        id_vars=['method'],
        value_vars=['mae', 'rmse', 'correlation', 'aupr'],
        var_name='Metric',
        value_name='Value'
    )
    
    # Create plot
    g = sns.FacetGrid(
        metrics_melted,
        col='Metric',
        col_wrap=2,
        height=4,
        aspect=1.5
    )
    
    g.map_dataframe(
        sns.barplot,
        x='method',
        y='Value'
    )
    
    # Rotate x labels
    for ax in g.axes:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main(snakemake):
    """
    Main function to create visualization plots.
    
    Args:
        snakemake: Snakemake object containing parameters
    """
    # Load data
    try:
        metrics = pd.read_csv(snakemake.input[0])
    except pd.errors.EmptyDataError:
        logging.error(f"Metrics file {snakemake.input[0]} is empty. Skipping plot generation.")
        return

    try:
        runtime = pd.read_csv(snakemake.input[1])
    except pd.errors.EmptyDataError:
        logging.error(f"Runtime file {snakemake.input[1]} is empty. Skipping plot generation.")
        return

    # Only create plots if data is not empty
    if not metrics.empty:
        if 'true_abundance' in metrics.columns and 'genome' in metrics.columns and 'predicted_abundance' in metrics.columns:
            plot_abundance_comparison(
                metrics['true_abundance'],
                metrics[['genome', 'predicted_abundance']],
                snakemake.output.abundance
            )
        plot_accuracy_metrics(
            metrics,
            snakemake.output.accuracy
        )
    else:
        logging.error("Metrics data is empty. Skipping abundance and accuracy plots.")

    if not runtime.empty:
        plot_performance_comparison(
            runtime,
            snakemake.output.performance
        )
    else:
        logging.error("Runtime data is empty. Skipping performance plot.")

    logging.info("Visualization completed successfully")

if __name__ == "__main__":
    main(snakemake)
