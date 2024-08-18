import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process data and perform analysis.')
    parser.add_argument('--input', type=str, required=True, help='Input CSV file name')
    parser.add_argument('--output', type=str, required=True, help='Output CSV file name')
    return parser.parse_args()

def load_and_preprocess_data(input_file):
    meta = pd.read_csv(input_file, index_col=0)
    meta = meta[meta['CT_final'] != "Artifact"]
    return meta

def calculate_frequencies(meta):
    ct_freq = pd.crosstab(index=meta['Image'], columns=meta['CT_final'], normalize="index")
    all_freqs = pd.crosstab(index=[meta['neighborhood10'], meta['Image']], columns=meta['CT_final'], normalize="index")
    return ct_freq, all_freqs

def calculate_group_counts(meta):
    df_unique_images = meta[['Image', 'Group']].drop_duplicates()
    return df_unique_images['Group'].value_counts()

def normalize(df):
    return pd.DataFrame(df.values / np.sum(df.values, axis=1, keepdims=True), index=df.index, columns=df.columns).fillna(0)

def perform_statistical_analysis(ct_freq, all_freqs, group):
    changes = {}
    cells = ct_freq.columns
    nbs = all_freqs.reset_index()['neighborhood10'].unique()

    for col in cells:
        for nb in nbs:
            X = pd.concat([ct_freq[col], group.astype('int'), pd.Series(np.ones(len(group)), index=group.index.values)], axis=1).values
            if f'{col}_{nb}' in all_freqs.columns:
                Y = all_freqs[f'{col}_{nb}'].values
                X, Y = X[~pd.isna(Y)], Y[~pd.isna(Y)]
                results = sm.OLS(Y, X).fit()
                changes[(col, nb)] = (results.pvalues[1], results.params[1])

    return pd.DataFrame(changes).loc[1].unstack(), pd.DataFrame(changes).loc[0].unstack()

def generate_heatmap(dat, pvals, output_filename):
    plt.figure(figsize=(12, 6))
    sns.heatmap(dat, cmap='coolwarm', vmin=-1, vmax=1, cbar=False)

    for a, b in zip(*np.where(pvals < 0.1)):
        plt.text(b + 0.5, a + 0.55, '*', fontsize=30, ha='center', va='center')
    
    plt.tight_layout()
    plt.savefig(output_filename, format="svg", dpi=400, bbox_inches="tight")

def main():
    args = parse_arguments()

    # Load and preprocess data
    meta = load_and_preprocess_data(args.input)
    
    # Calculate frequencies
    ct_freq, all_freqs = calculate_frequencies(meta)
    
    # Calculate group counts
    group_counts = calculate_group_counts(meta)
    
    # Map groups to images
    map_group = dict(zip(meta['Image'], meta['Group']))
    all_freqs['group'] = all_freqs.reset_index()['Image'].map(map_group).astype('category')
    all_freqs = all_freqs[all_freqs['group'] != 'Unknown']
    
    # Perform statistical analysis
    group = meta.set_index('Image')['Group'].map({'Responder': 1, 'Non_Responder': 0})
    dat, pvals = perform_statistical_analysis(ct_freq, all_freqs, group)
    
    # Generate heatmap
    output_svg = Path(args.output).with_suffix('.svg')
    generate_heatmap(dat, pvals, output_svg)

    # Save processed data
    all_freqs.to_csv(args.output)

if __name__ == "__main__":
    main()
