import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import logging
import time
import numpy as np
from sklearn.neighbors import NearestNeighbors

def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Spatial Transcriptomics Analysis Script")
    parser.add_argument("--input_file", required=True, help="Path to the input CSV file")
    parser.add_argument("--output_file_test", required=True, help="Path to the output CSV file for ML and survival analysis results")
    parser.add_argument("--output_file_metadata", required=True, help="Path to the output CSV file for metadata results")
    return parser.parse_args()

def create_output_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

def get_windows(job, n_neighbors):
    '''
    For each region and each individual cell in the dataset, return the indices of the nearest neighbors.

    job:  metadata containing the start time, index of region, region name, indices of region in original dataframe
    n_neighbors:  the number of neighbors to find for each cell
    '''
    start_time, idx, tissue_name, indices = job
    job_start = time.time()

    print("Starting:", str(idx + 1) + '/' + str(len(exps)), ': ' + exps[idx])

    # tissue_group: a grouped data frame with X and Y coordinates grouped by unique tissue regions
    tissue = tissue_group.get_group(tissue_name)

    to_fit = tissue.loc[indices][['x', 'y']].values

    # Unsupervised learner for implementing neighbor searches.
    fit = NearestNeighbors(n_neighbors=n_neighbors).fit(tissue[['x', 'y']].values)

    # Find the nearest neighbors
    m = fit.kneighbors(to_fit) 
    m = m[0], m[1]

    ## sort_neighbors
    args = m[0].argsort(axis=1)
    add = np.arange(m[1].shape[0]) * m[1].shape[1]
    sorted_indices = m[1].flatten()[args + add[:, None]]

    neighbors = tissue.index.values[sorted_indices]

    end_time = time.time()

    print("Finishing:", str(idx + 1) + "/" + str(len(exps)), ": " + exps[idx], end_time - job_start,
          end_time - start_time)
    return neighbors.astype(np.int32)

def plot_neighborhood_frequencies(data2, output_directory):
    fc = data2.groupby(['Image', 'Group']).apply(lambda x: x['neighborhood10'].value_counts(sort=False, normalize=True))
    melt = pd.melt(fc.reset_index(), id_vars=['Image', 'Group'])
    melt = melt.rename(columns={'variable': 'neighborhood', 'value': 'frequency of neighborhood'})
    melt['neighborhood'] = melt['neighborhood'].map({
        0: 'NB 0', 1: 'NB 1', 2: 'NB 2', 3: 'NB 3', 4: 'NB 4',
        5: 'NB 5', 6: 'NB 6', 7: 'NB 7', 8: 'NB 8', 9: 'NB 9'
    })

    f, ax = plt.subplots(figsize=(15, 7))
    sns.pointplot(data=melt, hue='Group', dodge=0.1, join=False, x='neighborhood', y='frequency of neighborhood')
    ax.legend(title="Groups", handletextpad=0, columnspacing=2, loc="upper left", ncol=3, frameon=True)
    plt.xticks(rotation=45)

    frequency_figure_path = os.path.join(output_directory, "NB_freq.svg")
    plt.savefig(frequency_figure_path, dpi=300, bbox_inches="tight")
    logging.info(f"Frequency figure saved to {frequency_figure_path}")

def save_results(df, output_path, description):
    df.to_csv(output_path, index=False)
    logging.info(f"{description} saved to {output_path}")

def main():
    setup_logging()
    args = parse_arguments()
    
    output_directory = create_output_directory("output")
    
    input_file_path = os.path.join(os.getcwd(), args.input_file)
    df = pd.read_csv(input_file_path)

    # Example of extracting features from df
    data2 = df.copy()  # Replace with the actual logic to process df and create data2

    # Grouping data by tissue regions (assuming 'Tissue' is the column name)
    global tissue_group
    tissue_group = df.groupby('Tissue')  # Adjust 'Tissue' to the actual column name

    # Simulate jobs and exps for the sake of example
    jobs = [(time.time(), idx, "Tissue", list(range(100))) for idx in range(5)]  # This is just an example
    global exps
    exps = ["Experiment 1", "Experiment 2", "Experiment 3", "Experiment 4", "Experiment 5"]  # Example experiment names
    
    # Get windows based on nearest neighbors for each job
    all_windows = []
    for job in jobs:
        neighbors = get_windows(job, n_neighbors=5)  # Example with 5 neighbors
        all_windows.append(neighbors)
    
    # Plotting neighborhood frequencies
    plot_neighborhood_frequencies(data2, output_directory)

    # Save output for test results
    output_file_test_path = os.path.join(output_directory, args.output_file_test)
    save_results(data2, output_file_test_path, "Test results")

    # Save metadata results
    data3 = data2[['Object ID', 'neighborhood10']]
    merged_df = pd.merge(df, data3, on='Object ID', how='inner')
    output_file_meta_path = os.path.join(output_directory, args.output_file_metadata)
    save_results(merged_df, output_file_meta_path, "Metadata results")

if __name__ == "__main__":
    main()
