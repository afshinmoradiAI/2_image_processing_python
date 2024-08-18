import os
import argparse
import pandas as pd
from anndata import read_h5ad
import scimap as sm

def calculate_spatial_pscore(adata, cells):
    """
    Calculates spatial p-scores between pairs of cell types.
    """
    for i_index, i_item in enumerate(cells):
        for j_index in range(i_index + 1, len(cells)):
            j_item = cells[j_index]
            label = f'spatial_pscore_{i_item}__{j_item}'
            sm.tl.spatial_pscore(
                adata, 
                proximity=[i_item, j_item], 
                score_by='Image', 
                x_coordinate='x',
                y_coordinate='y', 
                phenotype='CT_final', 
                method='knn', 
                radius=20, 
                knn=3,
                imageid='Image', 
                subset=None, 
                label=label
            )
    
    return adata

def save_spatial_scores(adata, cells, output_dir):
    """
    Saves spatial p-score results to CSV files.
    """
    for i_index, i_item in enumerate(cells):
        for j_index in range(i_index + 1, len(cells)):
            j_item = cells[j_index]
            key = f'spatial_pscore_{i_item}__{j_item}'
            
            if key in adata.uns:
                output_file = os.path.join(output_dir, f"{i_item}_{j_item}.csv")
                adata.uns[key].to_csv(output_file)
            else:
                print(f"Warning: {key} not found in adata.uns. Skipping.")

def process_csv_files(csv_files, directory):
    """
    Processes CSV files to extract and merge Proximity Volume and Density columns.
    """
    selected_data_frames_v = []
    selected_data_frames_d = []

    for csv_file in sorted(csv_files):
        file_path = os.path.join(directory, csv_file)
        
        try:
            df = pd.read_csv(file_path, index_col='Image')
        except pd.errors.EmptyDataError:
            print(f"Warning: Empty CSV file at {file_path}")
            continue
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            continue
        
        # Rename columns to include the file name for uniqueness
        df.rename(
            columns={
                'Proximity Volume': f"{csv_file}_Proximity Volume", 
                'Proximity Density': f"{csv_file}_Proximity Density"
            }, 
            inplace=True
        )
        
        # Select relevant columns
        selected_columns_v = [col for col in df.columns if col.endswith('Proximity Volume')]
        selected_columns_d = [col for col in df.columns if col.endswith('Proximity Density')]

        if selected_columns_v:
            selected_data_frames_v.append(df[selected_columns_v])
        else:
            print(f"No 'Proximity Volume' columns selected in {csv_file}")

        if selected_columns_d:
            selected_data_frames_d.append(df[selected_columns_d])
        else:
            print(f"No 'Proximity Density' columns selected in {csv_file}")

    return selected_data_frames_v, selected_data_frames_d

def merge_and_save_dataframes(adata, selected_data_frames_v, selected_data_frames_d, output_dir, output_volume_name, output_density_name):
    """
    Merges and saves Proximity Volume and Density dataframes with metadata.
    """
    merged_df_volume = pd.concat(selected_data_frames_v, axis=1)
    merged_df_density = pd.concat(selected_data_frames_d, axis=1)
    
    # Merge with metadata
    metadata = adata.obs[['Image','Group']].drop_duplicates(subset=['Image','Group']).set_index('Image')

    concatenated_df_v = merged_df_volume.join(metadata, how='inner')
    concatenated_df_d = merged_df_density.join(metadata, how='inner')

    # Save the results to CSV files
    concatenated_df_v.to_csv(os.path.join(output_dir, output_volume_name))
    concatenated_df_d.to_csv(os.path.join(output_dir, output_density_name))

def main(args):
    # Load AnnData object
    adata = read_h5ad(args.anndata)

    # Calculate spatial p-scores
    cells = list(set(adata.obs['CT_final']))
    adata = calculate_spatial_pscore(adata, cells)

    # Create output directory if it doesn't exist
    output_dir = os.path.join(os.getcwd(), "output_pscore")
    os.makedirs(output_dir, exist_ok=True)

    # Save the modified AnnData object
    adata.write(filename=os.path.join(output_dir, args.output_andata_name))

    # Save spatial p-scores to CSV files
    save_spatial_scores(adata, cells, output_dir)

    # Process and merge CSV files
    file_list = [f for f in os.listdir(output_dir) if f.endswith('.csv')]
    selected_data_frames_v, selected_data_frames_d = process_csv_files(file_list, output_dir)

    # Merge and save the results
    merge_and_save_dataframes(adata, selected_data_frames_v, selected_data_frames_d, output_dir, args.output_volume_name, args.output_density_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial data and calculate p-scores.")
    parser.add_argument("--anndata", type=str, required=True, help="Path to the Anndata file.")
    parser.add_argument("--output_andata_name", type=str, default="output_adata.h5ad", help="Output file name for processed Anndata file.")
    parser.add_argument("--output_volume_name", type=str, default="merged_volume.csv", help="Output file name for merged volume CSV.")
    parser.add_argument("--output_density_name", type=str, default="merged_density.csv", help="Output file name for merged density CSV.")
    args = parser.parse_args()
    
    main(args)
