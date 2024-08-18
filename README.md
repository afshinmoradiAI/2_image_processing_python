
# Image Processing Scripts Repository

This repository contains a collection of Python scripts and Jupyter notebooks dedicated to image processing and analysis, specifically focused on tasks such as quality control, cell typing, and spatial analysis. Each script is designed to handle different aspects of image processing and bioinformatics, particularly in the context of spatial transcriptomics.

## Repository Contents

### 1. **1_QC_normalization_image.py**
   - **Purpose:** This script performs quality control (QC) and normalization on image datasets. It includes functions to handle various preprocessing steps that are essential before conducting any further analysis on the images.
   - **Usage:**
     ```bash
     python 1_QC_normalization_image.py --input_dir /path/to/images --output_dir /path/to/output --normalize_method minmax
     ```
     - **Parameters:**
       - `--input_dir`: Path to the directory containing input images.
       - `--output_dir`: Path to save the normalized images.
       - `--normalize_method`: Method for normalization (`minmax`, `zscore`, etc.).

### 2. **2_Cell_typing.ipynb**
   - **Purpose:** This Jupyter notebook is dedicated to cell typing, leveraging machine learning techniques to classify cells based on extracted features. It's particularly suited for analyzing multiplex imaging data.
   - **Usage:**
     - Open the notebook in Jupyter:
       ```bash
       jupyter notebook 2_Cell_typing.ipynb
       ```
     - Follow the step-by-step cells within the notebook to:
       1. Load and preprocess your data.
       2. Extract features from the images.
       3. Train a classifier on the data.
       4. Evaluate the classifier's performance.
       5. Visualize the results.

### 3. **3_NB.py**
   - **Purpose:** This script performs normalization and background correction on image data, which is crucial for reducing noise and improving the accuracy of downstream analysis.
   - **Usage:**
     ```bash
     python 3_NB.py --input_file /path/to/input_image --output_file /path/to/output_image --method nb
     ```
     - **Parameters:**
       - `--input_file`: Path to the input image file.
       - `--output_file`: Path to save the processed image.
       - `--method`: Specify the normalization/background correction method.

### 4. **4_CellEnrich_InNB.py**
   - **Purpose:** This script is used for enrichment analysis, likely in the context of cell typing or classification results, refining the analysis by focusing on specific cell populations.
   - **Usage:**
     ```bash
     python 4_CellEnrich_InNB.py --input_data /path/to/classification_results --output_data /path/to/enrichment_results
     ```
     - **Parameters:**
       - `--input_data`: Path to the file containing classification results.
       - `--output_data`: Path to save the enrichment analysis results.

### 5. **5_Dist_intraction_scimap.py**
   - **Purpose:** This script calculates distance interactions between various cellular components within the image, leveraging the `scimap` library. It is useful for spatial analysis of cellular interactions within the tissue microenvironment.
   - **Usage:**
     ```bash
     python 5_Dist_intraction_scimap.py --image_data /path/to/image_data --output_dir /path/to/output --distance_threshold 50
     ```
     - **Parameters:**
       - `--image_data`: Path to the input image data file.
       - `--output_dir`: Directory to save the interaction analysis results.
       - `--distance_threshold`: Threshold distance for considering interactions (in pixels or micrometers, depending on your image scale).

### 6. **6_Pscore_scimap.py**
   - **Purpose:** This script computes P-scores using the `scimap` library, which is a measure used in spatial transcriptomics to quantify spatial relationships within the tissue.
   - **Usage:**
     ```bash
     python 6_Pscore_scimap.py --input_file /path/to/input_data --output_file /path/to/output_pscore
     ```
     - **Parameters:**
       - `--input_file`: Path to the input file containing spatial data.
       - `--output_file`: Path to save the P-score results.

## How to Use

1. **Installation:**
   - Ensure you have Python installed.
   - Install the required dependencies by running:
     ```bash
     pip install -r requirements.txt
     ```

2. **Running Scripts:**
   - Each script can be run individually from the command line with the specified parameters.

3. **Jupyter Notebooks:**
   - The `2_Cell_typing.ipynb` notebook can be opened and run in a Jupyter environment:
     ```bash
     jupyter notebook 2_Cell_typing.ipynb
     ```

4. **Dependencies:**
   - The necessary Python libraries for each script are specified within the script or notebook. Make sure to install these libraries before running the scripts.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please submit pull requests or raise issues to discuss any improvements or features you’d like to see.
