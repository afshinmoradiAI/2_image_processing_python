# **Image Processing**

This repository contains a Python script for quality control and normalization of Fusion images. The script is designed to preprocess images and ensure that they meet the required standards before further analysis.

**Table of Contents**

- Introduction
- Features
- Installation
- Usage
- Requirements
- Contributing
- License

**Introduction**

The 1_QC_normalization_image.py script provides functionality for performing quality control and normalization on a set of images. This is crucial in ensuring that all images are standardized and suitable for subsequent processing steps, such as feature extraction or model training.
Features

Quality Control: Checks the images for quality issues and flags those that do not meet specified criteria.
Normalization: Normalizes image data to ensure consistency across the dataset.
Logging: Provides detailed logging of the processing steps for traceability.


**Installation**
```
git clone https://github.com/your-username/image-processing.git
cd image-processing
pip install -r requirements.txt
```


**Usage**
Once the dependencies are installed, you can run the script by executing the following command:
```
python 1_QC_normalization_image.py --input data/input_data.csv --output results/processed_data.h5ad --exclude_core Core1,Core2 --exclude_marker Marker1,Marker2

```
Command Line Arguments
    --input_path: The directory containing the images to be processed.
    --output_path: The directory where the processed images will be saved.



**Requirements**
    Python 3.8+
    Required Python packages (listed in requirements.txt)



To install the necessary packages, use:
```
pip install -r requirements.txt
```
