# Mangroves_change_ML
Mangroves change detection using GEE 
# Multi-Temporal Mangrove Dynamics & GeoAI Frameworks: 

An advanced Earth Observation (EO) research repository demonstrating the integration of traditional remote sensing feature engineering with modern deep learning optimization pipelines for fragile coastal ecosystem monitoring.

## 🛰️ Project Methodology & Architecture

### 1. Cloud-Based Remote Sensing Engineering (GEE API)
The JavaScript engine (`gee_feature_engineering.js`) executes scalable cloud-computing algorithms across the Sentinel-2 Harmonized Surface Reflectance archive over Bhitarkanika, Odisha:
* **Atmospheric & Cloud Suppression:** Bitwise QA60 masking to isolate clear-sky targets.
* **High-Dimensional Feature Stacking:** Extraction of >35 custom spectral indices spanning red-edge vegetative health (TCARI, CCCI, CIRE), moisture variations (MNDWI, NDMI), and canopy structure (EVI, SAVI).
* **Multi-Criteria Topographic Masking:** Integration of SRTM Digital Elevation Models (DEM) to clean and isolate tidal-zone vegetative vectors (<80m elevation threshold).

### 2. Downstream Deep Learning Execution (PyTorch U-Net Core)
The Python module (`pytorch_segmentation_unet.py`) transitions the GEE feature maps into an advanced AI operational space:
* **Spatial Ingestion Engine:** Natively ingests multi-band GeoTIFF targets via `rasterio` and `GeoPandas`, breaking complex arrays into non-overlapping spatial patches for tensor operations.
* **Semantic Segmentation Head:** Implements a Deep Convolutional Architecture designed for high-resolution pixel-wise classification, bridging classical machine learning (Random Forest, SVM) with modern GeoAI.
