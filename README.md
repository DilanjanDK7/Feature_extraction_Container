# fMRI Feature Extraction Container

A comprehensive Docker-based solution for extracting analytical features from fMRI data using a Snakemake pipeline.

## Overview

This package provides a containerized workflow for calculating various analytical metrics from resting-state fMRI data, including:

- **ALFF** (Amplitude of Low-Frequency Fluctuation)
- **fALFF** (Fractional ALFF)
- **ReHo** (Regional Homogeneity)
- **Hurst Exponent**
- **Fractal Dimension**
- **Quantum Mechanical Fourier Transform (QM-FFT)**
- **Resting State Network (RSN) Analysis**

The pipeline uses Snakemake for workflow management and Docker for containerization, ensuring reproducibility and ease of deployment across different systems.

## Pipeline Workflow

The following diagram illustrates the basic workflow of the fMRI Feature Extraction Container:

```mermaid
graph TD
    A[Input: BIDS-formatted fMRI data] --> B[Docker Container]
    
    subgraph "Docker Container"
        B --> C[Snakemake Workflow]
        
        subgraph "Feature Extraction Pipeline"
            C --> D1[ALFF/fALFF]
            C --> D2[ReHo]
            C --> D3[Hurst Exponent]
            C --> D4[Fractal Dimension]
            C --> D5[QM-FFT]
            C --> D6[RSN Analysis]
        end
        
        D1 --> E[Output Generation]
        D2 --> E
        D3 --> E
        D4 --> E
        D5 --> E
        D6 --> E
    end
    
    E --> F[Output: Analytical Metrics]
```

Each feature extraction module operates independently, allowing you to select which analyses to run via the configuration file or command-line parameters.

## Documentation

Detailed documentation for using this pipeline is available in the following files:

- [**User Guide**](docs/USER_GUIDE.md): Comprehensive information about setting up and using the pipeline
- [**Quick Reference**](docs/QUICK_REFERENCE.md): Common commands and procedures for everyday use
- [**RSN Guide**](docs/RSN_GUIDE.md): Detailed information about the Resting State Network analysis
- [**Visualization Guide**](docs/VISUALIZATION_GUIDE.md): Instructions for visualizing and interpreting outputs
- [**Technical Details**](docs/TECHNICAL_DETAILS.md): Theoretical background and computational methods for core extracted features
- [**Advanced Technical Details**](docs/TECHNICAL_DETAILS_ADVANCED.md): Advanced features (RSN Analysis) and integrative approaches

For a complete documentation index, see the [Documentation README](docs/README.md).

## Requirements

- Docker (latest version recommended)
- At least 8GB RAM (16GB recommended for full analysis)
- Sufficient disk space for input data and results

## Quick Start

### Using the Simplified Pipeline Script

The easiest way to run the pipeline is using the simplified script:

```bash
./run_container_pipeline.sh --input /path/to/your/data --features alff,reho,qm_fft
```

#### Key Features of the Pipeline Script

- Automatically builds the Docker container if needed
- Places outputs in an "analytical_metrics" folder inside your input directory
- Allows selection of specific features to extract
- Supports custom parameter settings
- Can target specific subjects

#### Examples

Run all analyses for all subjects:
```bash
./run_container_pipeline.sh --input /path/to/your/data
```

Run ALFF and ReHo analyses with 4 cores:
```bash
./run_container_pipeline.sh --input /path/to/your/data --features alff,reho --cores 4
```

Run QM-FFT on a specific subject with custom parameters:
```bash
./run_container_pipeline.sh --input /path/to/your/data --subject sub-17017 --features qm_fft --param qm_fft_eps=1e-5
```

For more options:
```bash
./run_container_pipeline.sh --help
```

### Manual Approach

If you need more control, you can also manually build and run the container:

#### Building the Container

```bash
./run_container.sh build
```

For a clean build that doesn't use cache:

```bash
./run_container.sh build --no-cache
```

#### Running the Pipeline Manually

The basic command structure for running the pipeline is:

```bash
./run_container.sh run -v /path/to/input/data:/data/input -v /path/to/output/dir:/data/output -v $(pwd)/workflows:/app/workflows [command]
```

To run the entire Snakemake pipeline:

```bash
./run_container.sh run -v /path/to/input/data:/data/input -v /path/to/output/dir:/data/output -v $(pwd)/workflows:/app/workflows snakemake --snakefile /app/workflows/Snakefile -d /app/workflows --cores [N]
```

Replace `[N]` with the number of cores you want to use, or omit for maximum parallelization.

## Input Data Requirements

The pipeline expects input data organized according to BIDS format:

```
/path/to/input/data/
├── sub-<ID>/
│   └── func/
│       ├── sub-<ID>_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
│       └── sub-<ID>_task-rest_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz
```

## Output Structure

The pipeline generates output in the following structure:

```
/path/to/output/dir/
├── sub-<ID>/
│   └── func/
│       └── Analytical_metrics/
│           ├── ALFF/
│           │   ├── sub-<ID>_alff.nii.gz
│           │   └── sub-<ID>_falff.nii.gz
│           ├── ReHo/
│           │   └── sub-<ID>_reho.nii.gz
│           ├── Hurst/
│           │   └── sub-<ID>_hurst.nii.gz
│           ├── Fractal/
│           │   └── sub-<ID>_fractal.nii.gz
│           ├── QM_FFT/
│           │   └── sub-<ID>_qm_fft.h5
│           └── RSN/
│               └── sub-<ID>_rsn_activity.h5
```

## Configuration

The pipeline configuration is managed through the `workflows/config/config.yaml` file. Key parameters include:

- **Subject IDs**: Specify the subjects to process
- **I/O Paths**: Define input and output paths 
- **Processing Parameters**: Set specific parameters for each analysis
- **Core Usage**: Adjust parallel processing settings

Example configuration:

```yaml
# Input/Output settings
bids_derivatives_dir: "/data/input" 
output_dir: "/data/output"
overwrite_existing: true

# Subject settings
subjects:
  - "sub-17017"
  # - "sub-002"  # Uncomment to add more subjects

# Processing options
n_jobs: 4          # Number of parallel jobs to run
memory_limit: 8    # Memory limit in GB per job

# ALFF computation settings
compute_falff: true
alff_bandpass_low: 0.01
alff_bandpass_high: 0.08

# ReHo computation settings
reho_neighborhood: 27  # 27 for 3x3x3 cube

# ... additional settings
```

## Feature Details

### ALFF and fALFF

Measures the amplitude of BOLD signal fluctuations in the low-frequency range (typically 0.01-0.08 Hz), providing insight into regional spontaneous brain activity.

- **ALFF**: Amplitude of Low-Frequency Fluctuation directly measures the strength of low-frequency oscillations
- **fALFF**: Fractional ALFF represents the relative contribution of low-frequency oscillations to the entire frequency spectrum
- **mALFF**: Mean-normalized ALFF, providing improved sensitivity to neural signal
- **RSFA**: Resting State Fluctuation Amplitude, an additional metric for characterizing regional activity

Implementation leverages AFNI's `3dRSFC` for efficient and validated computation of all four metrics.

### ReHo (Regional Homogeneity)

Calculates similarity of the time series of a given voxel to those of its nearest neighbors, reflecting local synchronization of spontaneous brain activity.

- Kendall's coefficient of concordance (KCC) implementation
- Configurable neighborhood size (6, 18, or 26 neighbors)
- Optimized calculation for 3D volumes

### Hurst Exponent

Quantifies the long-range temporal dependence in the fMRI time series, indicating the predictability or persistence of the signal.

- **Two computation methods**:
  - Detrended Fluctuation Analysis (DFA): Robust to non-stationarity in the signal
  - Rescaled Range Analysis (R/S): Classical approach for long-memory processes
- Outputs both raw and normalized (z-score) maps
- Configurable parameters for scale range and trend removal
- Optimized implementation using the `nolds` library

### Fractal Dimension

Estimates the complexity and self-similarity of the fMRI time series, providing insight into the chaotic nature of brain activity.

- **Two computation methods**:
  - Higuchi Fractal Dimension (HFD): Direct time-domain estimate of signal complexity
  - Power Spectral Density (PSD) method: Frequency-domain approach using spectral slope
- Outputs both raw and normalized (z-score) maps
- Configurable parameters (kmax for Higuchi method, frequency range for PSD)
- Robust implementation with comprehensive error handling

### QM-FFT (Quantum Mechanical Fourier Transform)

Applies quantum mechanical principles to analyze fMRI data in the frequency domain, extracting several features:
- Magnitude
- Phase
- Temporal difference magnitude
- Temporal difference phase
- Local variance (optional)

### RSN (Resting State Network) Activity

Extracts time series from established resting-state networks using the Yeo atlas:
- 7-Network parcellation
- 17-Network parcellation

## Advanced Usage

### Running Individual Analyses

To run a specific analysis only (e.g., ALFF):

```bash
./run_container.sh run -v /path/to/input/data:/data/input -v $(pwd)/pipeline_outputs:/data/output -v $(pwd)/workflows:/app/workflows snakemake --snakefile /app/workflows/Snakefile -d /app/workflows --cores /data/output/sub-17017/func/Analytical_metrics/ALFF/sub-17017_alff.nii.gz
```

### Debugging and Testing

For testing with smaller sample sizes:

```bash
./run_container.sh run -v /path/to/input/data:/data/input -v $(pwd)/pipeline_outputs:/data/output -v $(pwd)/workflows:/app/workflows snakemake --snakefile /app/workflows/Snakefile -d /app/workflows --cores 1 /data/output/sub-17017/func/Analytical_metrics/QM_FFT/sub-17017_qm_fft.h5
```

Sample analysis uses reduced data sizes to speed up testing:
- QM-FFT uses spatial sampling (configurable in scripts/qm_fft_test.py)
- RSN analysis uses temporal sampling (configurable in config.yaml)

### Inspecting Outputs

To examine HDF5 outputs (QM-FFT, RSN analysis):

```bash
./run_container.sh run -v $(pwd)/pipeline_outputs:/data/output h5ls -r /data/output/sub-17017/func/Analytical_metrics/QM_FFT/sub-17017_qm_fft.h5
```

## Backing Up the Pipeline

Create a backup of the current configuration:

```bash
tar -czf feature_extraction_backup_$(date +%Y-%m-%d_%H%M%S).tar.gz workflows scripts docker run_container.sh environment.yml README.md requirements.txt
```

## Troubleshooting

- **Permission Issues**: If you encounter permission errors when accessing outputs, you may need to change ownership:
  ```bash
  sudo chown -R $(id -u):$(id -g) /path/to/output/dir
  ```

- **Docker Issues**: If the container fails to build, check Docker settings for sufficient memory allocation.

- **Snakemake Locks**: If Snakemake reports a locked directory, unlock it:
  ```bash
  ./run_container.sh run -v /path/to/input/data:/data/input -v $(pwd)/pipeline_outputs:/data/output -v $(pwd)/workflows:/app/workflows snakemake --snakefile /app/workflows/Snakefile -d /app/workflows --unlock
  ```


## Citation

If you use this pipeline in your research, please cite:

```
Author, Dilanjan DK (2025). fMRI Feature Extraction Container: A Comprehensive Pipeline for Analytical Metrics. ( Paper on the way)
```

## Copyright

Copyright (c) 2024 Dilanjan DK and BrainLab. All rights reserved.
For permission requests, contact: Dilanjan DK (ddiyabal@uwo.ca)

## Acknowledgments

This pipeline incorporates several open-source tools and packages:
- Snakemake for workflow management
- Nilearn and Nibabel for neuroimaging analysis
- AFNI for ReHo computation
- QM_FFT_Feature_Package for quantum mechanical analysis 

Developed at the BrainLab by Dilanjan DK ;-