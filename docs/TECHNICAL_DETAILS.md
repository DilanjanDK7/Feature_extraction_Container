# Technical Details: fMRI Analytical Features

This document provides comprehensive information about the analytical features extracted by the fMRI Feature Extraction Container pipeline, including theoretical background, computational methods, and interpretation guidelines.

## Table of Contents

1. [ALFF (Amplitude of Low-Frequency Fluctuation)](#alff)
2. [fALFF (Fractional ALFF)](#falff)
3. [ReHo (Regional Homogeneity)](#reho)
4. [Hurst Exponent](#hurst)
5. [Fractal Dimension](#fractal)
6. [Quantum Mechanical Fourier Transform (QM-FFT)](#qm-fft)
7. [Resting State Network (RSN) Analysis](#rsn)

<a name="alff"></a>
## 1. ALFF (Amplitude of Low-Frequency Fluctuation)

### Theoretical Background

ALFF quantifies the amplitude of spontaneous low-frequency oscillations (typically 0.01-0.08 Hz) in the BOLD signal, providing a measure of regional spontaneous brain activity. It was first proposed by Zang et al. (2007) as a metric for characterizing the local intensity of spontaneous fluctuations in the resting brain.

### Calculation Method

The ALFF computation follows these steps:
1. Convert the time-domain BOLD signal to the frequency domain using Fast Fourier Transform (FFT)
2. Extract the square root of power spectrum within the low-frequency range (0.01-0.08 Hz)
3. Calculate the sum of amplitudes within this frequency band

Mathematically, for a time series x(t), ALFF is defined as:

```
ALFF = √(∑(|FFT(x(t))|²)) for frequencies in the 0.01-0.08 Hz range
```

### Interpretation

Higher ALFF values suggest stronger spontaneous activity in a brain region. ALFF is sensitive to physiological noise (particularly in blood vessels, ventricles, and cisterns), making it useful for examining regional brain activity but requiring careful interpretation in specific regions.

ALFF has been widely used to detect abnormal brain activity in various disorders, including:
- Alzheimer's disease
- Depression
- Schizophrenia
- Attention deficit hyperactivity disorder

### Configuration Options

In this pipeline, ALFF calculation can be customized with the following parameters:
- `alff_bandpass_low`: Lower frequency bound (default: 0.01 Hz)
- `alff_bandpass_high`: Upper frequency bound (default: 0.08 Hz)
- `detrend_method`: Method for detrending the signal (options: "linear", "constant", or "none")
- `normalize_method`: Method for normalizing ALFF values (options: "zscore", "percent", or "none")

<a name="falff"></a>
## 2. fALFF (Fractional ALFF)

### Theoretical Background

fALFF, introduced by Zou et al. (2008), addresses some limitations of ALFF by measuring the ratio of power in the low-frequency range relative to the entire frequency range. This normalization reduces sensitivity to physiological noise.

### Calculation Method

The fALFF computation follows these steps:
1. Calculate the ALFF in the low-frequency range (0.01-0.08 Hz)
2. Calculate the amplitude sum across the entire frequency range (0-Nyquist frequency)
3. Compute the ratio of low-frequency amplitude to the entire frequency range amplitude

Mathematically:

```
fALFF = (∑|FFT(x(t))| for frequencies in 0.01-0.08 Hz) / (∑|FFT(x(t))| for all frequencies)
```

### Interpretation

fALFF represents the relative contribution of low-frequency oscillations to the entire detectable frequency range. It has greater specificity for detecting spontaneous brain activity than ALFF and better suppresses non-specific signal components.

Higher fALFF values indicate that low-frequency oscillations contribute more to the total frequency spectrum, often interpreted as more organized neural activity.

### Configuration Options

fALFF shares configuration options with ALFF but adds one unique parameter:
- `compute_falff`: Boolean to enable/disable fALFF computation (default: true)

<a name="reho"></a>
## 3. ReHo (Regional Homogeneity)

### Theoretical Background

Regional Homogeneity (ReHo), introduced by Zang et al. (2004), measures the similarity or synchronization of the time series of a given voxel to those of its nearest neighbors. It evaluates local synchronization of spontaneous brain activity, based on the assumption that brain activity occurs in clusters rather than in isolated voxels.

### Calculation Method

ReHo is calculated using Kendall's coefficient of concordance (KCC) to measure the similarity among time series of neighboring voxels:

1. For each voxel, identify its K nearest neighboring voxels (typically 26, forming a 3×3×3 cube)
2. Calculate Kendall's coefficient of concordance (KCC) for the time series of these K+1 voxels

The KCC calculation is defined as:

```
W = (∑(Ri)² - n(R̄)²) / (1/12 × K² × (n³-n))
```

Where:
- W is the KCC value
- Ri is the sum rank of the ith time point
- n is the number of time points
- K is the number of voxels included in the calculation
- R̄ is the mean of the Ri

### Interpretation

Higher ReHo values indicate greater local synchronization of spontaneous brain activity, suggesting more coordinated function among neighboring neurons. 

ReHo can identify changes in local connectivity and has been applied to various neurological and psychiatric conditions, including:
- Parkinson's disease
- Depression
- Autism spectrum disorders
- Epilepsy

### Configuration Options

In this pipeline, ReHo calculation can be customized with:
- `reho_neighborhood`: Number of neighbors to include (options: 27 for 3×3×3 cube, 19 for face+edge, 7 for face only)

<a name="hurst"></a>
## 4. Hurst Exponent

### Theoretical Background

The Hurst exponent (H), named after Harold Edwin Hurst, quantifies the long-range temporal dependence or autocorrelation in the fMRI time series. It measures the predictability or persistence of the signal over time.

### Calculation Method

This pipeline implements two main methods for calculating the Hurst exponent:

1. **Detrended Fluctuation Analysis (DFA)**:
   - Integrates the time series to convert it to a random walk
   - Divides the integrated time series into non-overlapping segments
   - Detrends each segment and calculates the fluctuation
   - Determines the scaling relationship between fluctuation and segment size

2. **Rescaled Range Analysis (R/S)**:
   - Divides the time series into segments
   - For each segment, calculates the range R (max - min) and standard deviation S
   - Plots log(R/S) against log(segment length)
   - Estimates H as the slope of this relationship

The Hurst exponent (H) typically falls in the range of 0 to 1:

### Interpretation

- H ≈ 0.5: Random walk (Brownian motion) - no long-range correlations
- 0 < H < 0.5: Anti-persistent behavior - increases tend to be followed by decreases
- 0.5 < H < 1: Persistent behavior - increases tend to be followed by increases

In the context of fMRI:
- Higher H values (>0.5) indicate more predictable, structured brain activity
- Lower H values (<0.5) indicate more irregular, anti-correlated activity
- H ≈ 0.5 suggests random, unpredictable activity

The Hurst exponent has been used to characterize alterations in brain dynamics in conditions like Alzheimer's disease, schizophrenia, and epilepsy.

### Configuration Options

In this pipeline, Hurst exponent calculation can be customized with:
- `compute_hurst`: Boolean to enable/disable calculation (default: true)
- `hurst_method`: Method for calculation (options: "dfa" or "rs")

<a name="fractal"></a>
## 5. Fractal Dimension

### Theoretical Background

Fractal dimension (FD) is a measure of the complexity and self-similarity of a time series. In fMRI analysis, it quantifies the irregularity or complexity of BOLD signal fluctuations, providing insight into the chaotic or self-similar nature of brain activity.

### Calculation Method

The pipeline implements two main methods for calculating fractal dimension:

1. **Higuchi Fractal Dimension (HFD)**:
   - Creates new time series by subsampling the original time series with different intervals
   - Calculates the "length" of each new time series
   - Determines the scaling relationship between the calculated lengths and the sampling intervals
   - Estimates FD as the slope of this relationship

2. **Power Spectral Density (PSD) method**:
   - Computes the power spectral density of the time series
   - Plots log(power) against log(frequency)
   - Estimates the fractal dimension based on the slope β of this plot using the relation FD = (5-β)/2

### Interpretation

Fractal dimension typically ranges from 1 to 2 for time series data:
- FD closer to 1: Smoother, more predictable signal
- FD closer to 2: More complex, irregular signal

In the context of fMRI:
- Higher FD values indicate more complex, random-like brain activity
- Lower FD values suggest more regular, predictable brain dynamics

Fractal dimension provides complementary information to the Hurst exponent and has been applied to study various brain states and conditions, including:
- Sleep states
- Aging
- Neurodegenerative disorders
- Consciousness levels

### Configuration Options

In this pipeline, fractal dimension calculation can be customized with:
- `compute_fractal`: Boolean to enable/disable calculation (default: true)
- `fractal_method`: Method for calculation (options: "higuchi" or "psd")

<a name="qm-fft"></a>
## 6. Quantum Mechanical Fourier Transform (QM-FFT)

### Theoretical Background

The Quantum Mechanical Fourier Transform (QM-FFT) is an advanced analytical method that applies principles from quantum mechanics to analyze fMRI data in the frequency domain. It interprets brain activity as a quantum mechanical system, providing unique insights into the spatiotemporal dynamics of brain function.

### Calculation Method

The QM-FFT analysis in this pipeline proceeds through several key steps:

1. **Spatial Normalization**: Normalizes voxel coordinates to a quantum-appropriate reference frame
2. **Forward FFT**: Transforms the time series data to k-space (frequency domain) using a non-uniform FFT
3. **Spherical Masking**: Applies spherical masks in k-space centered at specific frequency ranges
4. **Inverse Transformation**: Transforms masked data back to real space
5. **Feature Extraction**: Computes multiple features from the transformed data:
   - Magnitude: Strength of the signal components
   - Phase: Temporal relationships between signals
   - Temporal difference magnitude: Changes in signal strength over time
   - Temporal difference phase: Changes in signal phase over time
   - Local variance (optional): Spatial heterogeneity of the signal

The process results in a rich set of features that characterize both the spatial and temporal aspects of brain activity.

### Interpretation

QM-FFT features provide a multidimensional view of brain activity:

- **Magnitude maps**: Indicate the strength of specific frequency components, highlighting regions with prominent oscillatory activity
- **Phase maps**: Reveal temporal relationships between brain regions, potentially indicating communication or synchronization
- **Temporal difference maps**: Show how brain activity patterns change over time, capturing dynamic aspects that static measures miss
- **Local variance**: Quantifies the spatial heterogeneity of brain activity, potentially indicating functional specialization

These features are particularly valuable for capturing complex, non-linear aspects of brain function that traditional methods might miss.

### Configuration Options

In this pipeline, QM-FFT calculation can be customized with:
- `compute_qm_fft`: Boolean to enable/disable calculation (default: true)
- `qm_fft_eps`: Epsilon value for numerical stability (default: 1e-6)
- `qm_fft_radius`: Radius parameter for spherical masking (default: 0.6)
- `qm_fft_local_k`: Local k parameter for neighborhood analysis (default: 5)

<a name="rsn"></a>
## 7. Resting State Network (RSN) Analysis

### Theoretical Background

Resting State Networks (RSNs) are coherent patterns of functionally connected brain regions that consistently show synchronized activity during rest. These networks persist even in the absence of task performance and reflect fundamental organizational principles of the brain.

The pipeline leverages the established Yeo atlas, which defines:
- A 7-network parcellation: Dividing the cortex into seven major functional networks
- A 17-network parcellation: Providing a more fine-grained division of functional networks

### Calculation Method

The RSN analysis in this pipeline follows these steps:

1. **Atlas Registration**: Align the Yeo atlas (7 and 17 network parcellations) to the subject's fMRI space
2. **Time Series Extraction**: For each network in the atlas, extract the average time series across all voxels within that network
3. **Network Activity Calculation**: Compute summary statistics for each network's time series
4. **Connectivity Analysis**: Calculate correlations between different networks' time series to assess functional connectivity

### Interpretation

RSN analysis provides insights into:

- **Network Integrity**: How well each canonical network maintains its functional coherence
- **Between-Network Interactions**: How different functional systems communicate with each other
- **Network Dynamics**: How network activity and connectivity fluctuate over time

Alterations in RSN patterns have been associated with various neurological and psychiatric conditions, including:
- Alzheimer's disease (disrupted default mode network)
- Schizophrenia (altered connectivity between networks)
- Depression (changes in salience and executive control networks)
- Autism (atypical connectivity patterns)

### Configuration Options

In this pipeline, RSN analysis can be customized with:
- `compute_rsn`: Boolean to enable/disable RSN extraction (default: true)
- `rsn_sample_tp`: Number of time points to sample for testing mode (default: 100)

## References

1. Zang YF, He Y, Zhu CZ, et al. (2007). Altered baseline brain activity in children with ADHD revealed by resting-state functional MRI. Brain Dev, 29(2):83-91.

2. Zou QH, Zhu CZ, Yang Y, et al. (2008). An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: fractional ALFF. J Neurosci Methods, 172(1):137-141.

3. Zang Y, Jiang T, Lu Y, et al. (2004). Regional homogeneity approach to fMRI data analysis. Neuroimage, 22(1):394-400.

4. Mandelbrot BB, Van Ness JW. (1968). Fractional Brownian motions, fractional noises and applications. SIAM Review, 10(4):422-437.

5. Higuchi T. (1988). Approach to an irregular time series on the basis of the fractal theory. Physica D, 31(2):277-283.

6. Yeo BT, Krienen FM, Sepulcre J, et al. (2011). The organization of the human cerebral cortex estimated by intrinsic functional connectivity. J Neurophysiol, 106(3):1125-1165. 