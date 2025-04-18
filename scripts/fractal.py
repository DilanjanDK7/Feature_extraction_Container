#!/usr/bin/env python3
"""
Script to compute Fractal Dimension from fMRI data.
"""

import os
import argparse
import time
import numpy as np
import nibabel as nib
from scipy import signal
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed
import sys


def compute_higuchi_fd(ts, kmax=10):
    """
    Compute Higuchi Fractal Dimension of a time series.
    
    Parameters:
    -----------
    ts : numpy.ndarray
        Input time series
    kmax : int
        Maximum delay/lag. Default is 10.
        
    Returns:
    --------
    float
        Higuchi Fractal Dimension
    """
    n = len(ts)
    if n < kmax * 2:
        # Time series too short for reliable estimation with this kmax
        return np.nan
        
    lk = np.zeros(kmax)
    x_reg = np.array(range(kmax))
    y_reg = np.zeros(kmax)
    
    for k in range(1, kmax + 1):
        lm = np.zeros(k)
        
        for m in range(k):
            # Construct subsequence
            ll = 0
            
            # Number of subsequences
            n_m = int((n - m) / k)
            
            # Skip if too few points
            if n_m <= 1:
                lm[m] = np.nan
                continue
                
            for i in range(1, n_m):
                ll += abs(ts[m + i * k] - ts[m + (i - 1) * k])
            
            ll /= k  # Normalize with factor k
            # Length for subsequence starting at m
            lm[m] = ll * (n - 1) / (n_m * k)
        
        # Mean length for step k
        lk[k - 1] = np.nanmean(lm)
        y_reg[k - 1] = np.log(lk[k - 1]) if lk[k - 1] > 0 else np.nan
    
    # Filter out invalid values
    valid = ~np.isnan(y_reg)
    if np.sum(valid) < 2:
        return np.nan
        
    x_reg_valid = np.log(1.0 / np.array(range(1, kmax + 1)))[valid]
    y_reg_valid = y_reg[valid]
    
    # Perform linear regression
    try:
        slopes = np.polyfit(x_reg_valid, y_reg_valid, 1)
        
        # Return the slope (fractal dimension)
        return slopes[0]
    except:
        return np.nan


def compute_psd_fd(ts):
    """
    Compute Fractal Dimension using Power Spectral Density (PSD) method.
    
    Parameters:
    -----------
    ts : numpy.ndarray
        Input time series
        
    Returns:
    --------
    float
        Fractal Dimension from PSD slope
    """
    # Remove mean
    ts = ts - np.mean(ts)
    
    # Compute PSD using Welch's method
    try:
        freqs, psd = signal.welch(ts, nperseg=min(256, len(ts)//4))
        
        # Avoid zero frequency and use only positive frequencies
        mask = (freqs > 0)
        freqs = freqs[mask]
        psd = psd[mask]
        
        # Perform linear regression in log-log space
        if len(freqs) > 5:  # Ensure enough points for regression
            # Remove zeros from PSD to avoid log(0)
            valid = psd > 0
            if np.sum(valid) < 5:
                return np.nan
                
            log_freqs = np.log10(freqs[valid])
            log_psd = np.log10(psd[valid])
            
            slopes = np.polyfit(log_freqs, log_psd, 1)
            
            # Calculate fractal dimension from PSD slope
            # FD = (5 - beta) / 2, where beta is the negative of the slope
            beta = -slopes[0]
            fd = (5 - beta) / 2
            
            # Valid range for fractal dimension is typically 1-2
            if 0.5 < fd < 2.5:  # Slightly wider range to account for estimation errors
                return fd
    except:
        pass
        
    return np.nan


def compute_voxel_fractal(ts, method='higuchi', kmax=10):
    """
    Compute Fractal Dimension for a single voxel time series.
    
    Parameters:
    -----------
    ts : numpy.ndarray
        Input time series for a single voxel
    method : str
        Method for fractal dimension calculation ('higuchi' or 'psd')
    kmax : int
        Maximum lag parameter for Higuchi method
        
    Returns:
    --------
    float
        Fractal Dimension value or NaN if calculation fails
    """
    # Skip if time series has no variance
    if np.std(ts) <= 1e-6:
        return np.nan
    
    try:
        # Calculate fractal dimension using specified method
        if method == 'higuchi':
            return compute_higuchi_fd(ts, kmax)
        elif method == 'psd':
            return compute_psd_fd(ts)
        else:
            return compute_higuchi_fd(ts, kmax)  # Default to Higuchi
    except Exception:
        return np.nan


def compute_fractal(fmri_file, output_file, method='higuchi', kmax=10, mask_file=None, n_jobs=1, min_var=1e-6):
    """
    Compute Fractal Dimension from fMRI data.
    
    Parameters:
    -----------
    fmri_file : str
        Path to the input fMRI data
    output_file : str
        Path to save the Fractal Dimension output
    method : str, optional
        Method for fractal dimension calculation ('higuchi' or 'psd'). Default is 'higuchi'.
    kmax : int, optional
        Maximum lag parameter for Higuchi method. Default is 10.
    mask_file : str, optional
        Path to a brain mask. If not provided, a mask will be created based on signal variance.
    n_jobs : int, optional
        Number of parallel jobs to run. Default is 1.
    min_var : float, optional
        Minimum variance threshold for time series. Default is 1e-6.
    """
    print("Starting Fractal Dimension calculation...")
    start_time = time.time()
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load fMRI data
    try:
        print(f"Loading fMRI data from {fmri_file}...")
        img = nib.load(fmri_file)
        data = img.get_fdata()
        affine = img.affine
        header = img.header
    except Exception as e:
        print(f"Error loading fMRI data: {e}")
        sys.exit(1)
    
    # Get dimensions
    nx, ny, nz, nt = data.shape
    print(f"Data dimensions: {nx} x {ny} x {nz} x {nt}")
    
    # Create or load mask
    try:
        if mask_file:
            print(f"Loading mask from {mask_file}...")
            mask_img = nib.load(mask_file)
            mask = mask_img.get_fdata() > 0
        else:
            print("Creating mask from fMRI data...")
            # Simple mask: voxels with non-zero variance
            mask = np.std(data, axis=3) > min_var
    except Exception as e:
        print(f"Error creating/loading mask: {e}")
        sys.exit(1)
    
    # Ensure mask dimensions match
    if mask.shape[:3] != data.shape[:3]:
        print(f"Mask dimensions {mask.shape[:3]} do not match fMRI data dimensions {data.shape[:3]}")
        sys.exit(1)
    
    # Create coordinates for voxels to process
    coords = [(x, y, z) for x in range(nx) for y in range(ny) for z in range(nz) if mask[x, y, z]]
    print(f"Processing {len(coords)} voxels...")
    
    # Parallel processing function
    def process_voxel(coord):
        x, y, z = coord
        ts = data[x, y, z, :]
        return (x, y, z, compute_voxel_fractal(ts, method, kmax))
    
    # Set number of jobs
    actual_n_jobs = min(n_jobs, multiprocessing.cpu_count())
    print(f"Using {actual_n_jobs} parallel jobs")
    
    # Process voxels in parallel with progress bar
    results = Parallel(n_jobs=actual_n_jobs)(
        delayed(process_voxel)(coord) for coord in tqdm(coords, desc=f"Computing {method.capitalize()} FD")
    )
    
    # Create Fractal Dimension map
    fd_map = np.zeros((nx, ny, nz))
    
    # Fill the map with results
    for x, y, z, fd in results:
        fd_map[x, y, z] = fd
    
    # Replace NaNs with 0
    fd_map = np.nan_to_num(fd_map)
    
    # Normalize fractal dimension map for visualization (z-score)
    # Valid fractal dimension range for most physiological signals
    brain_mask = (fd_map > 1.0) & (fd_map < 2.0)
    if np.sum(brain_mask) > 0:
        mean_fd = np.mean(fd_map[brain_mask])
        std_fd = np.std(fd_map[brain_mask])
        print(f"Fractal dimension statistics: mean = {mean_fd:.4f}, std = {std_fd:.4f}")
        
        if std_fd > 0:
            fd_map_norm = np.zeros_like(fd_map)
            fd_map_norm[brain_mask] = (fd_map[brain_mask] - mean_fd) / std_fd
            
            # Save both raw and normalized maps
            print(f"Saving Fractal Dimension map to {output_file}...")
            fd_img = nib.Nifti1Image(fd_map, affine, header)
            nib.save(fd_img, output_file)
            
            norm_output = output_file.replace('.nii.gz', '_norm.nii.gz')
            print(f"Saving normalized Fractal Dimension map to {norm_output}...")
            fd_norm_img = nib.Nifti1Image(fd_map_norm, affine, header)
            nib.save(fd_norm_img, norm_output)
        else:
            print("Warning: Standard deviation is zero, cannot normalize.")
            fd_img = nib.Nifti1Image(fd_map, affine, header)
            nib.save(fd_img, output_file)
    else:
        print("Warning: No valid fractal dimensions calculated.")
        fd_img = nib.Nifti1Image(fd_map, affine, header)
        nib.save(fd_img, output_file)
    
    elapsed_time = time.time() - start_time
    print(f"Fractal Dimension calculation completed in {elapsed_time:.2f} seconds")
    
    return output_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute Fractal Dimension from fMRI data')
    parser.add_argument('--fmri', required=True, help='Input fMRI file')
    parser.add_argument('--output', required=True, help='Output Fractal Dimension file')
    parser.add_argument('--method', choices=['higuchi', 'psd'], default='higuchi',
                      help='Method for fractal dimension calculation (higuchi or psd). Default is higuchi.')
    parser.add_argument('--kmax', type=int, default=10, help='Maximum lag for Higuchi method. Default is 10.')
    parser.add_argument('--mask', help='Brain mask (optional, will be created if not provided)')
    parser.add_argument('--n-jobs', type=int, default=1, 
                      help='Number of parallel jobs for computation. Default is 1.')
    parser.add_argument('--min-var', type=float, default=1e-6,
                      help='Minimum variance threshold for time series. Default is 1e-6.')
    
    args = parser.parse_args()
    
    compute_fractal(args.fmri, args.output, args.method, args.kmax, args.mask, args.n_jobs, args.min_var) 