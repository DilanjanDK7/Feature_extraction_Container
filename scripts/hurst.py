#!/usr/bin/env python3
"""
Script to compute Hurst exponent from fMRI data.
"""

import os
import argparse
import time
import numpy as np
import nibabel as nib
import nolds  # For Hurst exponent calculation
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed
import sys

def compute_voxel_hurst(ts, method='dfa'):
    """
    Compute Hurst exponent for a single voxel time series.
    
    Parameters:
    -----------
    ts : numpy.ndarray
        Input time series for a single voxel
    method : str
        Method for Hurst exponent calculation ('dfa' or 'rs')
        
    Returns:
    --------
    float
        Hurst exponent value or NaN if calculation fails
    """
    # Skip if time series has no variance
    if np.std(ts) <= 1e-6:
        return np.nan
    
    try:
        # Set Hurst calculation method
        if method == 'dfa':
            return nolds.dfa(ts)
        elif method == 'rs':
            return nolds.hurst_rs(ts)
        else:
            return nolds.dfa(ts)  # Default to DFA
    except Exception:
        return np.nan

def compute_hurst(fmri_file, output_file, method='dfa', mask_file=None, n_jobs=1, min_var=1e-6):
    """
    Compute Hurst exponent from fMRI data.
    
    Parameters:
    -----------
    fmri_file : str
        Path to the input fMRI data
    output_file : str
        Path to save the Hurst exponent output
    method : str, optional
        Method for Hurst exponent calculation ('dfa' or 'rs'). Default is 'dfa'.
    mask_file : str, optional
        Path to a brain mask. If not provided, a mask will be created based on signal variance.
    n_jobs : int, optional
        Number of parallel jobs to run. Default is 1.
    min_var : float, optional
        Minimum variance threshold. Time series with variance below this will be skipped.
    """
    print("Starting Hurst exponent calculation...")
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
        return (x, y, z, compute_voxel_hurst(ts, method))
    
    # Set number of jobs
    actual_n_jobs = min(n_jobs, multiprocessing.cpu_count())
    print(f"Using {actual_n_jobs} parallel jobs")
    
    # Process voxels in parallel with progress bar
    results = Parallel(n_jobs=actual_n_jobs)(
        delayed(process_voxel)(coord) for coord in tqdm(coords, desc="Computing Hurst")
    )
    
    # Create Hurst exponent map
    hurst_map = np.zeros((nx, ny, nz))
    
    # Fill the map with results
    for x, y, z, h in results:
        hurst_map[x, y, z] = h
    
    # Replace NaNs with 0
    hurst_map = np.nan_to_num(hurst_map)
    
    # Normalize Hurst map for visualization (z-score)
    brain_mask = (hurst_map > 0) & (hurst_map < 2)  # Valid Hurst range
    if np.sum(brain_mask) > 0:
        mean_hurst = np.mean(hurst_map[brain_mask])
        std_hurst = np.std(hurst_map[brain_mask])
        print(f"Hurst exponent statistics: mean = {mean_hurst:.4f}, std = {std_hurst:.4f}")
        
        if std_hurst > 0:
            hurst_map_norm = np.zeros_like(hurst_map)
            hurst_map_norm[brain_mask] = (hurst_map[brain_mask] - mean_hurst) / std_hurst
            
            # Save both raw and normalized maps
            print(f"Saving Hurst exponent map to {output_file}...")
            hurst_img = nib.Nifti1Image(hurst_map, affine, header)
            nib.save(hurst_img, output_file)
            
            norm_output = output_file.replace('.nii.gz', '_norm.nii.gz')
            print(f"Saving normalized Hurst exponent map to {norm_output}...")
            hurst_norm_img = nib.Nifti1Image(hurst_map_norm, affine, header)
            nib.save(hurst_norm_img, norm_output)
        else:
            print("Warning: Standard deviation is zero, cannot normalize.")
            hurst_img = nib.Nifti1Image(hurst_map, affine, header)
            nib.save(hurst_img, output_file)
    else:
        print("Warning: No valid Hurst exponents calculated.")
        hurst_img = nib.Nifti1Image(hurst_map, affine, header)
        nib.save(hurst_img, output_file)
    
    elapsed_time = time.time() - start_time
    print(f"Hurst exponent calculation completed in {elapsed_time:.2f} seconds")
    
    return output_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute Hurst exponent from fMRI data')
    parser.add_argument('--fmri', required=True, help='Input fMRI file')
    parser.add_argument('--output', required=True, help='Output Hurst exponent file')
    parser.add_argument('--method', choices=['dfa', 'rs'], default='dfa',
                      help='Method for Hurst calculation (dfa or rs). Default is dfa.')
    parser.add_argument('--mask', help='Brain mask (optional, will be created if not provided)')
    parser.add_argument('--n-jobs', type=int, default=1, 
                      help='Number of parallel jobs for computation. Default is 1.')
    parser.add_argument('--min-var', type=float, default=1e-6,
                      help='Minimum variance threshold for time series. Default is 1e-6.')
    
    args = parser.parse_args()
    
    compute_hurst(args.fmri, args.output, args.method, args.mask, args.n_jobs, args.min_var) 