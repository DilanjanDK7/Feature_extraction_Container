#!/usr/bin/env python3
import argparse
import os
import tempfile
import subprocess
import nibabel as nib
import shutil
import sys

def run_3dRSFC_all(fmri, mask, band_low, band_high, metrics, prefix):
    """
    Run AFNI's 3dRSFC to compute specified metrics.
    metrics: list of strings from ['ALFF','mALFF','fALFF','RSFA'].
    """
    cmd = [
        "3dRSFC",
        "-input", fmri,
        "-band", str(band_low), str(band_high),
        "-prefix", prefix
    ]
    if mask:
        cmd += ["-mask", mask]
    for m in metrics:
        cmd += ["-"+m]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return prefix
    except subprocess.CalledProcessError as e:
        print(f"Error running 3dRSFC: {e}")
        print(f"Command: {' '.join(cmd)}")
        print(f"Stderr: {e.stderr.decode()}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Compute resting‐state metrics via AFNI 3dRSFC'
    )
    parser.add_argument('--fmri',    required=True,
                        help='Path to 4D BOLD NIfTI')
    parser.add_argument('--mask',    default=None,
                        help='Brain mask NIfTI (optional)')
    parser.add_argument('--output',  required=True,
                        help='Output directory for metric NIfTIs')
    parser.add_argument('--tr',      type=float, default=2.0,
                        help='TR in seconds (unused by AFNI)')
    # Metrics flags (all on by default)
    parser.add_argument('--no-alff',  dest='alff',  action='store_false',
                        help='Disable ALFF')
    parser.add_argument('--no-malff', dest='malff', action='store_false',
                        help='Disable mALFF')
    parser.add_argument('--no-falff', dest='falff', action='store_false',
                        help='Disable fALFF')
    parser.add_argument('--no-rsfa',  dest='rsfa',  action='store_false',
                        help='Disable RSFA')
    parser.set_defaults(alff=True, malff=True, falff=True, rsfa=True)
    parser.add_argument('--low',     type=float, default=0.01,
                        help='Lower bandpass frequency')
    parser.add_argument('--high',    type=float, default=0.08,
                        help='Upper bandpass frequency')
    args = parser.parse_args()

    metrics = []
    if args.alff:  metrics.append('ALFF')
    if args.malff: metrics.append('mALFF')
    if args.falff: metrics.append('fALFF')
    if args.rsfa:  metrics.append('RSFA')
    if not metrics:
        raise ValueError("No metrics enabled. At least one must be enabled.")

    # Ensure input file exists
    if not os.path.exists(args.fmri):
        print(f"Error: Input fMRI file {args.fmri} does not exist")
        sys.exit(1)
        
    # Check mask if provided
    if args.mask and not os.path.exists(args.mask):
        print(f"Error: Mask file {args.mask} does not exist")
        sys.exit(1)

    # Ensure output dir
    os.makedirs(args.output, exist_ok=True)

    # Create temporary directory for AFNI outputs
    temp_dir = tempfile.mkdtemp(prefix="alff_")
    prefix = os.path.join(temp_dir, "rsfc_")
    
    try:
        # Run AFNI
        run_3dRSFC_all(
            fmri=args.fmri,
            mask=args.mask,
            band_low=args.low,
            band_high=args.high,
            metrics=metrics,
            prefix=prefix
        )

        # Save each AFNI output
        for m in metrics:
            afni_head = f"{prefix}_{m}+orig.HEAD"
            nii_path   = os.path.join(args.output, f"{m.lower()}.nii.gz")
            
            # Check if AFNI output exists
            if not os.path.exists(afni_head):
                print(f"Warning: Expected output {afni_head} not found, skipping {m}")
                continue
                
            # Convert to NIfTI if needed
            try:
                if not os.path.exists(afni_head.replace('.HEAD','.nii')):
                    subprocess.run(
                        ["3dAFNItoNIFTI", afni_head],
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                
                if not os.path.exists(afni_head.replace('.HEAD','.nii')):
                    print(f"Warning: Conversion to NIfTI failed for {m}")
                    continue
                    
                img = nib.load(afni_head.replace('.HEAD','.nii'))
                nib.save(img, nii_path)
                print(f"[+] Saved {m} map to {nii_path}")
            except Exception as e:
                print(f"Error processing {m}: {e}")
                
    finally:
        # Clean up temporary files
        try:
            shutil.rmtree(temp_dir)
        except Exception as e:
            print(f"Warning: Failed to clean up temporary directory {temp_dir}: {e}")

if __name__ == '__main__':
    main()
