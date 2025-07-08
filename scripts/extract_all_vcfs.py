import os
import tarfile
from pathlib import Path

# Source directories where the tar files are stored
source_dirs = [
    "/path/to/ukb/stored/batch_01_to_batch_20",
    "/path/to/ukb/stored/batch_21_to_batch_50",
]

# Destination directory where all VCF files will be organized
destination_dir = "/path/to/ukb_vcfs/"
Path(destination_dir).mkdir(parents=True, exist_ok=True)

# Process each source directory
for source_dir in source_dirs:
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            # Process only batch_10.tar to batch_20.tar files
            if file.startswith("batch_") and file.endswith(".tar"):
                # Extract the batch number from the filename
                batch_number = int(file.split('_')[1].split('.')[0])
                if 10 <= batch_number <= 20:  
                    tar_path = os.path.join(root, file)
                    
                    # Decompress the tar file into a specific batch directory
                    with tarfile.open(tar_path, "r") as tar:
                        batch_dir = os.path.join(destination_dir, Path(tar_path).stem)
                        Path(batch_dir).mkdir(parents=True, exist_ok=True)
                        tar.extractall(batch_dir)
                        
                        # Process each "run_*" subdirectory within the batch directory
                        for run_root, run_dirs, run_files in os.walk(batch_dir):
                            for run_file in run_files:
                                # Process only VCF files (with or without .gz)
                                if run_file.endswith(".vcf") or run_file.endswith(".vcf.gz"):
                                    source_vcf = os.path.join(run_root, run_file)
                                    
                                    # Extract the sample ID (first part of the filename, before the first underscore)
                                    sample_id = run_file.split('_')[0]
                                    
                                    # Create a subdirectory based on the first two digits of the sample ID
                                    subdir = os.path.join(destination_dir, sample_id[:2])  # First two digits of sample ID
                                    
                                    # Create the sample ID subdirectory if it doesn't exist
                                    Path(subdir).mkdir(parents=True, exist_ok=True)
                                    
                                    # Move or symlink the VCF file into the corresponding sample subdirectory
                                    dest_vcf = os.path.join(subdir, run_file)
                                    
                                    # If you want to symlink the files instead of copying them (to save space):
                                    if not os.path.exists(dest_vcf):
                                        os.symlink(source_vcf, dest_vcf)
                                    else:
                                        # If the symlink exists, you could either skip or overwrite (based on your preference)
                                        print(f"File {dest_vcf} already exists, skipping.")

print("VCF files have been organized into sample-based subdirectories.")
