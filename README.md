# Field Retrieval and Rytov Reconstruction

This directory contains scripts and functions for performing field retrieval from experimental holographic data, followed by Rytov tomogram reconstruction.

## Directory Structure

```
rytov_solver-matlab/
├── run_reconstruction_pipeline.m   # Main script (Field Retrieval + Rytov Reconstruction)
├── configuration/
│   └── field_retrieval_config.json # Configuration file
├── utils/
│   ├── field_processing/
│   │   ├── extract_complex_field.m         # Core field retrieval function
│   │   └── visualize_field_results.m       # Visualization utilities
│   ├── ht_reconstruction/
│   │   ├── BACKWARD_SOLVER_RYTOV.m         # Rytov reconstruction solver
│   │   ├── BACKWARD_SOLVER.m               # Standard reconstruction solver
│   │   ├── DERIVE_OPTICAL_TOOL.m           # Optical tool derivation
│   │   ├── potential2RI.m                  # Convert potential to refractive index
│   │   └── update_struct.m                 # Utility for struct updates
│   └── io/
│       ├── load_bin_data.m                 # Load binary data
│       ├── load_png_stack_data.m           # Load PNG stacks
│       └── save_png_stack_data.m           # Save PNG stacks
└── README.md
```

## Pipeline Overview

```
PNG stacks (background + sample)
        ↓
Field Retrieval (extract_complex_field)
        ↓
Complex field data (NOT saved)
        ↓
Rytov Reconstruction (BACKWARD_SOLVER_RYTOV)
        ↓
Refractive Index (RI.mat) - SAVED
```

The pipeline performs two main steps:
1. **Field Retrieval**: Extracts complex field from background and sample PNG stacks
2. **Rytov Reconstruction**: Performs 3D tomographic reconstruction to obtain refractive index distribution

## Usage

### 1. Configure Parameters

Edit `configuration/field_retrieval_config.json` to specify:

- **data_path**: Path to merged PNG stack data (output from 00_preprocessing)
- **output_path**: Directory where results will be saved
- **optical_parameters**: System parameters
  - `wavelength`: Illumination wavelength in microns (e.g., 0.532)
  - `NA`: Numerical aperture (e.g., 1.2)
  - `RI_bg`: Background refractive index (e.g., 1.336)
  - `resolution`: [dx, dy, dz] spatial resolution in microns
  - `resolution_image`: [dx, dy] image pixel size in microns
  - `vector_simulation`: Enable vector simulation (true/false)
  - `use_abbe_sine`: Use Abbe sine condition (true/false)
- **processing_parameters**: Processing options
  - `cutout_portion`: Fourier space cropping (0-0.5, typically 0.333)
  - `use_GPU`: Enable GPU acceleration (true/false)
- **reconstruction_parameters**: Reconstruction options
  - `zsize`: Z-dimension size for 3D reconstruction volume
- **sample_pairs**: Array of background-sample pairs to process

### 2. Run the Script

Open MATLAB and run:

```matlab
cd('path/to/rytov_solver-matlab')
run_reconstruction_pipeline
```

### 3. Output

For each sample pair, the script creates a subdirectory in `output_path` containing:

- `RI.mat`: Refractive index tomogram (3D array)

Note: Intermediate field data is NOT saved to disk, only the final Rytov reconstruction result is saved.

## Algorithm

The reconstruction pipeline follows these steps:

### Stage 1: Field Retrieval (based on extract_complex_field.m)

1. **Fourier Transform**: Convert images to Fourier space
2. **Centering**: Find and center the main diffraction peak
3. **Resampling**: Resize to match desired spatial resolution
4. **NA Filtering**: Apply numerical aperture mask
5. **Inverse Transform**: Convert back to real space
6. **Phase Unwrapping**: Remove 2π phase discontinuities
7. **Phase Correction**: Remove linear phase tilts

### Stage 2: Rytov Reconstruction (BACKWARD_SOLVER_RYTOV.m)

1. **Forward Propagation**: Propagate field through sample volume
2. **Rytov Approximation**: Apply Rytov scattering theory
3. **Inverse Problem Solving**: Recover refractive index distribution
4. **3D Reconstruction**: Build 3D tomogram from multiple illumination angles

## Dependencies

- MATLAB R2019b or later
- Image Processing Toolbox
- (Optional) Parallel Computing Toolbox for GPU acceleration
- `load_png_stack_data.m` from 00_preprocessing/utils/io/

## Example Configuration

```json
{
    "data_path": "F:\\Data\\merged_stacks",
    "output_path": "F:\\Data\\reconstruction_results",
    "optical_parameters": {
        "wavelength": 0.532,
        "NA": 1.2,
        "RI_bg": 1.336,
        "resolution": [0.1, 0.1, 0.1],
        "resolution_image": [0.1, 0.1],
        "vector_simulation": false,
        "use_abbe_sine": true
    },
    "processing_parameters": {
        "cutout_portion": 0.333,
        "use_GPU": true
    },
    "reconstruction_parameters": {
        "zsize": 256
    },
    "sample_pairs": [
        {
            "background": "bg_IMAGING",
            "sample": "sample_IMAGING",
            "output_name": "sample_001"
        }
    ]
}
```

## Notes

- Ensure background and sample images have been preprocessed and merged using scripts in `00_preprocessing/`
- The script expects PNG stacks in `data3d/` subdirectories
- GPU acceleration requires a CUDA-compatible GPU
- For large datasets, processing may take several minutes per sample pair
- Intermediate field data is automatically cleared from memory to save disk space; only RI.mat (Rytov reconstruction result) is saved
- If RI.mat already exists for a sample, it will be skipped (enable rapid re-runs with modified parameters)
