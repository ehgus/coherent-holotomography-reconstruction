# Field Retrieval

This directory contains scripts and functions for performing field retrieval from experimental holographic data.

## Directory Structure

```
01_field_retrieval/
├── run_field_retrieval.m           # Main script
├── configuration/
│   └── field_retrieval_config.json # Configuration file
├── utils/
│   └── field_processing/
│       ├── process_field_pair.m         # Core field processing function
│       ├── save_field_results.m         # Save results to disk
│       └── visualize_field_results.m    # Visualization utilities
└── README.md
```

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
- **processing_parameters**: Processing options
  - `cutout_portion`: Fourier space cropping (0-0.5, typically 0.333)
  - `use_GPU`: Enable GPU acceleration (true/false)
  - `verbose`: Show visualization during processing (true/false)
- **sample_pairs**: Array of background-sample pairs to process

### 2. Run the Script

Open MATLAB and run:

```matlab
cd('path/to/01_field_retrieval')
run_field_retrieval
```

### 3. Output

For each sample pair, the script creates a subdirectory in `output_path` containing:

- `input_field.mat`: Background field (complex-valued)
- `output_field.mat`: Sample field (complex-valued)
- `parameters.mat`: Updated optical parameters
- `amplitude.mat`: Amplitude ratio (|output_field / input_field|)
- `phase.mat`: Phase difference (angle(output_field / input_field))
- `k0s.mat`: Peak positions in Fourier space
- `visualization.png`: (if verbose=true) Visual summary

## Algorithm

The field retrieval process follows these steps (based on FIELD_EXPERIMENTAL_RETRIEVAL.m):

1. **Fourier Transform**: Convert images to Fourier space
2. **Centering**: Find and center the main diffraction peak
3. **Resampling**: Resize to match desired spatial resolution
4. **NA Filtering**: Apply numerical aperture mask
5. **Inverse Transform**: Convert back to real space
6. **Phase Unwrapping**: Remove 2π phase discontinuities
7. **Phase Correction**: Remove linear phase tilts

## Dependencies

- MATLAB R2019b or later
- Image Processing Toolbox
- (Optional) Parallel Computing Toolbox for GPU acceleration
- `load_png_stack_data.m` from 00_preprocessing/utils/io/

## Example Configuration

```json
{
    "data_path": "F:\\Data\\merged_stacks",
    "output_path": "F:\\Data\\field_retrieval_results",
    "optical_parameters": {
        "wavelength": 0.532,
        "NA": 1.2,
        "RI_bg": 1.336,
        "resolution": [0.1, 0.1, 0.1],
        "resolution_image": [0.1, 0.1]
    },
    "processing_parameters": {
        "cutout_portion": 0.333,
        "use_GPU": true,
        "verbose": false
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
