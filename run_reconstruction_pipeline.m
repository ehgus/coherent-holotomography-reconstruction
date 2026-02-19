%% Integrated Field Retrieval and Tomogram Reconstruction
% Script to perform field retrieval followed by Tomogram reconstruction in a single pipeline
%
% This script:
% 1. Loads configuration from JSON file (e.g., smooth_TGV.json)
% 2. Processes each background-sample pair specified in the configuration
% 3. Performs field retrieval using extract_complex_field() (NOT saved)
% 4. Applies both Rytov and Born Tomogram reconstruction algorithms
% 5. Saves the final Tomograms (RI_rytov.mat and RI_born.mat) to output directory
%
% Configuration file should contain:
%   - data_path: Path to merged PNG stack data
%   - output_path: Path for saving Tomogram reconstruction results
%   - imaging_condition: Imaging parameters (wavelength, NA, RI_bg, resolution)
%   - field_generator_condition: Field processing options (cutout_portion, crop settings, etc.)
%   - tomography_generator_condition: Reconstruction options (resolution, zsize_micron)
%   - sample_pairs: Array of background-sample pairs to process
%
% Note: Field information is NOT saved to disk. Only the final Tomograms (RI_rytov.mat, RI_born.mat)
%       are saved to config.output_path/sample_name/ directory.
%
% See also: extract_complex_field, BACKWARD_SOLVER_RYTOV, BACKWARD_SOLVER_BORN

clc; clear; close all;

%% Add paths
current_dir = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(current_dir));

fprintf('Integrated Field Retrieval and Tomogram Reconstruction\n');
fprintf('Pipeline: PNG stacks -> extract_complex_field -> BACKWARD_SOLVER (Rytov & Born) -> RI_rytov.mat, RI_born.mat\n\n');

%% Load configuration
[config_file, config_path] = uigetfile('*.json', 'Select configuration file');

if config_file == 0
    fprintf('No configuration file selected. Exiting...\n');
    return;
end

config_filepath = fullfile(config_path, config_file);

config_text = fileread(config_filepath);
config = jsondecode(config_text);

% Validate configuration
required_fields = {'data_path', 'output_path', 'imaging_condition', ...
                   'field_generator_condition', 'tomography_generator_condition', 'sample_pairs'};
for i = 1:length(required_fields)
    if ~isfield(config, required_fields{i})
        error('Configuration missing required field: %s', required_fields{i});
    end
end

% Validate reconstruction parameters
if ~isfield(config.tomography_generator_condition, 'zsize_micron')
    error('Configuration missing required field: tomography_generator_condition.zsize_micron');
end

fprintf('Processing %d sample pair(s)...\n', length(config.sample_pairs));

%% Create output directory
if ~exist(config.output_path, 'dir')
    mkdir(config.output_path);
end

%% Process each sample pair
num_pairs = length(config.sample_pairs);

for pair_idx = 1:num_pairs
    current_pair = config.sample_pairs(pair_idx);
    fprintf('\n[%d/%d] %s', pair_idx, num_pairs, current_pair.output_name);
    
    % Check if Tomogram results already exist
    output_dir = fullfile(config.output_path, current_pair.output_name);
    tomogram_file_rytov = fullfile(output_dir, 'RI_rytov.mat');
    tomogram_file_born = fullfile(output_dir, 'RI_born.mat');
    
    if exist(tomogram_file_rytov, 'file') && exist(tomogram_file_born, 'file')
        fprintf(' - Skipped (RI_rytov.mat and RI_born.mat already exist)\n');
        continue;
    end
    
    fprintf(' - Loading stacks');
    
    % Load background and sample stacks
    bg_path = fullfile(config.data_path, current_pair.background, 'data3d');
    sample_path = fullfile(config.data_path, current_pair.sample, 'data3d');
    
    if ~exist(bg_path, 'dir')
        error('Background path does not exist: %s', bg_path);
    end
    background_stack = load_png_stack_data(bg_path);
    
    if ~exist(sample_path, 'dir')
        error('Sample path does not exist: %s', sample_path);
    end
    sample_stack = load_png_stack_data(sample_path);
    
    % Verify stacks have same dimensions
    if ~isequal(size(background_stack), size(sample_stack))
        error('Background and sample stacks must have the same dimensions');
    end
    
    % Crop to make square (X and Y dimensions equal)
    [height, width, num_images] = size(background_stack);
    offset = config.field_generator_condition.crop_offset_micron;
    fov = config.field_generator_condition.crop_fov_micron;
    img_resolution = config.imaging_condition.resolution;
    img_pixel_size = [size(sample_stack, 1); size(sample_stack, 2)];

    fov_min = round((offset - fov/2)./img_resolution + img_pixel_size/2)+1;
    fov_max = round((offset + fov/2)./img_resolution + img_pixel_size/2);

    background_stack = background_stack(fov_min(1):fov_max(1),fov_min(2):fov_max(2),:);
    sample_stack = sample_stack(fov_min(1):fov_max(1),fov_min(2):fov_max(2),:);
    
    % ========================================================================
    % STEP 1: FIELD RETRIEVAL (using extract_complex_field - NOT saved)
    % ========================================================================
    fprintf(' -> Processing fields');
    tic;
    field_generator_condition = config.field_generator_condition;
     if isfield(current_pair, 'refocus_z_micron') && current_pair.refocus_z_micron ~= 0
        field_generator_condition.refocus_z_micron = current_pair.refocus_z_micron;
     end
    % output_field: 3D [height x width x num_images] complex array
    % updated_params: updated optical parameters after field retrieval
    % k0s: wavenumber information for each illumination angle
    [output_field, updated_params, illum_k0] = extract_complex_field(...
        background_stack, ...
        sample_stack, ...
        config.imaging_condition, ...
        field_generator_condition, ...
        config.tomography_generator_condition);
    % Clear stacks to save memory
    clear background_stack sample_stack;
    
    % ========================================================================
    % STEP 2: PREPARE TOMOGRAM PARAMETERS
    % ========================================================================
    fprintf(' -> Preparing Tomogram reconstruction');
    
    % Create Tomogram parameters from updated_params
    tomo_params = struct();
    tomo_params.wavelength = updated_params.wavelength;
    tomo_params.NA = updated_params.NA;
    tomo_params.RI_bg = updated_params.RI_bg;
    tomo_params.field_resolution = updated_params.resolution;
    
    % Note: output_field XY dimensions are already resized to match tomography resolution
    % Only need to set Z dimension based on configuration
    tomo_params.tomogram_size = [size(output_field, 1); 
                                  size(output_field, 2); 
                                  round(config.tomography_generator_condition.zsize_micron / config.tomography_generator_condition.resolution(3))];                       
    tomo_params.tomogram_resolution = [updated_params.resolution(:); 
                                        config.tomography_generator_condition.resolution(3)];    

    % ========================================================================
    % STEP 3 & 4: PERFORM TOMOGRAM RECONSTRUCTION AND SAVE (Rytov and Born)
    % ========================================================================
    fprintf(' -> Executing Tomogram solvers');
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Resolution for saving
    resolution = tomo_params.tomogram_resolution;
    
    % Define solvers and output files
    solvers = {
        'BACKWARD_SOLVER_RYTOV', tomogram_file_rytov;
        'BACKWARD_SOLVER_BORN', tomogram_file_born
    };
    
    % Process each solver
    for solver_idx = 1:size(solvers, 1)
        solver_name = solvers{solver_idx, 1};
        solver_file = solvers{solver_idx, 2};
        
        % Skip if file already exists
        if exist(solver_file, 'file')
            continue;
        end
        
        fprintf(' (%s)', solver_name);
        
        % Create solver instance
        if strcmp(solver_name, 'BACKWARD_SOLVER_RYTOV')
            tomogram_solver = BACKWARD_SOLVER_RYTOV(tomo_params);
        else
            tomogram_solver = BACKWARD_SOLVER_BORN(tomo_params);
        end
        
        % Perform reconstruction
        [potential, fourier_mask] = tomogram_solver.solve(output_field, illum_k0);
        RIreal = real(potential2RI(potential*4*pi, tomo_params.wavelength, tomo_params.RI_bg));
        
        % Save results
        save(solver_file, 'potential', 'RIreal', 'resolution', 'fourier_mask', '-v7.3');
        
        % Clear variables
        clear potential RIreal fourier_mask tomogram_solver;
    end
    
    % Clear field data to save memory
    clear output_field updated_params illum_k0;
    
    elapsed_time = toc;
    
    fprintf(' - Done (%.1fs)\n', elapsed_time);
end

%% Summary
fprintf('\n========================================\n');
fprintf('Complete. Results saved to: %s\n', config.output_path);
fprintf('Files saved: RI_rytov.mat, RI_born.mat (Tomogram only)\n');
fprintf('========================================\n');
