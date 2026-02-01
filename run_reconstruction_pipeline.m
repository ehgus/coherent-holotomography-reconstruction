%% Integrated Field Retrieval and Rytov Tomogram Reconstruction
% Script to perform field retrieval followed by Rytov reconstruction in a single pipeline
%
% This script:
% 1. Loads configuration from field_retrieval_config.json
% 2. Processes each background-sample pair specified in the configuration
% 3. Performs field retrieval using extract_complex_field() (NOT saved)
% 4. Applies Rytov reconstruction algorithm
% 5. Saves only the final Rytov tomograms (RI.mat) to output directory
%
% Configuration file should contain:
%   - data_path: Path to merged PNG stack data
%   - output_path: Path for saving Rytov reconstruction results
%   - optical_parameters: Optical system parameters (wavelength, NA, RI, etc.)
%   - processing_parameters: Processing options (cutout_portion, GPU usage, etc.)
%   - reconstruction_parameters: Reconstruction options (zsize)
%   - sample_pairs: Array of background-sample pairs to process
%
% Note: Field information is NOT saved to disk. Only the final Rytov tomogram (RI.mat)
%       is saved to config.output_path/sample_name/ directory.
%
% See also: extract_complex_field, BACKWARD_SOLVER_RYTOV

clc; clear; close all;

%% Add paths based on pwd()
current_dir = pwd();
addpath(genpath(current_dir));

fprintf('Integrated Field Retrieval and Rytov Tomogram Reconstruction\n');
fprintf('Pipeline: PNG stacks -> extract_complex_field -> BACKWARD_SOLVER_RYTOV -> RI.mat\n\n');

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
required_fields = {'data_path', 'output_path', 'optical_parameters', ...
                   'processing_parameters', 'reconstruction_parameters', 'sample_pairs'};
for i = 1:length(required_fields)
    if ~isfield(config, required_fields{i})
        error('Configuration missing required field: %s', required_fields{i});
    end
end

% Validate reconstruction parameters
if ~isfield(config.reconstruction_parameters, 'zsize')
    error('Configuration missing required field: reconstruction_parameters.zsize');
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
    
    % Check if Rytov result already exists
    output_dir = fullfile(config.output_path, current_pair.output_name);
    rytov_file = fullfile(output_dir, 'RI.mat');
    
    if exist(rytov_file, 'file')
        fprintf(' - Skipped (RI.mat already exists)\n');
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
    if height ~= width
        min_size = min(height, width);
        if height > width
            % Crop height (rows)
            crop_amount = height - min_size;
            crop_top = floor(crop_amount / 2);
            crop_bottom = crop_amount - crop_top;
            background_stack = background_stack(crop_top+1:end-crop_bottom, :, :);
            sample_stack = sample_stack(crop_top+1:end-crop_bottom, :, :);
        else
            % Crop width (columns)
            crop_amount = width - min_size;
            crop_left = floor(crop_amount / 2);
            crop_right = crop_amount - crop_left;
            background_stack = background_stack(:, crop_left+1:end-crop_right, :);
            sample_stack = sample_stack(:, crop_left+1:end-crop_right, :);
        end
    end
    
    % ========================================================================
    % STEP 1: FIELD RETRIEVAL (using extract_complex_field - NOT saved)
    % ========================================================================
    fprintf(' -> Processing fields');
    tic;
    
    % output_field: 3D [height x width x num_images] complex array
    % updated_params: updated optical parameters after field retrieval
    % k0s: wavenumber information for each illumination angle
    [output_field, updated_params, k0s] = extract_complex_field(...
        background_stack, ...
        sample_stack, ...
        config.optical_parameters, ...
        config.processing_parameters);
    
    % Clear stacks to save memory
    clear background_stack sample_stack;
    

    
    % ========================================================================
    % STEP 2: PREPARE RYTOV PARAMETERS
    % ========================================================================
    fprintf(' -> Preparing Rytov reconstruction');
    
    % Create Rytov parameters from updated_params
    rytov_params = struct();
    rytov_params.size = updated_params.size;
    rytov_params.wavelength = updated_params.wavelength;
    rytov_params.NA = updated_params.NA;
    rytov_params.RI_bg = updated_params.RI_bg;
    rytov_params.resolution = updated_params.resolution;
    rytov_params.vector_simulation = updated_params.vector_simulation;
    rytov_params.use_abbe_sine = updated_params.use_abbe_sine;
    rytov_params.use_GPU = updated_params.use_GPU;
    
    % Set reconstruction z-size from configuration
    rytov_params.size(3) = config.reconstruction_parameters.zsize;

    % ========================================================================
    % STEP 3: PERFORM RYTOV RECONSTRUCTION
    % ========================================================================
    fprintf(' -> Executing Rytov solver');
    
    % Create Rytov solver and perform reconstruction
    rytov_solver = BACKWARD_SOLVER_RYTOV(rytov_params);
    [RI, ~] = rytov_solver.solve(output_field);
    
    % Clear field data to save memory
    clear output_field updated_params k0s;
    
    % ========================================================================
    % STEP 4: SAVE RESULTS
    % ========================================================================
    fprintf(' -> Saving results');
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Save Rytov reconstruction result ONLY (no field data saved)
    save(rytov_file, 'RI', '-v7.3');
    
    elapsed_time = toc;
    
    % Clear large variables to free memory
    clear RI;
    
    fprintf(' - Done (%.1fs)\n', elapsed_time);
end

%% Summary
fprintf('\n========================================\n');
fprintf('Complete. Results saved to: %s\n', config.output_path);
fprintf('Files saved: RI.mat (Rytov tomogram only)\n');
fprintf('========================================\n');
