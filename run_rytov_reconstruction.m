%% Rytov Tomogram Reconstruction from Experimental Data
% Script to perform Rytov reconstruction from background and sample image pairs
%
% This script:
% 1. Loads configuration from field_retrieval_config.json
% 2. Processes each background-sample pair specified in the configuration
% 3. Performs in-memory field retrieval using FIELD_EXPERIMENTAL_RETRIEVAL (not saved)
% 4. Applies Rytov reconstruction algorithm
% 5. Saves only the final Rytov tomograms (RI_rytov)
%
% Configuration file should contain:
%   - data_path: Path to merged PNG stack data
%   - output_path: Path for saving Rytov reconstruction results
%   - optical_parameters: Optical system parameters (wavelength, NA, RI, etc.)
%   - processing_parameters: Processing options (cutout_portion, GPU usage, etc.)
%   - reconstruction_parameters: Reconstruction options (zsize)
%   - sample_pairs: Array of background-sample pairs to process
%
% See also: FIELD_EXPERIMENTAL_RETRIEVAL, BACKWARD_SOLVER_RYTOV

clc; clear; close all;

%% Add paths
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir));

% Add preprocessing utilities for loading PNG stacks
preprocessing_path = fullfile(fileparts(current_dir), '00_preprocessing');
addpath(genpath(preprocessing_path));

fprintf('Rytov Tomogram Reconstruction from Experimental Data\n');

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

fprintf('Processing %d sample pairs...\n', length(config.sample_pairs));

%% Create output directory
if ~exist(config.output_path, 'dir')
    mkdir(config.output_path);
end

%% Process each sample pair
num_pairs = length(config.sample_pairs);

for pair_idx = 1:num_pairs
    current_pair = config.sample_pairs(pair_idx);
    fprintf('[%d/%d] %s... ', pair_idx, num_pairs, current_pair.output_name);

    % Check if results already exist
    output_dir = fullfile(config.output_path, current_pair.output_name);
    rytov_file = fullfile(output_dir, 'RI_rytov.mat');

    if exist(rytov_file, 'file')
        fprintf('Skipped (already exists)\n');
        continue;
    end

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

    % Setup field retrieval parameters following C_albicans.m pattern
    tic;
    field_retrieval_params = struct();
    field_retrieval_params.NA = config.optical_parameters.NA;
    field_retrieval_params.wavelength = config.optical_parameters.wavelength;
    field_retrieval_params.RI_bg = config.optical_parameters.RI_bg;
    field_retrieval_params.resolution = config.optical_parameters.resolution;
    field_retrieval_params.resolution_image = config.optical_parameters.resolution_image;
    field_retrieval_params.vector_simulation = config.optical_parameters.vector_simulation;
    field_retrieval_params.use_abbe_sine = config.optical_parameters.use_abbe_sine;
    field_retrieval_params.use_abbe_correction = config.optical_parameters.use_abbe_correction;
    field_retrieval_params.cutout_portion = config.processing_parameters.cutout_portion;
    field_retrieval_params.other_corner = config.processing_parameters.other_corner;
    field_retrieval_params.conjugate_field = config.processing_parameters.conjugate_field;
    field_retrieval_params.verbose = config.processing_parameters.verbose;
    field_retrieval_params.normalidx = config.processing_parameters.normalidx;
    field_retrieval_params.use_GPU = config.processing_parameters.use_GPU;

    % Create field retrieval object and process fields
    field_retrieval = FIELD_EXPERIMENTAL_RETRIEVAL(field_retrieval_params);
    [input_field, output_field, rytov_params] = field_retrieval.get_fields(background_stack, sample_stack);

    % Clear stacks to save memory
    clear background_stack sample_stack;

    % Configure Rytov solver parameters
    rytov_params.use_non_negativity = false;
    rytov_params.non_negativity_iteration = 100;

    % Set reconstruction z-size from configuration
    rytov_params.size(3) = config.reconstruction_parameters.zsize;

    % Create Rytov solver and perform reconstruction
    rytov_solver = BACKWARD_SOLVER_RYTOV(rytov_params);
    [RI_rytov, ORytov] = rytov_solver.solve(input_field, output_field);

    % Clear field data to save memory
    clear input_field output_field ORytov;

    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Save Rytov reconstruction result
    save(rytov_file, 'RI_rytov', '-v7.3');

    elapsed_time = toc;
    fprintf('Done (%.1fs)\n', elapsed_time);

    % Clear large variables to free memory
    clear RI_rytov;
end

%% Summary
fprintf('Complete. Results saved to: %s\n', config.output_path);
