%% Field Retrieval from Experimental Data
% Main script to perform field retrieval from background and sample image pairs
%
% This script:
% 1. Loads configuration from field_retrieval_config.json
% 2. Processes each background-sample pair specified in the configuration
% 3. Saves retrieved fields, amplitude, and phase to output directory
%
% Configuration file should contain:
%   - data_path: Path to merged PNG stack data
%   - output_path: Path for saving retrieval results
%   - optical_parameters: Optical system parameters (wavelength, NA, RI, etc.)
%   - processing_parameters: Processing options (cutout_portion, GPU usage, etc.)
%   - sample_pairs: Array of background-sample pairs to process
%
% See also: process_field_pair, save_field_results, visualize_field_results

clc; clear; close all;

%% Add paths
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir));

% Add preprocessing utilities for loading PNG stacks
preprocessing_path = fullfile(fileparts(current_dir), '00_preprocessing');
addpath(genpath(preprocessing_path));

fprintf('Field Retrieval from Experimental Data\n');

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
                   'processing_parameters', 'sample_pairs'};
for i = 1:length(required_fields)
    if ~isfield(config, required_fields{i})
        error('Configuration missing required field: %s', required_fields{i});
    end
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
        k0s_file = fullfile(output_dir, 'k0s.mat');
        output_field_file = fullfile(output_dir, 'output_field.mat');

        if exist(k0s_file, 'file') && exist(output_field_file, 'file')
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

        % Process field pair
        tic;
        [input_field, output_field, updated_params, k0s] = process_field_pair(...
            background_stack, ...
            sample_stack, ...
            config.optical_parameters, ...
            config.processing_parameters);

        % Save results
        save_field_results(...
            input_field, ...
            output_field, ...
            config.output_path, ...
            current_pair.output_name, ...
            updated_params);

        % Also save k0s if computed
        if exist('k0s', 'var')
            k0s_path = fullfile(config.output_path, current_pair.output_name, 'k0s.mat');
            save(k0s_path, 'k0s');
        end

        % Optional visualization
        if config.processing_parameters.verbose
            visualize_field_results(input_field, output_field);
            fig_path = fullfile(config.output_path, current_pair.output_name, 'visualization.png');
            saveas(gcf, fig_path);
        end

        elapsed_time = toc;
        fprintf('Done (%.1fs)\n', elapsed_time);

    % Clear large variables to free memory
    clear background_stack sample_stack input_field output_field k0s;
end

%% Summary
fprintf('Complete. Results saved to: %s\n', config.output_path);
