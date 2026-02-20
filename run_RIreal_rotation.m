%% Run RIreal Rotation
% This script applies 3D rotation to RI_rytov.mat files in all sample directories
% configured in a JSON configuration file, similar to run_tv_optimization_coherent.m
%
% The script:
% 1. Loads a configuration file (JSON)
% 2. Finds all sample directories based on data_dir_pattern
% 3. Rotates RIreal data from RI_rytov.mat files
% 4. Saves rotated results as RI_rytov_rotated.mat

clear; clc;
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, "utils")));

%% Load Configuration
[config_file, config_path] = uigetfile('*.json', 'Select Configuration File');

if isequal(config_file, 0)
    error('No configuration file selected. Exiting.');
end

config_fullpath = fullfile(config_path, config_file);

% Read and parse JSON configuration
config_text = fileread(config_fullpath);
config = jsondecode(config_text);

% Extract configuration parameters
data_dir_pattern = config.data_dir_pattern;

%% Get rotation parameters from user
prompt = {'Rotation angle (degrees):', 'Rotation axis (x/y/z):'};
dlg_title = 'Rotation Parameters';
num_lines = 1;
default_answer = {'-15', 'x'};

user_input = inputdlg(prompt, dlg_title, num_lines, default_answer);

if isempty(user_input)
    error('User cancelled input dialog. Exiting.');
end

angle_deg = str2double(user_input{1});
axis = lower(strtrim(user_input{2}));

% Validate inputs
valid_axes = {'x', 'y', 'z'};
if ~ismember(axis, valid_axes)
    error('Invalid rotation axis. Must be one of: x, y, z');
end

if isnan(angle_deg)
    error('Invalid rotation angle. Must be a numeric value.');
end

%% Find reconstruction directories and apply rotation

sample_dirs = dir(data_dir_pattern);
sample_dirs = sample_dirs([sample_dirs.isdir]);
sample_dirs = sample_dirs(~ismember({sample_dirs.name}, {'.', '..'}));

task_count = 0;
fprintf('Starting rotation of RI data files...\n');
fprintf('Angle: %.1f degrees, Axis: %s\n', angle_deg, upper(axis));
fprintf('=================================================\n');

for i = 1:length(sample_dirs)
    sample_path = fullfile(sample_dirs(i).folder, sample_dirs(i).name);
    
    % Check for existing reconstruction file
    reconstruction_file = fullfile(sample_path, 'RI_rytov.mat');
    if ~isfile(reconstruction_file)
        fprintf('[%d/%d] Skipping "%s" - RI_rytov.mat not found\n', i, length(sample_dirs), sample_dirs(i).name);
        continue;
    end
    
    % Define output file path
    rotated_output_file = fullfile(sample_path, 'RI_rytov_rotated.mat');
    if isfile(rotated_output_file)
        fprintf('[%d/%d] Skipping "%s" - rotated file already exists\n', i, length(sample_dirs), sample_dirs(i).name);
        continue;
    end
    
    % Load reconstruction data
    fprintf('[%d/%d] Processing "%s"...\n', i, length(sample_dirs), sample_dirs(i).name);
    fprintf('  Loading RI_rytov.mat...\n');
    recon_data = load(reconstruction_file);
    RIreal = recon_data.RIreal;
    
    fprintf('  Original RIreal size: %d x %d x %d\n', size(RIreal, 1), size(RIreal, 2), size(RIreal, 3));
    
    % Apply rotation
    fprintf('  Rotating RIreal by %.1f degrees around %s-axis...\n', angle_deg, upper(axis));
    task_start = tic;
    RIreal_rotated = rotate3D(RIreal, angle_deg, axis);
    elapsed_time = toc(task_start);
    fprintf('  Rotation completed in %.2f seconds\n', elapsed_time);
    
    fprintf('  Rotated RIreal size: %d x %d x %d\n', size(RIreal_rotated, 1), size(RIreal_rotated, 2), size(RIreal_rotated, 3));
    
    % Save rotated data
    fprintf('  Saving rotated RIreal to: %s\n', rotated_output_file);
    save(rotated_output_file, 'RIreal_rotated', '-v7.3');
    task_count = task_count + 1;
    fprintf('  Completed successfully\n');
end

fprintf('=================================================\n');
fprintf('Total files processed: %d\n', task_count);
fprintf('Rotation operation completed.\n');


