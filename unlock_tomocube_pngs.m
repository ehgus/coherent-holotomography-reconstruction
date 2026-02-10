%% Unlock Tomocube .dat -> PNG stacks using JSON config
% Select a JSON config used by run_reconstruction_pipeline.m.

clc; clear;

current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir));

[config_file, config_path] = uigetfile('*.json', 'Select configuration file');
if config_file == 0
    fprintf('No configuration file selected. Exiting...\n');
    return;
end

config_filepath = fullfile(config_path, config_file);
config_text = fileread(config_filepath);
config = jsondecode(config_text);

required_fields = {'data_path', 'sample_pairs'};
for i = 1:length(required_fields)
    if ~isfield(config, required_fields{i})
        error('Configuration missing required field: %s', required_fields{i});
    end
end

num_pairs = length(config.sample_pairs);
fprintf('Unlocking %d sample pair(s)...\n', num_pairs);

for pair_idx = 1:num_pairs
    current_pair = config.sample_pairs(pair_idx);

    sample_dir = fullfile(config.data_path, current_pair.sample);
    bg_dir = fullfile(config.data_path, current_pair.background);

    data_dat = fullfile(sample_dir, 'images.dat');
    data_png = fullfile(sample_dir, 'data3d');

    bg_dat = fullfile(bg_dir, 'calibration.dat');
    bg_png = fullfile(bg_dir, 'data3d');

    fprintf('[%d/%d] sample: %s\n', pair_idx, num_pairs, current_pair.sample);
    if has_png_stack(data_png)
        fprintf('  - Skipped (PNG stack already exists)\n');
    else
        unlock_tomocube_dat(data_dat, data_png);
    end

    fprintf('[%d/%d] background: %s\n', pair_idx, num_pairs, current_pair.background);
    if has_png_stack(bg_png)
        fprintf('  - Skipped (PNG stack already exists)\n');
    else
        unlock_tomocube_dat(bg_dat, bg_png);
    end
end

function tf = has_png_stack(dir_path)
    if ~exist(dir_path, 'dir')
        tf = false;
        return;
    end
    png_files = dir(fullfile(dir_path, '*.png'));
    tf = ~isempty(png_files);
end
