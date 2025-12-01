function save_field_results(input_field, output_field, output_path, sample_name, updated_params)
% SAVE_FIELD_RESULTS Save field retrieval results to disk
%
% Syntax:
%   save_field_results(input_field, output_field, output_path, sample_name, updated_params)
%
% Inputs:
%   input_field - Processed background field [height x width x 1 x num_images]
%   output_field - Processed sample field [height x width x 1 x num_images]
%   output_path - Base path for saving results
%   sample_name - Name identifier for this sample
%   updated_params - Structure with updated optical parameters
%
% Outputs:
%   Saves the following files in output_path/sample_name/:
%       - input_field.mat: Background field data
%       - output_field.mat: Sample field data
%       - parameters.mat: Updated optical parameters
%       - amplitude.mat: Amplitude (abs(output_field/input_field))
%       - phase.mat: Phase (angle(output_field/input_field))
%
% See also: process_field_pair

    % Create output directory
    sample_output_path = fullfile(output_path, sample_name);
    if ~exist(sample_output_path, 'dir')
        mkdir(sample_output_path);
    end

    % Save field data
    save(fullfile(sample_output_path, 'input_field.mat'), 'input_field', '-v7.3');
    save(fullfile(sample_output_path, 'output_field.mat'), 'output_field', '-v7.3');

    % Save parameters
    save(fullfile(sample_output_path, 'parameters.mat'), 'updated_params');

    % Calculate and save amplitude and phase
    field_ratio = output_field ./ input_field;
    amplitude = abs(field_ratio);
    phase = angle(field_ratio);

    save(fullfile(sample_output_path, 'amplitude.mat'), 'amplitude', '-v7.3');
    save(fullfile(sample_output_path, 'phase.mat'), 'phase', '-v7.3');
end
