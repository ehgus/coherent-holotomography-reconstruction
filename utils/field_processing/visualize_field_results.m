function visualize_field_results(input_field, output_field, slice_idx)
% VISUALIZE_FIELD_RESULTS Display amplitude and phase of field retrieval results
%
% Syntax:
%   visualize_field_results(input_field, output_field, slice_idx)
%
% Inputs:
%   input_field - Processed background field [height x width x 1 x num_images]
%   output_field - Processed sample field [height x width x 1 x num_images]
%   slice_idx - (Optional) Index of slice to visualize. If not provided, shows middle slice.
%
% Description:
%   Creates a figure showing the amplitude and phase of the field ratio
%   (output_field / input_field) for a specified slice.
%
% See also: process_field_pair, save_field_results

    if nargin < 3
        slice_idx = round(size(input_field, 4) / 2);
    end

    % Validate slice index
    num_slices = size(input_field, 4);
    if slice_idx < 1 || slice_idx > num_slices
        error('Slice index must be between 1 and %d', num_slices);
    end

    % Extract the specified slice
    input_slice = squeeze(input_field(:, :, 1, slice_idx));
    output_slice = squeeze(output_field(:, :, 1, slice_idx));

    % Calculate amplitude and phase
    field_ratio = output_slice ./ input_slice;
    amplitude = abs(field_ratio);
    phase = angle(field_ratio);

    % Create visualization
    figure('Name', sprintf('Field Results - Slice %d/%d', slice_idx, num_slices), 'NumberTitle', 'off');

    subplot(2, 2, 1);
    imagesc(abs(input_slice));
    axis image;
    colorbar;
    title('Input Field Amplitude');
    colormap(gca, 'gray');

    subplot(2, 2, 2);
    imagesc(abs(output_slice));
    axis image;
    colorbar;
    title('Output Field Amplitude');
    colormap(gca, 'gray');

    subplot(2, 2, 3);
    imagesc(amplitude);
    axis image;
    colorbar;
    title('Ratio Amplitude');
    colormap(gca, 'gray');

    subplot(2, 2, 4);
    imagesc(phase);
    axis image;
    colorbar;
    title('Ratio Phase');
    colormap(gca, 'parula');

    % Add overall title
    sgtitle(sprintf('Field Retrieval Results - Slice %d/%d', slice_idx, num_slices));
end
