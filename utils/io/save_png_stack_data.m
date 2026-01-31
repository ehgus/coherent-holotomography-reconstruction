function save_png_stack_data(stack, output_dir)
% SAVE_PNG_STACK_DATA Save a 3D image stack as a series of PNG files
%
% Syntax:
%   save_png_stack_data(stack, output_dir)
%
% Inputs:
%   stack      - 3D uint8 array with dimensions [height, width, num_images]
%   output_dir - String or char array specifying the output directory path
%                where PNG files will be saved
%
% Outputs:
%   Creates PNG files named 'image-0001.png', 'image-0002.png', etc.
%   in the specified output directory
%
% Description:
%   This function takes a 3D image stack and saves each slice as a
%   separate PNG file with zero-padded sequential numbering.
%
% Example:
%   stack = load_bin_data('F:\Data\TGV imaging\sample_IMAGING_001');
%   output_path = 'F:\Data\TGV imaging\sample_IMAGING_001\data3d';
%   save_png_stack_data(stack, output_path);
%
% See also: load_png_stack_data, unwrap_bin_to_png, imwrite

    % Validate input stack
    if ~isa(stack, 'uint8')
        error('The image stack should be uint8');
    end

    if ndims(stack) ~= 3
        error('Unexpected number of dimensions in stack: %d', ndims(stack));
    end

    % Create output directory if it doesn't exist
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    % Get the number of images in the stack
    % Assuming the last dimension is the image index
    num_images = size(stack, ndims(stack));

    % Save each slice as a PNG file
    fprintf('Saving %d images to %s\n', num_images, output_dir);

    for i = 1:num_images
        img = stack(:, :, i);

        % Generate filename with zero-padded numbering
        filename = sprintf('image-%04d.png', i);
        filepath = fullfile(output_dir, filename);

        % Write the image
        imwrite(img, filepath);
    end

    fprintf('Successfully saved all images to: %s\n', output_dir);
end
