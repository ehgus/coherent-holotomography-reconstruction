function stack = load_png_stack_data(base_path)
% LOAD_PNG_STACK_DATA Load PNG files from a directory into a 3D stack
%
% Syntax:
%   stack = load_png_stack_data(base_path)
%
% Inputs:
%   base_path - String or char array specifying the directory path
%               containing PNG files (e.g., 'data3d' subdirectory created
%               by unwrap_bin_to_png)
%
% Outputs:
%   stack - 3D uint8 array containing the loaded image stack with
%           dimensions [height, width, num_images]
%
% Description:
%   This function reads all PNG files from the specified directory and
%   combines them into a 3D stack. When used on a 'data3d' directory
%   created by unwrap_bin_to_png, it should return the same stack as
%   load_bin_data would return from the parent directory's bin files.
%
% Example:
%   data_path = 'F:\Data\TGV imaging\sample_IMAGING_001\data3d';
%   images = load_png_stack_data(data_path);
%
% See also: load_bin_data, unwrap_bin_to_png, imread

    % Get list of PNG files in the directory
    png_files = dir(fullfile(base_path, '*.png'));

    if isempty(png_files)
        error('No PNG files found in directory: %s', base_path);
    end

    % Sort files by name to ensure correct order
    [~, idx] = sort({png_files.name});
    png_files = png_files(idx);

    num_images = length(png_files);

    % Read first image to get dimensions
    first_img = imread(fullfile(base_path, png_files(1).name));
    [height, width, channels] = size(first_img);

    % Check if image is grayscale
    if channels ~= 1
        error('Expected grayscale PNG images, but found %d channels', channels);
    end

    % Check data type
    if ~isa(first_img, 'uint8')
        error('Expected uint8 PNG images, but found %s', class(first_img));
    end

    % Pre-allocate 3D stack
    stack = zeros(height, width, num_images, 'uint8');

    % Load first image
    stack(:, :, 1) = first_img;

    for i = 2:num_images
        img = imread(fullfile(base_path, png_files(i).name));

        % Verify dimensions match
        if ~isequal(size(img), [height, width])
            error('Image %s has inconsistent dimensions', png_files(i).name);
        end

        stack(:, :, i) = img;
    end
end
