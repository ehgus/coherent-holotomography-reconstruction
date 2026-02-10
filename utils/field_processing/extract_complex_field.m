function [field, optical_params, illum_k0] = extract_complex_field(background_stack, sample_stack, imaging_condition, field_generator_condition, tomography_generator_condition)
% extract_complex_field Process background and sample field pairs for field retrieval
%
% Syntax:
%   [field, updated_params, illum_k0] = extract_complex_field(background_stack, sample_stack, imaging_condition, field_generator_condition, tomography_generator_condition)
%
% Inputs:
%   background_stack - 3D array of background images [height, width, num_images]
%   sample_stack - 3D array of sample images [height, width, num_images]
%   imaging_condition - Structure containing imaging parameters:
%       .wavelength - Wavelength in microns
%       .NA - Numerical aperture
%       .RI_bg - Background refractive index
%       .resolution - [dx, dy] image resolution in microns
%   field_generator_condition - Structure containing field processing parameters:
%       .cutout_portion - Portion to cut out for centering (0 to 0.5)
%       .other_corner - Boolean to use other corner
%       .conjugate_field - Boolean to conjugate fields
%       .normalidx - Index of normal image
%   tomography_generator_condition - Structure containing tomography parameters:
%       .resolution - Target resolution for tomogram [dx, dy, dz]
%       .zsize_micron - Z-size in microns
%
% Outputs:
%   field - Processed sample field (resized to tomography resolution)
%   optical_params - Updated optical parameters
%   illum_k0 - Peak positions in Fourier space [2 x num_images]
%
% Description:
%   This function performs the core field retrieval processing based on
%   FIELD_EXPERIMENTAL_RETRIEVAL.m. It converts images to Fourier space,
%   centers them, resizes to match desired resolution, crops according to NA,
%   and converts back to real space with phase unwrapping.
%   Finally, fields are resized to match tomography generator resolution.
%
% See also: FIELD_EXPERIMENTAL_RETRIEVAL

    % Input validation
    assert(isequal(size(background_stack), size(sample_stack)), ...
        'Background and sample fields must be of same size');

    % Convert to single for processing
    input_field = single(background_stack);
    output_field = single(sample_stack);

    % Check for overexposure
    if isinteger(background_stack)
        maximum_value = single(intmax(class(background_stack)));
        is_overexposed = max(max(input_field(:)), max(output_field(:))) > maximum_value;
        if is_overexposed
            warning('Images are overexposed');
        end
    end

    % Step 1: Convert to Fourier space
    input_field = fft2(input_field);
    output_field = fft2(output_field);
    [xsize, ysize, zsize] = size(input_field);

    % Step 2: Center the field in Fourier space
    assert(1 <= field_generator_condition.normalidx && field_generator_condition.normalidx <= zsize, ...
        'Normal index should be a valid z index');

    search_band_1 = round(xsize*(1/2 - field_generator_condition.cutout_portion)):round(xsize/2);

    normal_bg = zeros(xsize, ysize);
    normal_bg(search_band_1, :) = input_field(search_band_1, :, field_generator_condition.normalidx);

    [~, linear_index] = max(abs(normal_bg(:)));
    [center_pos_1, center_pos_2] = ind2sub(size(normal_bg), linear_index);
    peak2origin = [1 - center_pos_1, 1 - center_pos_2, 0];

    input_field = circshift(input_field, peak2origin);
    output_field = circshift(output_field, peak2origin);

    % Step 3: Create NA circle and crop

    % Create frequency coordinates
    kmax = imaging_condition.NA / imaging_condition.wavelength ;
    dx = imaging_condition.resolution(1);
    dy = imaging_condition.resolution(2);

    fx = (-xsize/2:xsize/2-1) / (xsize * dx);
    fy = (-ysize/2:ysize/2-1) / (ysize * dy);
    [FX, FY] = meshgrid(fy, fx);
    coorxy = sqrt(FX.^2 + FY.^2);
    NA_circle = coorxy < kmax;

    shifted_NA_circle = ifftshift(NA_circle);
    input_field = input_field .* shifted_NA_circle;
    output_field = output_field .* shifted_NA_circle;

    % Find peaks for illum_k0
    illum_k0 = zeros(2, size(input_field, 3));
    [~, max_indices] = max(abs(input_field), [], [1, 2], 'linear');
    for jj = 1:size(input_field, 3)
        [y_pos, x_pos, ~] = ind2sub(size(input_field(:, :, jj)), max_indices(jj));
        illum_k0(1, jj) = mod(y_pos + floor(xsize/2) - 1, xsize) - floor(xsize/2);
        illum_k0(2, jj) = mod(x_pos + floor(ysize/2) - 1, ysize) - floor(ysize/2);
    end

    % Step 4.5: Resize fourier images to match tomography generator resolution
    % (before converting to real space)
    if nargin >= 5 && ~isempty(tomography_generator_condition)
        old_xsize = xsize;
        old_ysize = ysize;
        
        % Calculate target size based on tomography resolution and current FOV
        current_fov_x = xsize * imaging_condition.resolution(1);
        current_fov_y = ysize * imaging_condition.resolution(2);
        
        % Calculate new size based on tomography resolution
        xsize = 2 * round(current_fov_x / tomography_generator_condition.resolution(1) / 2);
        ysize = 2 * round(current_fov_y / tomography_generator_condition.resolution(2) / 2);
        
        % Perform padding/cropping if size changed
        if xsize ~= old_xsize || ysize ~= old_ysize
            old_field = {input_field, output_field};
            new_field = {zeros(xsize, ysize, zsize), zeros(xsize, ysize, zsize)};
            half_xsize = floor(min([old_xsize xsize])/2);
            half_ysize = floor(min([old_ysize ysize])/2);
            
            % Copy quadrants: top-left, top-right, bottom-left, bottom-right
            for i = 1:2
                new_field{i}(1:half_xsize, 1:half_ysize, :) = old_field{i}(1:half_xsize, 1:half_ysize, :);
                new_field{i}(1:half_xsize, end-half_ysize+1:end, :) = old_field{i}(1:half_xsize, end-half_ysize+1:end, :);
                new_field{i}(end-half_xsize+1:end, 1:half_ysize, :) = old_field{i}(end-half_xsize+1:end, 1:half_ysize, :);
                new_field{i}(end-half_xsize+1:end, end-half_ysize+1:end, :) = old_field{i}(end-half_xsize+1:end, end-half_ysize+1:end, :);
            end
            
            input_field = new_field{1};
            output_field = new_field{2};
            
            % Update imaging condition resolution to reflect actual resolution after resize
            imaging_condition.resolution(1) = current_fov_x / xsize;
            imaging_condition.resolution(2) = current_fov_y / ysize;
        end
    end

    % Step 5: Convert to real space
    input_field = ifft2(input_field);
    output_field = ifft2(output_field);

    % Step 6: Phase correction
    field = output_field ./ input_field;
    if field_generator_condition.conjugate_field
        field = conj(field);
    end

    % Subpixel phase shift correction
    phase = angle(field);
    % note: ad-hoc to work with current unwrapping function
    padded_phase = padarray(phase, [max(0,size(field,2)-size(field,1)), max(0,size(field,1)-size(field,2)), 0], ...
                        'replicate','post');
    slice_step = 10;
    for i = 1:ceil(size(padded_phase,3)/slice_step)
        zslice_view = (1+slice_step*(i-1)):min(slice_step*i,size(padded_phase,3));
        padded_phase(:,:,zslice_view) = gather(unwrapp2_gpu(gpuArray(single(padded_phase(:,:,zslice_view)))));
    end
    phase = padded_phase(1:size(field,1),1:size(field,2),:);
    for jj = 1:size(phase, 3)
        phase(:, :, jj) = shift_phi(phase(:, :, jj),1,1);
    end
    field = abs(field) .* exp(1i * phase);
    
    % Set output optical parameters
    optical_params = imaging_condition;
end
