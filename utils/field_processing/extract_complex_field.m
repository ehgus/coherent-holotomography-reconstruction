function [field, optical_params, illum_k0] = extract_complex_field(background_stack, sample_stack, imaging_condition, field_generator_condition)
% extract_complex_field Process background and sample field pairs for field retrieval
%
% Syntax:
%   [field, updated_params, illum_k0] = extract_complex_field(background_stack, sample_stack, imaging_condition, field_generator_condition)
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
%
% Outputs:
%   field - Processed sample field
%   optical_params - Updated optical parameters
%   illum_k0 - Peak positions in Fourier space [2 x num_images]
%
% Description:
%   This function performs the core field retrieval processing based on
%   FIELD_EXPERIMENTAL_RETRIEVAL.m. It converts images to Fourier space,
%   centers them, resizes to match desired resolution, crops according to NA,
%   and converts back to real space with phase unwrapping.
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
    optical_params = imaging_condition;

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

    % Step 5: Convert to real space
    input_field = ifft2(input_field);
    output_field = ifft2(output_field);

    % Step 6: Phase correction
    field = output_field ./ input_field;
    if field_generator_condition.conjugate_field
        field = conj(field);
    end

    % Subpixel phase shift correction
    for jj = 1:size(field, 3)
        field(:, :, jj) = remove_abs_phase(field(:, :, jj));
    end
end

function complex_phase = remove_abs_phase(complex_phase)
    % Consider background position (negligible angle change)
    illum_angle_map_x = angle(circshift(complex_phase,1,1)./complex_phase);
    illum_angle_map_y = angle(circshift(complex_phase,1,2)./complex_phase);

    % Correct absolute phase
    naive_bg_mask = abs(illum_angle_map_x) < std(illum_angle_map_x,1,'all')/3 & ...
                    abs(illum_angle_map_y) < std(illum_angle_map_y,1,'all')/3 & ...
                    abs(complex_phase) < 1.2 & ...
                    abs(complex_phase) > 0.8;
    abs_angle = zeros(1,1,size(complex_phase,3));
    for idx = 1:length(abs_angle)
        abs_angle(idx) = angle(mean(complex_phase(naive_bg_mask(:,:,idx))));
    end
    complex_phase = complex_phase .* exp(-1i * abs_angle);
end
