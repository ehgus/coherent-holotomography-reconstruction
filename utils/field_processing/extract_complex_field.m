function [field, optical_params, illum_k0] = extract_complex_field(background_stack, sample_stack, optical_params, processing_params)
% extract_complex_field Process background and sample field pairs for field retrieval
%
% Syntax:
%   [field, updated_params, illum_k0] = extract_complex_field(background_stack, sample_stack, optical_params, processing_params)
%
% Inputs:
%   background_stack - 3D array of background images [height, width, num_images]
%   sample_stack - 3D array of sample images [height, width, num_images]
%   optical_params - Structure containing optical parameters:
%       .wavelength - Wavelength in microns
%       .NA - Numerical aperture
%       .RI_bg - Background refractive index
%       .resolution - [dx, dy, dz] spatial resolution
%       .resolution_image - [dx, dy] image resolution
%       .vector_simulation - Boolean for vector field simulation
%       .use_abbe_sine - Boolean for Abbe sine correction
%       .use_abbe_correction - Boolean for Abbe correction
%   processing_params - Structure containing processing parameters:
%       .cutout_portion - Portion to cut out for centering (0 to 0.5)
%       .other_corner - Boolean to use other corner
%       .conjugate_field - Boolean to conjugate fields
%       .normalidx - Index of normal image
%       .use_GPU - Boolean to use GPU acceleration
%
% Outputs:
%   field - Processed sample field
%   updated_params - Updated optical parameters
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
    assert(optical_params.resolution_image(1) == optical_params.resolution_image(2), ...
        'Image resolution must be isotropic');
    assert(optical_params.resolution(1) == optical_params.resolution(2), ...
        'Output resolution must be isotropic');
    assert(0 < processing_params.cutout_portion && processing_params.cutout_portion < 0.5, ...
        'Cutout portion should be in (0, 0.5)');

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
    assert(1 <= processing_params.normalidx && processing_params.normalidx <= zsize, ...
        'Normal index should be a valid z index');

    search_band_1 = round(xsize*(1/2 - processing_params.cutout_portion)):round(xsize/2);
    if processing_params.other_corner
        search_band_1 = round(xsize/2):round(xsize*(1/2 + processing_params.cutout_portion));
    end
    search_band_2 = round(ysize*(1/2 - processing_params.cutout_portion)):round(ysize/2);

    normal_bg = zeros(xsize, round(ysize/2));
    normal_bg(search_band_1, search_band_2) = input_field(search_band_1, search_band_2, processing_params.normalidx);

    [~, linear_index] = max(abs(normal_bg(:)));
    [center_pos_1, center_pos_2] = ind2sub(size(normal_bg), linear_index);
    peak2origin = [1 - center_pos_1, 1 - center_pos_2, 0];

    input_field = circshift(input_field, peak2origin);
    output_field = circshift(output_field, peak2origin);

    % Step 3: Resize to match desired resolution
    old_xsize = xsize;
    old_ysize = ysize;
    resolution_ratio = optical_params.resolution(1:2) ./ optical_params.resolution_image(1:2);
    xsize = 2 * round(old_xsize / resolution_ratio(1) / 2);
    ysize = 2 * round(old_ysize / resolution_ratio(2) / 2);

    if xsize ~= old_xsize || ysize ~= old_ysize
        old_field = {input_field, output_field};
        new_field = {zeros(xsize, ysize, zsize), zeros(xsize, ysize, zsize)};
        half_xsize = floor(min([old_xsize, xsize]) / 2);
        half_ysize = floor(min([old_ysize, ysize]) / 2);

        for i = 1:2
            new_field{i}(1:half_xsize, 1:half_ysize, :) = old_field{i}(1:half_xsize, 1:half_ysize, :);
            new_field{i}(1:half_xsize, end-half_ysize+1:end, :) = old_field{i}(1:half_xsize, end-half_ysize+1:end, :);
            new_field{i}(end-half_xsize+1:end, 1:half_ysize, :) = old_field{i}(end-half_xsize+1:end, 1:half_ysize, :);
            new_field{i}(end-half_xsize+1:end, end-half_ysize+1:end, :) = old_field{i}(end-half_xsize+1:end, end-half_ysize+1:end, :);
        end

        clear old_field

        input_field = new_field{1};
        output_field = new_field{2};
        optical_params.resolution(1) = optical_params.resolution_image(1) * old_xsize / xsize;
        optical_params.resolution(2) = optical_params.resolution_image(2) * old_ysize / ysize;
    end

    % Step 4: Create NA circle and crop
    optical_params.size = [xsize, ysize, zsize];

    % Create frequency coordinates
    kmax = optical_params.NA / optical_params.wavelength;
    dx = optical_params.resolution(1);
    dy = optical_params.resolution(2);

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
    for jj = 1:size(input_field, 3)
        shifted_input_field = fftshift(input_field(:, :, jj));
        [~, max_idx] = max(abs(shifted_input_field(:)));
        [y_pos, x_pos] = ind2sub([xsize, ysize], max_idx);
        illum_k0(1, jj) = y_pos - floor(xsize/2) - 1;
        illum_k0(2, jj) = x_pos - floor(ysize/2) - 1;
    end

    % Step 5: Convert to real space
    input_field = ifft2(input_field);
    output_field = ifft2(output_field);

    % Step 6: Phase correction
    field = output_field ./ input_field;
    if processing_params.conjugate_field
        field = conj(field);
    end

    % Crop edges
    field = field(3:(end-2), 3:(end-2), :);
    optical_params.size(1) = size(field, 1);
    optical_params.size(2) = size(field, 2);

    % Subpixel phase shift correction
    for jj = 1:size(field, 3)
        field(:, :, jj) = remove_abs_phase(field(:, :, jj));
    end

    % Update parameters for output
    optical_params.use_GPU = processing_params.use_GPU;
end

function complex_phase = remove_abs_phase(complex_phase)
    % Consider background position (negligible angle change)
    illum_angle_map_x = angle(circshift(complex_phase,1,1)./complex_phase);
    illum_angle_map_y = angle(circshift(complex_phase,1,2)./complex_phase);

    % Correct absolute phase
    naive_bg_mask = abs(illum_angle_map_x) < std(illum_angle_map_x,1,'all') & ...
                    abs(illum_angle_map_y) < std(illum_angle_map_y,1,'all');
    abs_angle = angle(mean(complex_phase(naive_bg_mask),'all'));
    complex_phase = complex_phase .* exp(-1i * abs_angle);
end
