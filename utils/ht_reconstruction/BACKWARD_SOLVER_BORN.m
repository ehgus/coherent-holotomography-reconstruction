classdef BACKWARD_SOLVER_BORN < handle
    properties
        wavelength          % wavelength in micron
        NA                  % numerical aperture
        RI_bg               % background refractive index
        field_resolution    % field image resolution in micron [dx, dy]
        tomogram_size       % reconstruction volume size [nx, ny, nz]
        tomogram_resolution % reconstruction resolution in micron [dx, dy, dz]
    end
    methods
        function h=BACKWARD_SOLVER_BORN(params)
            h.wavelength = params.wavelength;
            h.NA = params.NA;
            h.RI_bg = params.RI_bg;
            h.field_resolution = params.field_resolution;
            h.tomogram_size = params.tomogram_size;
            h.tomogram_resolution = params.tomogram_resolution;
        end
        
        function field_fft = extractField(h, output_field)
            % Extract Born field (can be overridden in subclasses)
            field_fft = gather(fft2(gpuArray(output_field)));
        end
        
        function [potential, Count]=solve(h,output_field,illum_k0_pixel)
            % check fields and parameters
            assert(ndims(output_field) == 3, 'You need to provide the field with 3 dimenssion : dim1 x dim2 x illuminationnumber')

            % Common information
            k_res = 1 ./ (h.tomogram_resolution .* h.tomogram_size);
            k0_nm = h.RI_bg / h.wavelength;
            kmax = h.NA / h.wavelength;
            illum_k0_3d = zeros(3,size(illum_k0_pixel,2),'single');
            illum_k0_3d(1:2,:) = illum_k0_pixel .* reshape(k_res(1:2),2,1);
            illum_k0_3d(3,:) = real(sqrt((k0_nm)^2 - illum_k0_3d(1,:).^2 - illum_k0_3d(2,:).^2));

            % field-related information
            xsize_2d = size(output_field,1);
            ysize_2d = size(output_field,2);
            % Calculate fourier space coordinates for NA circle
            xcoords_2d = [0:floor((xsize_2d-1)/2) -floor(xsize_2d/2):-1]';
            ycoords_2d = [0:floor((ysize_2d-1)/2) -floor(ysize_2d/2):-1];
            kxy_coords_2d = sqrt((xcoords_2d *k_res(1)).^2 + (ycoords_2d * k_res(2)).^2);
            NA_circle_2d = kxy_coords_2d < kmax;
            valid_2d_indices = find(NA_circle_2d);
            % Calculate kz (z-component of wavenumber)
            kz = sqrt(max(0,(k0_nm).^2 - (kxy_coords_2d).^2));

            % tomogram-related information
            xsize_3d = h.tomogram_size(1);
            ysize_3d = h.tomogram_size(2);
            zsize_3d = h.tomogram_size(3);
            % 3D Fourier space indices
            [xcoords_3d_pixel, ycoords_3d_pixel] = ind2sub([xsize_2d, ysize_2d], valid_2d_indices);
            xcoords_3d_pixel(xcoords_3d_pixel > xsize_2d/2) = xcoords_3d_pixel(xcoords_3d_pixel > xsize_2d/2) + (xsize_3d - xsize_2d);
            ycoords_3d_pixel(ycoords_3d_pixel > ysize_2d/2) = ycoords_3d_pixel(ycoords_3d_pixel > ysize_2d/2) + (ysize_3d - ysize_2d);
            kxcoords_3d = (xcoords_3d_pixel - 1) * k_res(1);
            kycoords_3d = (ycoords_3d_pixel - 1) * k_res(2);
            kzcoords_3d = kz(valid_2d_indices);

            % Map 2D field to 3D Fourier space
            potential=zeros(xsize_3d,ysize_3d,zsize_3d,'single');
            Count=zeros(xsize_3d,ysize_3d,zsize_3d,'single');
            % Extract field (delegated to subclass if needed)
            field_fft = h.extractField(output_field);
            for i = 1:size(output_field,3)
                % Acquire 2D Fourier space indices
                % Rescale field
                field_slice=circshift(field_fft(:,:,i),[illum_k0_pixel(1:2,i)]);
                field_slice=kzcoords_3d/1i.*field_slice(valid_2d_indices);% unit: (um^1) % kz is spatial frequency, so 2pi is multiplied for wave vector
                % Acquire 3D Fourier space indices
                Fx_3d=mod(round((kxcoords_3d-illum_k0_3d(1,i))/k_res(1)),xsize_3d)+1;
                Fy_3d=mod(round((kycoords_3d-illum_k0_3d(2,i))/k_res(2)),ysize_3d)+1;
                Fz_3d=mod(round((kzcoords_3d-illum_k0_3d(3,i))/k_res(3)),zsize_3d)+1;
                Kzp_3d=sub2ind(size(Count),Fx_3d,Fy_3d,Fz_3d);
                % Accumulate into 3D Fourier space
                potential(Kzp_3d)=potential(Kzp_3d)+field_slice;
                Count(Kzp_3d)=Count(Kzp_3d)+1;
            end
            potential(Count>0)=potential(Count>0)./Count(Count>0)*h.tomogram_resolution(3); % should be (um^-2)*(px*py*pz), so (px*py*pz/um^3) should be multiplied.
            potential=fftshift(gather(ifftn(gpuArray(potential))),3);
        end
    end
end