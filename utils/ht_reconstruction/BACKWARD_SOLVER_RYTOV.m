classdef BACKWARD_SOLVER_RYTOV < handle
    properties
        wavelength          % wavelength in micron
        NA                  % numerical aperture
        RI_bg               % background refractive index
        field_resolution    % field image resolution in micron [dx, dy]
        tomogram_size       % reconstruction volume size [nx, ny, nz]
        tomogram_resolution % reconstruction resolution in micron [dx, dy, dz]
    end
    methods
        function h=BACKWARD_SOLVER_RYTOV(params)
            h.wavelength = params.wavelength;
            h.NA = params.NA;
            h.RI_bg = params.RI_bg;
            h.field_resolution = params.field_resolution;
            h.tomogram_size = params.tomogram_size;
            h.tomogram_resolution = params.tomogram_resolution;
        end
        function [potential, Count]=solve(h,output_field,illum_k0)
            % check fields and parameters
            assert(ndims(output_field) == 3, 'You need to provide the field with 3 dimenssion : dim1 x dim2 x illuminationnumber')

            % Common information
            k_res = 1 ./ (h.tomogram_resolution .* h.tomogram_size);
            k0_nm = h.RI_bg / h.wavelength;
            kmax = h.NA / h.wavelength;
            illum_k0_3d = zeros(3,size(illum_k0,2),'single');
            illum_k0_3d(1:2,:) = illum_k0;
            illum_k0_3d(3,:) = round(real(sqrt((k0_nm)^2 - (illum_k0(1,:)*k_res(1)).^2 - (illum_k0(2,:)*k_res(2)).^2))/k_res(3));

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
            kz=kz(valid_2d_indices);

            % tomogram-related information
            xsize_3d = h.tomogram_size(1);
            ysize_3d = h.tomogram_size(2);
            zsize_3d = h.tomogram_size(3);
            % 3D Fourier space indices
            [xcoords_3d, ycoords_3d] = ind2sub([xsize_2d, ysize_2d], valid_2d_indices);
            xcoords_3d(xcoords_3d > xsize_2d/2) = xcoords_3d(xcoords_3d > xsize_2d/2) + (xsize_3d - xsize_2d);
            ycoords_3d(ycoords_3d > ysize_2d/2) = ycoords_3d(ycoords_3d > ysize_2d/2) + (ysize_3d - ysize_2d);
            zcoords_3d=round(kz/k_res(3));

            % Map 2D field to 3D Fourier space
            potential=zeros(xsize_3d,ysize_3d,zsize_3d,'single');
            Count=zeros(xsize_3d,ysize_3d,zsize_3d,'single');
            % Extract rytov field
            phase = angle(output_field);
            % note: ad-hoc to work with current unwrapping function
            padded_phase = padarray(phase, [max(0,size(output_field,2)-size(output_field,1)), max(0,size(output_field,1)-size(output_field,2)), 0], ...
                                'replicate','post');
            slice_step = 10;
            for i = 1:ceil(size(padded_phase,3)/slice_step)
                zslice_view = (1+slice_step*(i-1)):min(slice_step*i,size(padded_phase,3));
                padded_phase(:,:,zslice_view) = gather(unwrapp2_gpu(gpuArray(single(padded_phase(:,:,zslice_view)))));
            end
            phase = padded_phase(1:size(output_field,1),1:size(output_field,2),:);
            amp = abs(output_field);
            UsRytov=log(amp)+1i*phase;
            UsRytov=gather(fft2(gpuArray(UsRytov)));
            for i = 1:size(output_field,3)
                % Acquire 2D Fourier space indices
                % Rescale field
                UsRytov_slice=circshift(UsRytov(:,:,i),[illum_k0_3d(1:2,i)]);
                
                UsRytov_slice=kz/1i.*UsRytov_slice(valid_2d_indices);% unit: (um^1) % kz is spatial frequency, so 2pi is multiplied for wave vector
                % Acquire 3D Fourier space indices
                Fx_3d=xcoords_3d-illum_k0_3d(1,i);
                Fy_3d=ycoords_3d-illum_k0_3d(2,i);
                Fz_3d=zcoords_3d-illum_k0_3d(3,i);
                Fx_3d=mod(Fx_3d,xsize_3d)+1;
                Fy_3d=mod(Fy_3d,ysize_3d)+1;
                Fz_3d=mod(Fz_3d,zsize_3d)+1;
                Kzp_3d=sub2ind(size(Count),Fx_3d,Fy_3d,Fz_3d);
                % Accumulate into 3D Fourier space
                potential(Kzp_3d)=potential(Kzp_3d)+UsRytov_slice;
                Count(Kzp_3d)=Count(Kzp_3d)+1;
            end
            potential(Count>0)=potential(Count>0)./Count(Count>0)/k_res(3); % should be (um^-2)*(px*py*pz), so (px*py*pz/um^3) should be multiplied.
            potential=fftshift(gather(ifftn(gpuArray(potential))),3);
        end
    end
end