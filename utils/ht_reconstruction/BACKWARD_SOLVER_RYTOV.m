classdef BACKWARD_SOLVER_RYTOV < handle
    properties
        wavelength  % wavelength in micron
        NA          % numerical aperture
        RI_bg       % background refractive index
        volume_size % reconstruction volume size [nx, ny, nz]
        resolution  % reconstruction resolution in micron [dx, dy, dz]
        use_GPU     % boolean to use GPU acceleration
    end
    methods
        function h=BACKWARD_SOLVER_RYTOV(params)
            h.wavelength = params.wavelength;
            h.NA = params.NA;
            h.RI_bg = params.RI_bg;
            h.volume_size = params.volume_size;
            h.resolution = params.resolution;
            h.use_GPU = params.use_GPU;
        end
        function [RI, ORytov]=solve(h,output_field,illum_k0)
            % check fields and parameters
            assert(ndims(output_field) == 3, 'You need to provide the field with 3 dimenssion : dim1 x dim2 x illuminationnumber')
            % ISSUE: mismatch size between output_field and RI -> resize output_field

            % Extract parameters
            xsize = h.volume_size(1);
            ysize = h.volume_size(2);
            zsize = h.volume_size(3);

            % Calculate fourier space resolution
            k_res = 1 ./ (h.resolution .* h.volume_size);

            % Calculate wavenumber parameters
            k0_nm = h.RI_bg / h.wavelength;
            kmax = h.NA / h.wavelength;

            % Calculate fourier space coordinates for NA circle
            fx_coords = single((1:xsize) - (floor(xsize/2) + 1)) * k_res(1);
            fy_coords = single((1:ysize) - (floor(ysize/2) + 1)) * k_res(2);
            [FX, FY] = meshgrid(fy_coords, fx_coords);
            coorxy = sqrt(FX.^2 + FY.^2);
            NA_circle = coorxy < kmax;

            % Calculate k3 (z-component of wavenumber)
            k3 = (k0_nm).^2 - (coorxy).^2;
            k3(k3 < 0) = 0;
            k3 = sqrt(k3);

            %find angle
            f_dx = illum_k0(1,:);
            f_dy = illum_k0(2,:);
            f_dz=round(real(sqrt((k0_nm)^2-(f_dx*k_res(1)).^2-(f_dy*k_res(2)).^2))/k_res(3));

            NA_circle_shift = ifftshift(NA_circle);
            xind=find(NA_circle_shift);
            kz=ifftshift(reshape(k3,xsize,ysize));
            fx=[0:floor((xsize-1)/2) -floor((xsize)/2):-1];
            fy=[0:floor((ysize-1)/2) -floor((ysize)/2):-1];
            fz=round(kz/k_res(3));
            fx=fx(rem(xind-1,xsize)+1)';
            fy=fy(floor((xind-1)/xsize)+1)';
            fz=fz(xind);
            kz=kz(xind);

            ORytov=gpuArray(zeros(xsize,ysize,zsize,'single'));
            Count=gpuArray(zeros(xsize,ysize,zsize,'single')); 
            for i = 1:size(output_field,3)
                % Extract rytov field
                phase = unwrap_phase(angle(output_field(:,:,i)));
                amp = abs(output_field(:,:,i));
                UsRytov=squeeze(log(amp)+1i*phase);
                UsRytov=fft2(UsRytov);
                % Rescale field
                UsRytov=circshift(UsRytov,[f_dx(i) f_dy(i)]);
                Fx=f_dx(i)-fx;Fy=f_dy(i)-fy;Fz=f_dz(i)-fz;
                Uprime=kz/1i.*UsRytov(xind);% unit: (um^1) % kz is spatial frequency, so 2pi is multiplied for wave vector
                
                Fx=mod(Fx,xsize)+1;
                Fy=mod(Fy,ysize)+1;
                Fz=mod(Fz,zsize)+1;
                Kzp=sub2ind(size(Count),Fx,Fy,Fz);

                ORytov(Kzp)=ORytov(Kzp)+Uprime;
                Count(Kzp)=Count(Kzp)+(Uprime~=0);
            end
            ORytov(Count>0)=ORytov(Count>0)./Count(Count>0)/k_res(3); % should be (um^-2)*(px*py*pz), so (px*py*pz/um^3) should be multiplied.
            potential=gather(fftshift(ifftn(ORytov),3));
            RI = potential2RI(potential*4*pi,h.wavelength,h.RI_bg);
        end
    end
end