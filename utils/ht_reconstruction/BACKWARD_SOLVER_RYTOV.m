classdef BACKWARD_SOLVER_RYTOV < handle
    properties (SetAccess = protected, Hidden = true)
        parameters = struct(...
            ... % device information 
            'wavelength', 0.532, ...
            'NA',1.2, ...
            'RI_bg',1.336, ...
            ... % calculation option
            'size',[100 100 100], ...
            'resolution',[0.1 0.1 0.1], ...
            'use_GPU',true ...
        );
        utility;
    end
    methods
        function h=BACKWARD_SOLVER_RYTOV(params)
            if nargin==1
                warning ('off','all');
                h.parameters=update_struct(h.parameters, params);
                warning ('on','all');
            end
        end
        function [RI, ORytov]=solve(h,output_field,illum_k0)
            warning('off','all');
            h.utility=DERIVE_OPTICAL_TOOL(h.parameters);
            warning('on','all');
            % check fields and parameters
            assert(ndims(output_field) == 3, 'You need to provide the field with 3 dimenssion : dim1 x dim2 x illuminationnumber')
            % ISSUE: mismatch size between output_field and RI -> resize output_field

            %preset variables
            kx_res = h.utility.fourier_space.res{1};
            ky_res = h.utility.fourier_space.res{2};
            kz_res = h.utility.fourier_space.res{3};
            xsize = h.parameters.size(1);
            ysize = h.parameters.size(2);
            zsize = h.parameters.size(3);

            %find angle
            f_dx = illum_k0(1,:);
            f_dy = illum_k0(2,:);
            f_dz=round(real(sqrt((h.utility.k0_nm)^2-(f_dx*kx_res).^2-(f_dy*ky_res).^2))/kz_res);

            NA_circle = ifftshift(h.utility.NA_circle);
            xind=find(NA_circle);
            kz=ifftshift(reshape(h.utility.k3,xsize,ysize));
            fx=[0:floor((xsize-1)/2) -floor((xsize)/2):-1];
            fy=[0:floor((ysize-1)/2) -floor((ysize)/2):-1];
            fz=round(kz/kz_res);
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
            ORytov(Count>0)=ORytov(Count>0)./Count(Count>0)/kz_res; % should be (um^-2)*(px*py*pz), so (px*py*pz/um^3) should be multiplied.
            potential=gather(fftshift(ifftn(ORytov),3));
            RI = potential2RI(potential*4*pi,h.parameters.wavelength,h.parameters.RI_bg);
        end
    end
end