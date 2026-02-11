function field_fourier = refocus_field_fourier_domain(field_fourier, z_shift, imaging_condition)
% refocus_field_fourier_domain Refocus a complex field by a given axial shift using
% angular spectrum method.
%   refocused_field = refocus_field_fourier_domain(field, z_shift, imaging_condition)
%   refocuses the input complex field 'field' by 'z_shift' micrometers.
    k0_nm = imaging_condition.RI_immersion / imaging_condition.wavelength;
    resolution = imaging_condition.resolution;
    
    Nx = size(field_fourier, 1);
    Ny = size(field_fourier, 2);

    % Frequency coordinates
    fx = [0:floor((Nx-1)/2) -floor(Nx/2):-1]' / (Nx * resolution(1));
    fy = [0:floor((Ny-1)/2) -floor(Ny/2):-1]  / (Ny * resolution(2));

    % Transfer function for propagation
    H = exp(2i*pi * z_shift * sqrt(k0_nm^2 - fx.^2 - fy.^2));
    H(fx.^2 + fy.^2 > k0_nm^2) = 0; % Evanescent waves

    % Perform angular spectrum propagation
    field_fourier = field_fourier .* H;
end