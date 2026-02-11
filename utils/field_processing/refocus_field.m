function refocused_field = refocus_field(field, z_shift, imaging_condition)
% REFOCUS_FIELD Refocus a complex field by a given axial shift using
% angular spectrum method.
%   refocused_field = REFOCUS_FIELD(field, z_shift, imaging_condition)
%   refocuses the input complex field 'field' by 'z_shift' micrometers.
    % Perform angular spectrum propagation
    field_F = fft2(field);
    refocused_field_F = refocus_field_fourier_domain(field_F, z_shift, imaging_condition);
    refocused_field = ifft2(refocused_field_F);
end