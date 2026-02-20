function output_volume = rotate3D(volume, angle_deg, axis, varargin)
% ROTATE3D  3D rotation using imrotate for 2D slices
%
% Usage:
%   output_volume = rotate3D(volume, angle_deg, axis)
%   output_volume = rotate3D(volume, angle_deg, axis, 'crop')
%
% Inputs:
%   volume      - 3D array (real-valued image/volume)
%   angle_deg   - Rotation angle in degrees (positive is counter-clockwise)
%   axis        - Rotation axis: 'x', 'y', or 'z' only
%   crop        - (Optional) 'crop' to remove padding, 'loose' for full size (default)
%
% Outputs:
%   output_volume - 3D rotated volume
%
% Description:
%   This function rotates a 3D volume by applying 2D rotations (imrotate) to 2D slices
%   perpendicular to the specified rotation axis.
%   
%   Rotation axes:
%   - 'x': Rotate around X-axis (rotate YZ slices)
%   - 'y': Rotate around Y-axis (rotate XZ slices)
%   - 'z': Rotate around Z-axis (rotate XY slices)

    % Parse inputs
    if nargin < 3
        error('At least 3 inputs required: volume, angle_deg, axis');
    end
    
    crop_flag = 'loose';  % default
    if nargin > 3
        crop_flag = varargin{1};
    end
    
    % Validate inputs
    if ~ismatrix(volume) && ndims(volume) ~= 3
        error('Volume must be a 2D or 3D array');
    end
    
    % Convert axis to lowercase and validate
    axis = lower(axis);
    valid_axes = {'x', 'y', 'z'};
    if ~ismember(axis, valid_axes)
        error('axis must be one of: x, y, z');
    end
    
    % Ensure 3D
    if ndims(volume) == 2
        volume = reshape(volume, [size(volume), 1]);
    end
    
    % Get volume dimensions
    [nx, ny, nz] = size(volume);
    
    % Perform rotation based on axis
    switch axis
        case 'z'
            % Rotate around Z-axis: rotate each XY slice
            output_volume = zeros(size(volume), 'like', volume);
            for k = 1:nz
                output_volume(:, :, k) = imrotate(volume(:, :, k), angle_deg, crop_flag);
            end
            
        case 'y'
            % Rotate around Y-axis: rotate each XZ slice
            output_volume = zeros([nx, nz, ny], 'like', volume);
            for j = 1:ny
                % Extract XZ slice
                slice = squeeze(volume(:, j, :));  % size: [nx, nz]
                rotated_slice = imrotate(slice, angle_deg, crop_flag);
                % Check if size changed
                [new_nx, new_nz] = size(rotated_slice);
                if new_nx ~= nx || new_nz ~= nz
                    % Resize to original size
                    rotated_slice = imresize(rotated_slice, [nx, nz]);
                end
                output_volume(:, j, :) = rotated_slice;
            end
            % Permute back to [nx, ny, nz]
            output_volume = permute(output_volume, [1, 3, 2]);
            
        case 'x'
            % Rotate around X-axis: rotate each YZ slice
            output_volume = zeros([ny, nz, nx], 'like', volume);
            for i = 1:nx
                % Extract YZ slice
                slice = squeeze(volume(i, :, :));  % size: [ny, nz]
                rotated_slice = imrotate(slice, angle_deg, crop_flag);
                % Check if size changed
                [new_ny, new_nz] = size(rotated_slice);
                if new_ny ~= ny || new_nz ~= nz
                    % Resize to original size
                    rotated_slice = imresize(rotated_slice, [ny, nz]);
                end
                output_volume(:, :, i) = rotated_slice;
            end
            % Permute back to [nx, ny, nz]
            output_volume = permute(output_volume, [3, 1, 2]);
        
        otherwise
            error('Invalid axis: %s', axis);
    end
end


% ========== EXAMPLE USAGE ==========
% 
% Example 1: Rotate around Z-axis
% volume = rand(64, 64, 64);  % Create random 3D volume
% rotated_z = rotate3D(volume, 45, 'z');
%
% Example 2: Rotate around X-axis
% rotated_x = rotate3D(volume, 30, 'x');
%
% Example 3: Rotate around Y-axis with cropping
% rotated_y = rotate3D(volume, 90, 'y', 'crop');
%
% Example 4: Sequential rotations for different axes
% rotated = rotate3D(volume, 90, 'z');
% rotated = rotate3D(rotated, 45, 'y');
% rotated = rotate3D(rotated, -15, 'x');
