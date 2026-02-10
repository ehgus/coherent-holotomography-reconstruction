function unlock_tomocube_dat(src_file, dst_dir)
% UNLOCK_TOMOCUBE_DAT Convert Tomocube HDF5 images.dat to PNG stack.
%
% Syntax:
%   unlock_tomocube_dat(src_file, dst_dir)
%
% Inputs:
%   src_file - Path to images.dat or calibration.dat (HDF5)
%   dst_dir  - Output directory to store PNG files
%
    if ~exist(dst_dir, 'dir')
        mkdir(dst_dir);
    end

    img_cnt = length(h5info(src_file, '/images').Datasets);
    for i = 0:img_cnt-1
        h5loc = sprintf('/images/%06d', i);
        img_binary = h5read(src_file, h5loc);
        dst_file = fullfile(dst_dir, sprintf('%06d.png', i));
        file_id = fopen(dst_file, 'w');
        fwrite(file_id, img_binary);
        fclose(file_id);
    end
end
