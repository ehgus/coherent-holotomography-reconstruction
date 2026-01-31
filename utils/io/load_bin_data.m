function stack = load_bin_data(base_path)
% LOAD_BIN_DATA Load binary data from result.bin and specs.bin files
%
% Syntax:
%   stack = load_bin_data(base_path)
%
% Inputs:
%   base_path - String or char array specifying the directory path
%               containing 'result.bin' and 'specs.bin' files
%
% Outputs:
%   stack - Multi-dimensional array containing the loaded image stack
%           with dimensions specified in specs.bin
%
% Example:
%   data_path = 'F:\Data\TGV imaging\bg_IMAGING_000';
%   images = load_bin_data(data_path);
%
% See also: fopen, fread, reshape

    % Open and read result.bin file (contains image data)
    fileID = fopen(fullfile(base_path, 'result.bin'), 'r');
    if fileID == -1
        error('Could not open result.bin file at: %s', base_path);
    end
    buffer = fread(fileID, '*uint8');
    fclose(fileID);

    % Open and read specs.bin file (contains dimension information)
    fileID2 = fopen(fullfile(base_path, 'specs.bin'), 'r');
    if fileID2 == -1
        error('Could not open specs.bin file at: %s', base_path);
    end
    buffer_size = fread(fileID2, 'uint32');
    fclose(fileID2);

    % Reshape buffer according to dimensions specified in specs.bin
    stack = reshape(buffer(:), buffer_size(:)');

end
