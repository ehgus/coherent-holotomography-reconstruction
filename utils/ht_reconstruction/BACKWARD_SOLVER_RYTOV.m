classdef BACKWARD_SOLVER_RYTOV < BACKWARD_SOLVER_BORN
    methods
        function field_fft = extractField(h, output_field)
            % Extract Rytov field by phase unwrapping
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
            UsRytov = log(amp) + 1i*phase;
            field_fft = gather(fft2(gpuArray(UsRytov)));
        end
    end
end