% This function loads in an SDF model generated by the SDFGen tool
function [sdf] = read_sdf(filename)

    sdf.size   = dlmread(filename,' ',[0 0 0 2])'; % Dimensions of the SDF in number of cells
    sdf.origin = dlmread(filename,' ',[1 0 1 2])'; % Origin of the SDF
    sdf.dx     = dlmread(filename,' ',[2 0 2 0]); % Grid size of the SDF
    sdf.d      = zeros(sdf.size(1),sdf.size(2),sdf.size(3)); % Gridded distance values
    current = 1;
    
    data_temp = dlmread(filename,' ',[3 0 3+sdf.size(1)*sdf.size(2)*sdf.size(3)-1 0]);
    for iz=1:sdf.size(3)
%         disp(['   Progress ',num2str(100*iz/sdf.size(3)),'%']);
        for iy=1:sdf.size(2)
            for ix=1:sdf.size(1)
                sdf.d(ix,iy,iz) = data_temp(current);
                current = current + 1;
            end
        end
    end

end