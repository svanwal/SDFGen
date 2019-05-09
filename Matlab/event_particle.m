function [value,isterminal,direction] = event_particle(t,X,info)

    isterminal = 1;
    direction = 0;
    sdf_min = info.sdf.origin;
    sdf_max = info.sdf.origin + info.sdf.dx*info.sdf.size;
    if sum(X(1:3)<sdf_min)>0 || sum(X(1:3)>sdf_max)>0
        value = 999;
    else
        value = sample_sdf(X(1:3),info.sdf);
    end

end