function [value,isterminal,direction] = event_particle(t,X,info)

    isterminal = 1;
    direction = 0;
    value = sample_sdf(X(1:3),info.sdf);

end