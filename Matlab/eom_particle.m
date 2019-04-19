function [dX] = eom_particle(t,X,info)

    dX = zeros(6,1);
    dX(1:3) = X(4:6);
    r = X(1:3);
    g = -info.mu/(norm(r)^3)*r;
    dX(4:6) = g - cross(info.Omg,cross(info.Omg,X(1:3))) - 2*cross(info.Omg,X(4:6));

end