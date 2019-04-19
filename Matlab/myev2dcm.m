function M  = myev2dcm (phi, e)
%MYEV2DCM transfers a principal rotation vector into a rotation matrix

    c = cos(phi);
    s = sin(phi);
    Z = 1 - c;

    M = [e(1)^2*Z + c,          e(1)*e(2)*Z + e(3)*s,   e(1)*e(3)*Z - e(2)*s;
        e(2)*e(1)*Z - e(3)*s,   e(2)^2*Z + c,           e(2)*e(3)*Z + e(1)*s;
        e(3)*e(1)*Z + e(2)*s,   e(3)*e(2)*Z - e(1)*s,   e(3)^2*Z + c];

end