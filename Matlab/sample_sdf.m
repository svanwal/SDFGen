function [d,N,v] = sample_sdf(r,sdf)

    sdf_min = sdf.origin;
    sdf_max = sdf.origin + sdf.dx*sdf.size;
    
    if sum(r<sdf_min)>0 || sum(r>sdf_max)>0
        
        d = sdf.dx;
        N = [1;0;0];
        
    else

        rho = r - sdf.origin;
        v = min(sdf.size-1,max(1,floor(rho/sdf.dx) + 1)); % Clamped between 1 and sdf.size
        r_bar = (rho - (v-1)*sdf.dx)/sdf.dx;

        % This works now
        d000 = sdf.d(v(1),   v(2),   v(3));
        d100 = sdf.d(v(1)+1, v(2),   v(3));
        d010 = sdf.d(v(1),   v(2)+1, v(3));
        d001 = sdf.d(v(1),   v(2),   v(3)+1);
        d110 = sdf.d(v(1)+1, v(2)+1, v(3));
        d101 = sdf.d(v(1)+1, v(2),   v(3)+1);
        d011 = sdf.d(v(1),   v(2)+1, v(3)+1);
        d111 = sdf.d(v(1)+1, v(2)+1, v(3)+1);
        a0 = d000;
        a1 = d100 - d000;
        a2 = d010 - d000;
        a3 = d001 - d000;
        a4 = d000 - d100 - d010 + d110;
        a5 = d000 - d010 - d001 + d011;
        a6 = d000 - d100 - d001 + d101;
        a7 = -d000 + d100 + d010 + d001 - d110 - d011 - d101 + d111;
        d = a0 + a1*r_bar(1) + a2*r_bar(2) + a3*r_bar(3) ...
               + a4*r_bar(1)*r_bar(2) + a5*r_bar(2)*r_bar(3) + a6*r_bar(1)*r_bar(3) ...
               + a7*r_bar(1)*r_bar(2)*r_bar(3);
        grad = 1/sdf.dx*[a1 + a4*r_bar(2) + a6*r_bar(3) + a7*r_bar(2)*r_bar(3);
                         a2 + a5*r_bar(3) + a4*r_bar(1) + a7*r_bar(1)*r_bar(3);
                         a3 + a6*r_bar(1) + a5*r_bar(2) + a7*r_bar(1)*r_bar(2)];
        N = grad/norm(grad);
        
    end
    
end