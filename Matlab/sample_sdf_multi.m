function [d,N,bool_outside] = sample_sdf_multi(r,sdf)
    
    rho = r - sdf.origin;
    v0 = floor(rho/sdf.dx) + 1;
    outside3 = (v0<1) + (v0>sdf.size-1);
    outside1 = sum(outside3,1);
    v = min(sdf.size-1,max(1,v0));
    r_bar = (rho - (v-1)*sdf.dx)/sdf.dx;
    % v
    v000 = sub2ind(sdf.size,v(1,:),   v(2,:),   v(3,:));
    v100 = sub2ind(sdf.size,v(1,:)+1, v(2,:),   v(3,:));
    v010 = sub2ind(sdf.size,v(1,:),   v(2,:)+1, v(3,:));
    v001 = sub2ind(sdf.size,v(1,:),   v(2,:),   v(3,:)+1);
    v110 = sub2ind(sdf.size,v(1,:)+1, v(2,:)+1, v(3,:));
    v101 = sub2ind(sdf.size,v(1,:)+1, v(2,:),   v(3,:)+1);
    v011 = sub2ind(sdf.size,v(1,:),   v(2,:)+1, v(3,:)+1);
    v111 = sub2ind(sdf.size,v(1,:)+1, v(2,:)+1, v(3,:)+1);
    % d
    d000 = sdf.d(v000);
    d100 = sdf.d(v100);
    d010 = sdf.d(v010);
    d001 = sdf.d(v001);
    d110 = sdf.d(v110);
    d101 = sdf.d(v101);
    d011 = sdf.d(v011);
    d111 = sdf.d(v111);
    % a
    a0 = d000;
    a1 = d100 - d000;
    a2 = d010 - d000;
    a3 = d001 - d000;
    a4 = d000 - d100 - d010 + d110;
    a5 = d000 - d010 - d001 + d011;
    a6 = d000 - d100 - d001 + d101;
    a7 = -d000 + d100 + d010 + d001 - d110 - d011 - d101 + d111;
    % interp
    d = a0 + a1.*r_bar(1,:) + a2.*r_bar(2,:) + a3.*r_bar(3,:) ...
           + a4.*r_bar(1,:).*r_bar(2,:) + a5.*r_bar(2,:).*r_bar(3,:) + a6.*r_bar(1,:).*r_bar(3,:) ...
           + a7.*r_bar(1,:).*r_bar(2,:).*r_bar(3,:);
    grad = 1/sdf.dx*[a1 + a4.*r_bar(2,:) + a6.*r_bar(3,:) + a7.*r_bar(2,:).*r_bar(3,:);
                     a2 + a5.*r_bar(3,:) + a4.*r_bar(1,:) + a7.*r_bar(1,:).*r_bar(3,:);
                     a3 + a6.*r_bar(1,:) + a5.*r_bar(2,:) + a7.*r_bar(1,:).*r_bar(2,:)];
    % replacing outside points
    d(outside1>0) = sdf.d(v000(outside1>0)) - sdf.dx;
    grad(:,outside1>0) = sdf.origin + sdf.dx*v000(:,outside1>0);
    % finishing the normal
    Nn = sqrt(grad(1,:).^2 + grad(2,:).^2 + grad(3,:).^2);
    N = [grad(1,:)./Nn;grad(2,:)./Nn;grad(3,:)./Nn];
    bool_outside = outside1>0;
    
end