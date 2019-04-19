% This function samples a polyhedron shape model at field point r
% It computes the minimum distance d between r and the polyhedron
% and provides the normal vector of the corresponding closest facet
function [d,N] = sample_poly(r,poly)

    r1 = r - poly.pts(:,poly.tri(:,1));
    r2 = r - poly.pts(:,poly.tri(:,2));
    r3 = r - poly.pts(:,poly.tri(:,3));
    r21 = r2 - r1;
    r32 = r3 - r2;
    r13 = r1 - r3;
    n = cross(r21,r13);
    u1 = min(0,max(dot(r21,r1)/dot(r21,r21),1))*r21 - r1;
    u2 = min(0,max(dot(r32,r2)/dot(r32,r32),1))*r32 - r2;
    u3 = min(0,max(dot(r13,r3)/dot(r13,r13),1))*r13 - r3;
    q = abs(sign(dot(r1,cross(r21,n))) + sign(dot(r2,cross(r32,n))) + sign(dot(r3,cross(r13,n))));
    d1 = min([dot(u1,u1);dot(u2,u2);dot(u3,u3)]);
    d2 = (dot(n,r1).^2)./dot(n,n);
    df = d2;
    df(q<2) = d1(q<2);
    [d2min,imin] = min(df);
    d = sign(dot(r1(:,imin),poly.fN(:,imin)))*sqrt(d2min);
    N = poly.fN(:,imin);

end