function [ nodes ] = genpntsph( N )
%generate points on sphere for meshing triangles, generate by fibonacci seq
% http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere

% evenly space point on unit circle
dlong = pi*(3-sqrt(5)) ; % step along golden ratio
dz    = 2.0/N ; % 
long = 0 ;
z    = 1 - dz/2 ;
nodes=zeros(3,N);
for k = 1 : N
    r    = sqrt(1-z*z);
    nodes(:,k) = [ cos(long)*r, sin(long)*r, z ];
    z    = z - dz;
    long = long + dlong;
end

end