function [R,Q]=rq(A)
% function [R Q]=rq(A)
% Author: Koen Muller
% Last Update: 27/06/2017 
% note: This is the simplest version of that.

[Q,R]=qr(flipud(A)');
R=rot90(R',2);
Q=flipud(Q');

end