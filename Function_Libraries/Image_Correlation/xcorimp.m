function R = xcorimp(s1,s2)
%pwXCorrelate cross-correlate two images
%   R = pwXCorrelate(S1,S2) computes the cross-correlation of the images in S1 and S2.

% written by Jerry Westerweel 2007
% $Id: pwXCorrelate.m,v 1.1 2015/07/15 16:06:01 jwe Exp $

% koen: multidimensional extension to include 3d corrolation

% optional: add inclusion in power-2 image size

[isize,jsize,ksize] = size(s1);

f1 = fftn(s1);
f2 = fftn(s2);

Q = real(fftn(f1.*conj(f2)));

Q = Q/(isize*jsize*ksize-1)/isize/jsize/ksize;

% R = pwSwapPanels(Q);
R = fftshift(Q);

end