function !this function is not in use and to be removed! 

[ X, I, q ] = contriang( c, P )
%contriang Conic triangulation performs a point triangulation based on
%   ellipses imaged from multiple views using the forward projected light
%   cone. This triangulation method maximizes the forward scaled light
%   intensity, which scales quadratic in the depth of field (I~dof^2).
%
%   Input
%       c conic vector for multiple views 6xCxN
%       P camera matrices for each view {3x4}xC
%
%   Output
%       X triangulated coordinate 3xN
%       q conic matrix resulting from forward cones 10xN

% correct input
if iscell(c)
    c=permute(cat(3,c{:}),[1 3 2]);
end
if iscell(P)
    P=repmat(cat(3,P{:}),1,1,1,size(c,3));
end

% data size ...

% x=cvec2eshp(c);
% x=reshape(x,2,size(c,2),[]);
% I=c(1,:,:).*x(1,:,:).^2 ...
%     +c(2,:,:).*x(1,:,:).*x(2,:,:) ...
%     +c(3,:,:).*x(2,:,:).^2 ...
%     +c(4,:,:).*x(1,:,:) ...
%     +c(5,:,:).*x(2,:,:) ...
%     +c(6,:,:) ;
% figure; plot(I(:),'.')

% normalize conic
c=contrans(c,[],'normalize');

% compute forward cone quadrics
C=cvec2cmat(c);
C=reshape(C,size(C,1),size(C,2),1,size(c,2),size(c,3));
P=reshape(P,size(P,1),size(P,2),1,size(P,3),size(P,4));
Qco=squeeze(sum(repmat(permute(P,[2 1 3 4 5]),1,1,4).* ...
    repmat(permute(sum(repmat(C,1,1,4).* ...
    repmat(permute(P,[3 1 2 4 5]),3,1),2),[2 1 3 4 5]),4,1),2));
qco=qmat2qvec(Qco);
qco=reshape(qco,size(qco,1),size(c,2),size(c,3));



X=inhc2homc(x,1);
I=qco(1,:,:).*X(1,:,:).^2 + ...
    qco(2,:,:).*X(1,:,:).*X(2,:,:) + ...
    qco(3,:,:).*X(1,:,:).*X(3,:,:) + ...
    qco(4,:,:).*X(2,:,:).^2 + ...
    qco(5,:,:).*X(2,:,:).*X(3,:,:) + ...
    qco(6,:,:).*X(3,:,:).^2 + ...
    qco(7,:,:).*X(1,:,:) + ...
    qco(8,:,:).*X(2,:,:) + ...
    qco(9,:,:).*X(3,:,:) + ...
    qco(10,:,:) ;

figure; plot(I(:),'.')

qco=qco./I;

% compute triangulation quadric
q=squeeze(sum(qco,2)); % summation

% compute midpoint triangulated from forward quadric
X=qvec2eshp(q);

% compute intensity
I=q(1,:).*X(1,:).^2 + ...
    q(2,:).*X(1,:).*X(2,:) + ...
    q(3,:).*X(1,:).*X(3,:) + ...
    q(4,:).*X(2,:).^2 + ...
    q(5,:).*X(2,:).*X(3,:) + ...
    q(6,:).*X(3,:).^2 + ...
    q(7,:).*X(1,:) + ...
    q(8,:).*X(2,:) + ...
    q(9,:).*X(3,:) + ...
    q(10,:) ;

% compute depth of field position
Xc=f(P,X);

% keep valid
v=I>0 & min(Xc(3,:,:),[],2)>0;

end

