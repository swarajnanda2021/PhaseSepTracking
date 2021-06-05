function [ Y ] = dispeval( coef, fdat , der )
%dispeval Polynomial displacement field from a set of trajectories
%
%   Input
%       coef coefficients [oX oY oZ oT cX cY cZ]
%       fdat trajectory data [Ni Xi Yi Zi]
%       der derivative of interest
%
%   Output
%       Y any ouput desired, e.g. displacement seeding points  
%
%   Note: 
%       Supporting a grid indexing could impose simultaneous solution of
%       multiple vectors and impose consistentency contraints on the nodes.
%

% correct input
if nargin < 2
    fdat=zeros(4,1); % midpoint [n X]
end
if nargin < 3
    der =[0 0 0 1]; % give only time derivative, output velocity
end
if length(der)==3 % 2D
    der=[der(1:2) 0 der(3)]; % add z=0 plane
end
if size(fdat,1)==3 % 2D
    fdat=[fdat
        zeros(1,size(fdat,2))]; % add z=0 plane
end

% data size
sizX=size(fdat,1)-1; % always 3D

% assemble n
Ni=fdat(1,:); % time from initial data

% assemble X
Xi=fdat(2:4,:); % initial data / always 3D

% order indexing coeffients
Ord=coef(1:4,:); % always 3D

% coefficient data
Cdat=coef(5:7,:); % always 3D

% get derivative order and coefficients
OrdD=Ord-repmat(der',1,size(Ord,2));

% valid coefficients in derivatives
g_keep=min(OrdD,[],1)>=0;
Ord=Ord(:,g_keep);
Cdat=Cdat(:,g_keep);
OrdD=OrdD(:,g_keep);

% compute weight by differentiation of higher order terms
W=factorial(Ord)./factorial(OrdD);

% coeffient vector
c=reshape(Cdat',[],1); % always 3D

% Compute basis function expansion
phi=prod( permute(W,[3 2 1]).*permute([Xi' Ni'],[1 3 2]).^permute(OrdD,[3 2 1]) , 3 );

% Assemble blockdiagonal for each physical dimension
Phi=kron(speye(sizX),phi);

% evaluate output
y=Phi*c;

% format output
Y=reshape(y,[],sizX)';

end

