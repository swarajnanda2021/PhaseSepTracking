function [dis,cor,wav,spe,nor] = spatial2pointcor(grd,dat1,dat2,typ)
%spatial2pointcor Spatial 2-Point Corrolation, computer the 2-point
%   corrolation based on a fourrier transform for gridded vector data.
%
%   Input
%       grd Indexed uniform rectangular grid for the data
%       dat1 Indexed reference data [(1+dim1)xN2]
%       dat2 Indexed shifted data [(1+dim2)xN2] 
%           (dim2=>dim1, as we only export uppertriang of symmetric matrix)
%           (when dat2=dat1, autocorrolation)
%       typ direct | normalized (default), use normalized when missing data,
%          otherwise direct gives better quality spectrum.
%       
%   Output
%       dis Distance [3x3xfloor(N1,N2)]
%       cor Symmetric corrolation matrix [dim1xdim2xfloor(N1,N2)]
%       wav Wavenumber [3x3xfloor(N1,N2)]
%       psp Symmetric powerspectrum matrix [dim1xdim2xfloor(N1,N2)]
%
%   Note: 
%       Matrix element numbering by symmetry
%           [a1 a2 a3 a4       ... aj
%             a2 a(j+1) a(j+2) ... etc
%             ...
%             ai               ... end]
%       The Normalization variable can also be used to flag empty
%           corrolations performed.
%       We export the the whole domain, including its symmetry for higher
%           dimensional corrolations.
%       Work for 1D to 3D data.
%
%   TO BE CHECKED consistent shift and inverse

% Correct data size
if nargin<4
    typ='direct'; % simplest
end
if size(grd,1)-1==2
    grd=[grd 
        zeros(1,size(grd,2))];
end

% % Find shared index
% cind=dat1(1,ismember(dat1(1,:),dat2(1,:))); % grid indexing corrolation

% Data dimensions
dim1=size(dat1,1)-1;
dim2=size(dat2,1)-1;

% error message incorrect input dimensions
if dim2<dim1
    error('spatial2pointcor.m: dim1<dim2 gives wrong export results as dim2>=dim1.')
end

% Domain weight
Wght1=zeros(1,size(grd,2)); % don't need indexing anymore
Wght2=zeros(1,size(grd,2));
Wght1(:,dat1(1,:))=1;
Wght2(:,dat2(1,:))=1;

% Assemble data
Dat1=zeros(dim1,size(grd,2)); % don't need indexing anymore
Dat2=zeros(dim2,size(grd,2));
Dat1(:,dat1(1,:))=dat1(2:end,:);
Dat2(:,dat2(1,:))=dat2(2:end,:);

% get grid span
x=unique(grd(2,:));
y=unique(grd(3,:));
z=unique(grd(4,:));

% Get grid sampling
Sx=unique(gradient(x));
Sy=unique(gradient(y));
Sz=unique(gradient(z));

% Data size
isiz= length(x);
jsiz= length(y);
ksiz= length(z);

% Reshape data in cells
DAT1=cell(1,dim1);
for k=1:dim1
    DAT1{k}=reshape(Dat1(k,:),isiz,jsiz,ksiz);
end
DAT2=cell(1,dim2);
for k=1:dim2
    DAT2{k}=reshape(Dat2(k,:),isiz,jsiz,ksiz);
end

% Reshape domain weight
WGHT1=reshape(Wght1,isiz,jsiz,ksiz);
WGHT2=reshape(Wght2,isiz,jsiz,ksiz);

%figure; surf(X(:,:,ceil(ksiz/2)),Y(:,:,ceil(ksiz/2)),DAT1{1}(:,:,ceil(ksiz/2))); view(2); shading interp
%figure; surf(X(:,:,ceil(ksiz/2)),Y(:,:,ceil(ksiz/2)),DAT2{1}(:,:,ceil(ksiz/2))); view(2); shading interp

%figure; surf(X(:,:,ceil(ksiz/2)),Y(:,:,ceil(ksiz/2)),WGHT1(:,:,ceil(ksiz/2))); view(2); shading interp
%figure; surf(X(:,:,ceil(ksiz/2)),Y(:,:,ceil(ksiz/2)),WGHT2(:,:,ceil(ksiz/2))); view(2); shading interp

%-- Compute Fourrier Transforms

% Define sampling frequencies
kx=(-(isiz-1)/2:(isiz-1)/2)*(isiz/Sx); % scaling [0-1]*maxwavnumber
ky=(-(jsiz-1)/2:(jsiz-1)/2)*(jsiz/Sy);
kz=(-(ksiz-1)/2:(ksiz-1)/2)*(ksiz/Sz);

% Define grid
[kX,kY,kZ]=ndgrid(kx,ky,kz);

% Fourrier transform data
FOU1=cell(1,dim1);
for k=1:dim1
    FOU1{k}=fftshift(fftn(DAT1{k})); % centered fft
end
FOU2=cell(1,dim2);
for k=1:dim2
    FOU2{k}=fftshift(fftn(DAT2{k})); % centered fft
end

% Fourrier transform weight
WGHT1=fftshift(fftn(WGHT1));
WGHT2=fftshift(fftn(WGHT2));

%figure; surf(kX(:,:,ceil(ksiz/2)),kY(:,:,ceil(ksiz/2)),real(FOU1{1}(:,:,ceil(ksiz/2)))); view(2); shading interp
%figure; surf(kX(:,:,ceil(ksiz/2)),kY(:,:,ceil(ksiz/2)),imag(FOU1{1}(:,:,ceil(ksiz/2)))); view(2); shading interp

%-- Compute Spectrum and Corrolation

% Compute spectrum
SPE=cell(dim1,dim2);
for i=1:dim1
    for j=i:dim2 % from i speed up by symmetry
        SPE{i,j}=FOU2{i}.*conj(FOU1{j});
    end
end
WGHT=WGHT2.*conj(WGHT1);

%figure; surf(kX(:,:,ceil(ksiz/2)),kY(:,:,ceil(ksiz/2)),real(SPE{1,1}(:,:,ceil(ksiz/2)))); view(2); shading interp
%figure; surf(kX(:,:,ceil(ksiz/2)),kY(:,:,ceil(ksiz/2)),imag(SPE{1,1}(:,:,ceil(ksiz/2)))); view(2); shading interp

% Define sampling frequencies
dx=x-median(x); % scaling [0-1]*maxwavnumber
dy=y-median(y);
dz=z-median(z);

% Define grid
[dX,dY,dZ]=ndgrid(dx,dy,dz);

% Compute corrolation function
COR=cell(dim1,dim2);
for i=1:dim1
    for j=i:dim2 % from i speed up by symmetry
        COR{i,j}=fftshift(ifftn(ifftshift(SPE{i,j}))); %fftshift(real(fftn(ifftshift(PSP{i,j}))));%
    end
end

%figure; surf(dX(:,:,ceil(ksiz/2)),dY(:,:,ceil(ksiz/2)),COR{1,1}(:,:,ceil(ksiz/2))); view(2); shading interp

%-- Compute normalization (data might be empty in places)

% Compute normalization
NOR=fftshift(ifftn(ifftshift(WGHT)));
NOR(NOR<0.5)=inf; % are not part of the domain

%figure; plot(unique(NOR(:)),'.')
%figure; surf(dX(:,:,ceil(ksiz/2)),dY(:,:,ceil(ksiz/2)),NOR(:,:,ceil(ksiz/2))); view(2); shading interp

% Compute normalized corrolation function and spectrum
if strcmp(typ,'normalized')
    for i=1:dim1
        for j=i:dim2 % from i speed up by symmetry
            COR{i,j}=COR{i,j}./NOR;
            SPE{i,j}=fftshift(fftn(ifftshift(COR{i,j})));
        end
    end
end

%figure; surf(dX(:,:,ceil(ksiz/2)),dY(:,:,ceil(ksiz/2)),COR{1,1}(:,:,ceil(ksiz/2))); view(2); shading interp

%figure; surf(kX(:,:,ceil(ksiz/2)),kY(:,:,ceil(ksiz/2)),real(SPE{1,1}(:,:,ceil(ksiz/2)))); view(2); shading interp
%figure; surf(kX(:,:,ceil(ksiz/2)),kY(:,:,ceil(ksiz/2)),imag(SPE{1,1}(:,:,ceil(ksiz/2)))); view(2); shading interp

%-- Export Variables

% number number corroaltions
ncor=dim2*dim1-dim1*(dim1-1)/2; % n*k?k/2*(k+1) with k=dim1<=n=dim2

% Define wavenumbers
wav=[kX(:)' 
    kY(:)' 
    kZ(:)'];

% Spectrum data
spe=zeros(ncor,size(grd,2));
k=0; % count elements
for i=1:dim1
    for j=i:dim2 % from i speed up by symmetry
        k=k+1;
        spe(k,:)=reshape(SPE{i,j},1,[]);
    end
end
% wav=wav(:,cind);
% psp=psp(:,cind);

% Define wavenumbers
dis=[dX(:)' 
    dY(:)' 
    dZ(:)'];

% Corrolation data
cor=zeros(ncor,size(grd,2));
k=0; % count elements
for i=1:dim1
    for j=i:dim2 % from i speed up by symmetry
        k=k+1;
        cor(k,:)=reshape(COR{i,j},1,[]);
    end
end
% dis=dis(:,cind);
% cor=cor(:,:,cind);

% Normalization
nor=reshape(NOR,1,[]);

end


% % convolute in fourrier space
% Q = real(fftn(PSP));
% Q = Q/(isize*jsize*ksize-1)/isize/jsize/ksize; % 1/N*dxdydz some scalar normalization
% 
% % shift
% cor = fftshift(Q); % center the fourrier transform
