function [ coef ] = SGfilt_coef( dat, fobj, ord, crt )%, woi
%filt_coef_SG Up to 3-dimensional multivariate n-order convolution based
% Savitzky-Golay least squares polynomial fitting filter on uniform grid
%
%   $Author: Koen Muller$
%
% Input,
%   dat: n-dimensional input data x->dim1, y->dim2, z->dim3
%   fobj: filter object, either strel of window, note: uneven dimensions
%   ord: filter order in each dimension
%   crt: order of cross terms
%       none: 0
%       ord dim: 1 | 2 | 3
%       all: 4
%       array: custom
%   
% Output,
%   coef: local window SG-coefficients as gridded cell
%
% Note: @ midpoint eval the 1D coefficients coincide with the 0 to nth ord
%   derivative up to multiplication. When ord=win-1 & crt=4, fit is exact.
%
% 

% 0) correct input
% if nagin<5
%     woi=reshape([1+size(dat)*0 
%         size(dat)],1,[]); % [begin_dim1 end_dim1 ... etc]
% end

% 1) properties input data
switch class(fobj)
    case 'double'
        win=fobj;
        nse=ones(win)==1;
    case 'strel'
        nse=fobj.Neighborhood;
        win=ones(1,3);
        win(1:length(size(nse)))=size(nse);
end
hwin=floor(win/2); % half flank filter window

% 2) padding domain and data at boundaries
dat=2*padarray(dat,hwin+1,'replicate')...
    -padarray(dat,hwin+1,'symmetric');

% - remove replicate boundary points
dat( hwin(1)+1   , :           , :           )=[];
dat( end-hwin(1) , :           , :           )=[];
dat( :           , hwin(2)+1   , :           )=[];
dat( :           , end-hwin(2) , :           )=[];
dat( :           , :           , hwin(3)+1   )=[];
dat( :           , :           , end-hwin(3) )=[];

% remove paddaing outside window of interest

% 3) Terms polynomial model
if numel(crt)==1
    if crt==0 || crt==4 % selected terms to multivariate polynomial fit
        T=ones(ord+1)*(crt~=0);
    else
        [oy,ox,oz]=meshgrid(0:ord(2),0:ord(1),0:ord(3));
        T=(ox+oy+oz)<=ord(crt);
    end
    T(:,1,1)=1; T(1,:,1)=1; T(1,1,:)=1;
else %custom
    T=crt;
end
T=reshape(T,[],1); % set of terms

% 4) Defined local window indexing, with midpoint centered @stencil
W=cell(size(hwin));
[W{2},W{1},W{3}]=meshgrid((-hwin(2):hwin(2)),...
    (-hwin(1):hwin(1)),...
    (-hwin(3):hwin(3)));

% 4) (limited order) van der Monde matrix
% - Find the vd Monde matrix in each dimension separate
M1D=cell(size(W));
for d=1:3
    M1D{d}=bsxfun(@power, W{d}(nse) ,(0:ord(d))); %here neighborhood shape matters
end

% - Find all combinations of cross terms for each window index
M=zeros([nnz(nse) ord+1]);
for w=1:size(M,1)
    [ty,tx,tz]=meshgrid(M1D{2}(w,:),M1D{1}(w,:),M1D{3}(w,:));
    M(w,:,:,:)=tx.*ty.*tz;
end
M=reshape(M,size(M,1),[]); % reshape to matrix
M(:,T==0)=[]; % remove excluded terms

% - Defined final discrete operator by vd Monde matrix
A=(M'*M)\M'; % here all relevant xyz term couple by inversion, conv zero can therfore not be seen as altering domain

% 5) filter coefficients
coef=cell(1,length(T)); % all possible coefficients
T=find(T); % convert to indices
h=zeros(win);
for t=1:length(T) % loop over indices
    % filter object
    h(nse)=A(t,:);
    
    % coefficients
    c=convn(dat,flip(flip(flip(h,1),2),3),'valid'); 
    
    % write LS solution without boundary padding
    coef{T(t)}=c;
    
end

% reshape coef to easily workable cell struct
coef=reshape(coef,ord+1);

end