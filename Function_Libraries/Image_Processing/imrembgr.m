function [imd,scl] = imrembgr( imd , imb , ord , sig , iter )
%imfitbgr Image fit background using polynomial regression
%   spatial scaling intensity values to best fit background
%   this allows correcting flickering background
%   sign sig back/for grnd

% correct input
if nargin < 3 || isempty(ord)
    ord=[0 0 0];
end
if numel(ord)==1
    ord=ord*[1 1 0];
end
if numel(ord)==2
    ord=[ord 0];
end
if nargin<4 || isempty(sig)
    sig=3;
end
if nargin<5 || iter==0 || isempty(iter)
    iter=1;
end

% define an image grid
[y,x,t]=meshgrid(1:size(imd,2),1:size(imd,1),1:size(imd,3));

% order poly-coefficient space time
[oy,ox,ot]=meshgrid(0:ord(2),0:ord(1),0:ord(3));
T=(ox+oy+ot)<=max(ord) ;
ox=ox(T)';
oy=oy(T)';
ot=ot(T)';

% x reference
Xref=(x(:).^ox).*(y(:).^oy).*(t(:).^ot); % base

% divided data for bgr fit
div=imd./imb;

% mask include 0 for data (bound) exclude 0 background (div)
bnd=imd>0 & imb>0;

% loop 
for i=1:iter
    % scaling substraction: fit / regress polynomial
    Y=div(bnd); % intensity
    X=(x(bnd).^ox).*(y(bnd).^oy).*(t(bnd).^ot); % base
    
    % estimate coefficients
    [coef,~,res]=regress(Y,X);
    
    % deviations
    devp=std(res(res>0));
    devm=std(res(res<0));
    
    % evaluate scaling
    scl=reshape(Xref*coef,size(imd)); % image normalization
        
    % explicit removal data from image
    div(div>(scl+sig*devp))=scl(div>(scl+sig*devp))+sig*devp;
    div(div<(scl-sig*devm))=scl(div<(scl-sig*devm))-sig*devm;
    
end

% residual
imd=imd-scl.*imb; % segmentation

end
