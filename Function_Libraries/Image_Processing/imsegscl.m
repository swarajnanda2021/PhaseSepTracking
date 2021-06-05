function [ mskp , mskm , imd ] = imsegscl( imd , ord , sig, iter)
%imsegscl advanced segmentation script to segment image intensity values at
%   scale define by fobj, plus and min || we focus on spatial segmentation
%   can be generalized to temporal as well by fitting a temporal surface

% correct input
if nargin<2 || isempty(ord)
    ord=[0 0 0];
end
if numel(ord)==1
    ord=ord*[1 1 0];
end
if numel(ord)==2
    ord=[ord 0];
end
if nargin<3 || isempty(sig)
    sig=3;
end
if nargin<4 || iter==0 || isempty(iter)
    iter=1;
end

% define an image grid
[y,x,t]=meshgrid(1:size(imd,2),1:size(imd,1),1:size(imd,3));

% order poly-coefficient space time
[oy,ox,ot]=meshgrid(0:ord(2),0:ord(1),0:ord(3));
T=(ox+oy+ot)<=max(ord);
ox=ox(T)';
oy=oy(T)';
ot=ot(T)';

% x reference
Xref=(x(:).^ox).*(y(:).^oy).*(t(:).^ot); % base

% segmentation plus, include 0 for generality
segp=imd>=0;

% loop
for i=1:iter
    % fit a surface to positive part image
    Y=imd(segp); % intensity
    X=(x(segp).^ox).*(y(segp).^oy).*(t(segp).^ot); % base
    
    % estimate coefficients
    [coef,~,res]=regress(Y,X);
    
    % deviation by std for segmentation 
    devp=sqrt(mean(res.^2));%std(res);
    
    % evalute fit on reference grid
    fit=reshape(Xref*coef,size(imd));
    
    % explicit removal data from image
    imd(imd>(fit+sig*devp))=fit(imd>(fit+sig*devp))+sig*devp;
    
end

% define mask
mskp=imd>=(fit+sig*devp);

if nargout>=2
    % segmentation min, include 0 for generality
    segm=imd<=0;
    
    % loop
    for i=1:iter
        % fit a surface to positive part image
        Y=imd(segm); % intensity
        X=(x(segm).^ox).*(y(segm).^oy).*(t(segm).^ot); % base
        
        % estimate coefficients
        [coef,~,res]=regress(Y,X);
        devm=sqrt(mean(res.^2));%std(res);
        
        % evalute fit on reference grid
        fit=reshape(Xref*coef,size(imd));
        
        % explicit removal data from image
        imd(imd<=(fit-sig*devm))=fit(imd<=(fit-sig*devm))-sig*devm;
        
    end
    
    % define mask
    mskm=imd<=(fit-sig*devm);
end

end

