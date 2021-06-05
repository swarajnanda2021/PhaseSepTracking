function [ output_args ] = regpoldec( x,y,imd )
%regpoldec Regional Polynomial (Image) Decomposition
%   This is an iterative segmentation approach with multiple segmentation
%   level set by the domain of the residual, i.e. regions defined by what
%   is left.
%
%   28-10-2017: first idea, xy second order fit to image data
%   future: different order xy-t, xyz(-t) etc. as a general data
%   decomposition

X=reshape(x,[],1);
Y=reshape(y,[],1);

obs=[ones(size(X)) X Y X.*Y X.^2 Y.^2];

dat=reshape(imd,[],1);

fit=zeros(size(dat)); % could also ma a residual/ normalized fit for conic sections

lab=zeros(size(dat));
    
conn=8;
trs=nnz(ctrl.fobj.Neighborhood);
fobj=strel('disk',50);
% while
    
    reg=zeros(size(dat));
    for j=unique(lab)'
        disp(num2str(j));
        ind=lab==j;
        
        dom=reshape(ind,size(imd));
        
        c=regress(dat(ind),obs(ind,:));
        
        fit(ind)=obs(ind,:)*c;
        
        seg=imd>reshape(fit,size(imd));
        
        seg=imopen(seg & dom,strel('square',3));
        seg=imclose(seg & dom,strel('square',3));
        
        seg=bwareaopen(  seg & dom ,trs,conn); % remove small aere
        seg=bwareaopen(  ~seg & dom ,trs,conn); % remove small aere
        
%         seg=bwmorph(seg & dom,'erode',3);
%         seg=bwmorph(seg,'dilate',3);

        seg=~seg & dom;
        
        segp=  seg & dom;
        segm= ~seg & dom;
                
        labp=reshape(bwlabel(segp,conn),[],1);
        labm=reshape(bwlabel(segm,conn),[],1);
        
        indp=reshape(segp,[],1);
        indm=reshape(segm,[],1);
        
        reg(indp)=max(reg)+labp(indp);
        reg(indm)=max(reg)+labm(indm);
        
    end
    lab=reg;
    % seg=imd>imfilter(imd,double(fobj.Neighborhood))/nnz(fobj.Neighborhood);
% end

end