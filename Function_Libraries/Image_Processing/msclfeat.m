function [ X, U, C, I ] = msclfeat( imd, fobj , typ , foi, nscl )
%msclfeat overlay for features from image over multiple scales

% fine and back scale define by filter object
sig=(size(fobj.Neighborhood,1)-1)./(2*3);
sobj=double(getnhood(strel('disk',1*sig,0))); sobj=sobj/nnz(sobj);
bobj=double(getnhood(strel('disk',7*sig,0))); bobj=bobj/nnz(bobj);
rbox=2*[1 1];

% define image grid
[y,x]=meshgrid(1:size(imd,2),1:size(imd,1));

% ellipse function
efun=@(c,x,y)c(1)*x.^2+c(2)*x.*y+c(3)*y.^2+c(4)*x+c(5)*y+c(6);

% loop scales
count=1;
while count<nscl %&& numel(obj(:,:,1))>prod(!fobj)
    % define scalefspecial('average', sbox)
    scl=imfilter(imd,sobj,'symmetric'); % filter noise to scale
    
    % define backgroundfspecial('average', bbox)[~,~,~,bgr]=imbox(scl,bbox,'avg'); % background
    bgr=imfilter(imd,bobj,'symmetric'); % filter noise to scale
    
    % define objects to interest
    obj=scl-bgr; % object[ive] scale
    
    % detect new features ,int
    [pos,vel,con,int]=featpoints(obj,fobj,typ,foi); % ,lset,L
    
    ime=imevallset( x , y , efun , con );
    
    % unit ellipse distance between point
    E=con*([pos(:,1).^2 pos(:,1).*pos(:,2) pos(:,2).^2 pos(:,1) pos(:,2) ones(size(pos,1),1)]');
    E=sqrt(abs(1-E./repmat(diag(E),1,size(E,1))));
    
    % back forth within ellipse adjacency
    A=E<1 & E'<1;
    
    
    % pixel to scale transformation matrix
    T=[diff(x(1:2,1)) 0 x(1,1)-diff(x(1:2,1))
        0 diff(y(1,1:2)) y(1,1)-diff(y(1,1:2))
        0 0 1]; % -1 pix shift (T=Tcal*Tpix, did the math)
    
    % box the average image to desired scale
    imd=imbox(bgr,rbox,'avg'); % average to lower resolution
    [x,y] = boxgrid( x,y,[], [] ,rbox);
    
    % count next number scale
    count=count+1;
    
end

end



%     % padarray to expand regional labelling for interpolating to reference
%     L=padarray(L,[1 1],'replicate');
%     X=2*padarray(X,[2 2],'replicate')-padarray(X,[2 2],'symmetric');
%     Y=2*padarray(Y,[2 2],'replicate')-padarray(Y,[2 2],'symmetric');
%     X(2,:)=[]; X(:,2)=[]; X(end-1,:)=[]; X(:,end-1)=[];
%     Y(2,:)=[]; Y(:,2)=[]; Y(end-1,:)=[]; Y(:,end-1)=[];
%     


% grid to manipulate
% Lmap=zeros(size(imd,1),size(imd,2));
% X=x;
% Y=y;


% 
%     % interpolate label matrix to reference grid
%     L=interp2(Y,X,L,y,x,'nearest',0); % outside X and Y is not in estimation
%     
%     % stack label map
%     L(L>0)=L(L>0)+max(Lmap(:));
%     Lmap=cat(3,L,Lmap);
%     

% 
%     % transform levelset
%     lset(:,1:6)=contrans(lset(:,1:6)',T)';
%     lset(!)=transform here;
%     


%     ! revote concutive scales??
% 
% ! revote here?