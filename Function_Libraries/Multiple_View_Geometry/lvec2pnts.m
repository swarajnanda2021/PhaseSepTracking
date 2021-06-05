function [ x ] = lvec2pnts( l , boun )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% correct input
if size(l,1)==2
    l=inhc2homc(l);
end
if sum(boun(:,1)==boun(:,end))~=2
   boun(:,end+1)=boun(:,1); % close boundary 
end
if size(boun,1)==2
    boun=inhc2homc(boun);
end

% loop boundary points
nx=[1 0]';
ny=[0 1]';
xmin=zeros([2 size(l,2)]);
xmax=zeros([2 size(l,2)]);
for i=1:size(boun,2)-1
    % compute intersection with polygon vertex
    vboun=[boun(:,i),boun(:,i+1)]; % boundary vertex
    lboun=cross(vboun(:,1),vboun(:,2)); % boundary line
    xboun=homc2inhc(cross(repmat(lboun,1,size(l,2)),l)); % boundary coordinates
    
    % rotate (because of defined vertical / horizontal lines)
    nl=[lboun(1) lboun(2)]'/sqrt(lboun(1)^2+lboun(2)^2);
    R=rotmz(-sign(dot(nl,nx))*acosd(dot(nl,ny))); R=R(1:2,1:2);
    xdat=nx'*(R\xboun); % w.r.t. base line frame
    vlim=sort(nx'*(R\homc2inhc(vboun)));
    valx=xdat(1,:)>=vlim(1) & xdat(1,:)<=vlim(2);
    
    % write
    sel=sum(xmin==0,1)==2;
    xmin(:,sel & valx)=xboun(:,sel & valx);
    xmax(:,~sel & valx)=xboun(:,~sel & valx);
end
x=[xmin;xmax];


end