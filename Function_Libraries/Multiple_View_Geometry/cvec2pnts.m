function [ Y ] = cvec2pnts( c , X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% generate points on cirkel if not as input
if nargin<2
    X=genpntcir(100);
elseif numel(X)==1
    X=genpntcir(X);
end

% vector to levelset
[ H,g,s ] = cvec2lset(c);

% validity
g_val=ones(1,size(s,3))==1;

% make sure we plot the equiv ellipse in case flipped axis or saddles
for i=1:size(H,3)
    if ~isnan(s(:,:,i))
        
        % sign
        [~,p]=chol(H(:,:,i));
        sgn=double(p==0) - double(p>0);
        H(:,:,i)=sgn*H(:,:,i);
        g(:,:,i)=sgn*g(:,:,i);
        s(:,:,i)=sgn*s(:,:,i);
        
        % validity
        [~,p]=chol(H(:,:,i));
        if p~=0
            g_val(i)=0;
        end
        
%         % pos def
%         [R,D,~]=svd(H(:,:,i));
%         H(:,:,i)=R*D*(R');
        
    else
        
        % validity
        g_val(i)=0;
        
%         % available
%         H(:,:,i)=eye(size(H(:,:,i)));

    end
end

% triangulation vector
x=nan*zeros(size(g));
x(:,:,g_val)=cell2mat(cellfun(@(H,g)-H\g,num2cell(H(:,:,g_val),[1 2]),num2cell(g(:,:,g_val),[1 2]),'UniformOutput',false));

% scale (abs to ensure non imagionary numerics)
s(:,:,g_val)=cell2mat(cellfun(@(s,g,x)sqrt(abs(-2*s-g'*x)),num2cell(s(:,:,g_val),[1 2]),num2cell(g(:,:,g_val),[1 2]),num2cell(x(:,:,g_val),[1 2]),'UniformOutput',false));

% coordinates transformation
A=nan*zeros(size(H));
A(:,:,g_val)=cell2mat(cellfun(@(H)chol(H),num2cell(H(:,:,g_val),[1 2]),'UniformOutput',false)); % pos def only

% position ellipsoid
Y=nan*zeros(size(X,1),size(X,2),size(A,3));
Y(:,:,g_val)=cell2mat(cellfun(@(s,A,x)s*(A\X)+x,num2cell(s(:,:,g_val),[1 2]),num2cell(A(:,:,g_val),[1 2]),num2cell(x(:,:,g_val),[1 2]),'UniformOutput',false));

end


% % sign positive
% c=c./c(end,:);

% % make axis positive
% % H=cell2mat(cellfun(@(H)sign(trace(H))*H,num2cell(H,[1 2]),'Uniformoutput',false));
% % s=cell2mat(cellfun(@(H,s)sign(trace(H))*s,num2cell(H,[1 2]),num2cell(s,[1 2]),'Uniformoutput',false));
% 
% % positive definite: better do this before data into function
% P=zeros(1,size(H,3));
% for i=1:size(H,3)
%     [~,P(i)]=chol(H(:,:,i));
% end
% H=H(:,:,P==0);
% g=g(:,:,P==0);
% s=s(:,:,P==0);
% x=x(:,:,P==0);
