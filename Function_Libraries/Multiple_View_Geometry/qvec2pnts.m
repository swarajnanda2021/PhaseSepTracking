function [ Y,g_val ] = qvec2pnts( q , X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% generate points on cirkel if not as input
if nargin<2
    X=genpntsph(100);
elseif numel(X)==1
    X=genpntsph(X);
end

% vector to levelset
[ H,g,s ] = qvec2lset(q);

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

% scale
s(:,:,g_val)=cell2mat(cellfun(@(s,g,x)sqrt(-2*s-g'*x),num2cell(s(:,:,g_val),[1 2]),num2cell(g(:,:,g_val),[1 2]),num2cell(x(:,:,g_val),[1 2]),'UniformOutput',false));

% coordinates transformation
A=nan*zeros(size(H));
A(:,:,g_val)=cell2mat(cellfun(@(H)chol(H),num2cell(H(:,:,g_val),[1 2]),'UniformOutput',false)); % pos def only

% position ellipsoid
Y=nan*zeros(size(X,1),size(X,2),size(A,3));
Y(:,:,g_val)=cell2mat(cellfun(@(s,A,x)s*(A\X)+x,num2cell(s(:,:,g_val),[1 2]),num2cell(A(:,:,g_val),[1 2]),num2cell(x(:,:,g_val),[1 2]),'UniformOutput',false));

end


% % sign positive
% q=q./q(end,:);
