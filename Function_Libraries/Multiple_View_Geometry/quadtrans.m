function [ q ] = quadtrans( q , H , typ )
%quadtrans Quadric transformation
%   input 1: quadric vector
%   input 2: transformation input
%   input 3: type of transform see code below
%
%	analog to contrans.m

% correct input
if nargin<3
    typ='transmat';
end

% assemble transformation matrix
switch typ 
    case 'transmat' % transformation matrix - 3x3xN
				
		% make quadric matrix
		Q=qvec2qmat(q);

		% make same size if needed
		if size(H,3)<size(Q,3)
			H=repmat(H,1,1,size(Q,3));
		end
		if size(Q,3)<size(H,3)
			Q=repmat(Q,1,1,size(H,3));
		end

		% transform conic
		if ~isempty(Q) && ~isempty(H)
			Q=cell2mat(cellfun(@(Q,H)inv(H)'*Q*inv(H),...
				num2cell(Q,[1 2]),num2cell(H,[1 2]),'UniformOutput',false));
		end

		% write as vector
		q=qmat2qvec(Q);

    case 'displace' % displace new position
	
%         T=reshape(T,size(T,1),1,size(T,2));
%         T=cell2mat(cellfun(@(x)[1 0 x(1)
%             0 1 x(2)
%             0 0 1],num2cell(T,[1 2]),'UniformOutput',false)); % shift conic

    case 'rotate' % rotate
	
        % add
		
    case 'rigidbody'
	
        % i.e. transformation matrix
		
    case 'resize' % increase size
	
        % add
		
    case 'expand' % expand on size
	
        % add
		
    case 'normalize' % normalize peak value positive - []xN
			
		% compute center coordinate ellipse
		Xc=qvec2eshp(q);
		
		% evaluate polynomial basis
		xc=[Xc(1,:).^2
			Xc(1,:).*Xc(2,:)
			Xc(1,:).*Xc(3,:)
			Xc(2,:).^2
			Xc(2,:).*Xc(3,:)
			Xc(3,:).^2
			Xc(1,:)
			Xc(2,:)
			Xc(3,:)
			ones(size(Xc(1,:)))];
		I=sum(q.*xc,1);
		
		q=q./I;
        
    case 'invert'
        
        % compute inverse
        Q=qvec2qmat(q);
        iQ=zeros(size(Q));
        iQ(:,:,:)=cell2mat(cellfun(@(x)inv(x),num2cell(Q,[1 2]),'UniformOutput',false));
        
        % overwrite
        q=qmat2qvec(iQ);
        
end

end

