function [ c ] = contrans( c,T,typ )
%contrans Transform a Conic to new [reference] grid
%   input 1: conic vector 6xN
%   input 2: transformation input see code below
%   input 3: type of transform see code below

% correct input
if nargin<3
    typ='transmat';
end

% assemble transformation matrix
switch typ 
    case 'transmat' % transformation matrix - 3x3xN
        
        % make conic matrix
        C=cvec2cmat(c);
        
        % make same size if needed
        if size(T,3)<size(C,3)
            T=repmat(T,1,1,size(C,3));
        end
        if size(C,3)<size(T,3)
            C=repmat(C,1,1,size(T,3));
        end
        
        % transform conic
        if ~isempty(C) && ~isempty(T)
            C=cell2mat(cellfun(@(C,T)inv(T)'*C*inv(T),...
                num2cell(C,[1 2]),num2cell(T,[1 2]),'UniformOutput',false));
        end
        
        % write as vector
        c=cmat2cvec(C);
        
    case 'displace' % displace new position - 2xN
        
        % make conic matrix
        C=cvec2cmat(c);
        
        % compose transformation matrix
        T=reshape(T,size(T,1),1,size(T,2));
        T=cell2mat(cellfun(@(x)[1 0 x(1)
            0 1 x(2)
            0 0 1],num2cell(T,[1 2]),'UniformOutput',false)); % shift conic
        
        % make same size if needed
        if size(T,3)<size(C,3)
            T=repmat(T,1,1,size(C,3));
        end
        if size(C,3)<size(T,3)
            C=repmat(C,1,1,size(T,3));
        end
        
        % transform conic
        if ~isempty(C) && ~isempty(T)
            C=cell2mat(cellfun(@(C,T)inv(T)'*C*inv(T),...
                num2cell(C,[1 2]),num2cell(T,[1 2]),'UniformOutput',false));
        end
        
        % write as vector
        c=cmat2cvec(C);
        
    case 'rotate' % rotate - 1xN
	
        % add
		
    case 'rigidbody'
	
        % i.e. transformation matrix
		
    case 'resize' % increase size - 1xN
        
        % find level set
        [H,g,s]=cvec2lset(c);
        
        % find midpoint
        x=cvec2eshp(c);
        x=reshape(x,2,1,[]);
        
        % compose transformation matrix
        T=reshape(T,size(T,1),1,size(T,2));
        
        % make same size if needed
        if size(T,3)<size(H,3)
            T=repmat(T,1,1,size(H,3));
        end
        if size(H,3)<size(T,3)
            H=repmat(H,1,1,size(T,3));
        end
        
        % change scale, therefore and zero-levelset
        s=T.^2.*s - (1-T.^2)./2.*sum(g.*x,1);
        
        % write vector
        c=lset2cvec( H,g,s );
        
    case 'expand' % expand on size - 1/2xN
        
        % find level set
        [H,g,s]=cvec2lset(c);
        
        % find midpoint
        x=cvec2eshp(c);
        x=reshape(x,2,1,[]);
        
        % correct transformation input
        if size(T,1)==1
            T=repmat(T,2,1);
        end
        
        % compose transformation matrix
        T=reshape(T,size(T,1),1,size(T,2));
        
        % make same size if needed
        if size(T,3)<size(H,3)
            T=repmat(T,1,1,size(H,3));
        end
        if size(H,3)<size(T,3)
            H=repmat(H,1,1,size(T,3));
        end
        
        % regularize
        for i=1:size(H,3)
            [R,D]=schur(H(:,:,i)/2); % decompose C2x2
            d=s(:,:,i)+1/2*g(:,:,i)'*x(:,:,i);
            D=-D/d;
            D=inv((sqrt(inv(D))+diag(T(:,:,i))).^2);
            iH=inv([R x(:,:,i) ; zeros(1,2) 1]);
            c(:,i)=-cmat2cvec(iH'*[D zeros(2,1) ; zeros(1,2) -1]*iH);
        end
        
    case 'normalize' % normalize peak value positive - []xN
        
        x=cvec2eshp(c);
        
        x=reshape(x,size(x,1),size(c,2),size(c,3));
        
        I=c(1,:,:).*x(1,:,:).^2 ...
            +c(2,:,:).*x(1,:,:).*x(2,:,:) ...
            +c(3,:,:).*x(2,:,:).^2 ...
            +c(4,:,:).*x(1,:,:) ...
            +c(5,:,:).*x(2,:,:) ...
            +c(6,:,:) ; 
        
        c=c./I;
        
    case 'invert'
        
        % compute inverse conic
        C=cvec2cmat(reshape(c,6,[])); % not all be valid
        iC=zeros(size(C));
        iC(:,:,:)=cell2mat(cellfun(@(x)inv(x),num2cell(C,[1 2]),'UniformOutput',false));
        
        % overwrite
        c=cmat2cvec(iC);
                
end

end

