function [ tdat ] = polyeval( coef, dat, der )
%polyeval Polynomial trajectory evaluation
%   Detailed explanation goes here, can be optimize and further vectorized
%   Work with 2 and 3 and N D

% correct input
if nargin<3
    der=0;
end

% define trajectory data without indexing (assume same, no need to copy)
tdat=zeros(size(coef,1)-2,size(dat,2));

% loop time
for t=unique(dat(2,:))
    % time in trajectory
    T=dat(2,:)==t;
    
    % loop coefficients
    for c=unique(coef(2,:))
        % number in coefficient
        C=coef(2,:)==c;
        
        % order coefficient
        k=c-1;
        
        % if term exists by differentiation
        if k>=der
            
            % index relations
            [ui,~,ki]=unique([coef(1,C) dat(1,T)]);
            ic=ki(1:nnz(C));
            it=ki(nnz(C)+(1:nnz(T)));
            I=sparse(ic,(1:length(ic))',ones(size(ic)),length(ui),length(ic));
            J=sparse(it,(1:length(it))',ones(size(it)),length(ui),length(it));
            A=double(I'*J > 0.5);
            
            % differentiation weight
            w=factorial(k)/factorial(k-der);
            
            % append to trajectory information
            tdat(:,T)=tdat(:,T)+w.*(coef(3:end,C)*A).*( dat(2,T).^(k-der) );
            
        end
    end
end

end

