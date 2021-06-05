function residual = inv_distortn(params,dewarped,warped)

% A function that takes the dewarped coordinates and warps them back to the
% original pixel coordinates, then, computes the residual between the
% warped coordinates and the 


% Take cubic polynomial
x = dewarped(1,:)';
y = dewarped(2,:)';



P = [ones(size(x)) x y x.*y x.^2 y.^2 ...
         (x.^2).*y x.^3 (y.^2).*x ...
        y.^3];


    
warped_est = [P*params(1:10,1) P*params(11:20,1)]' ;



residual = (norm((warped_est - warped)'));





end

