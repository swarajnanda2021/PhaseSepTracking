function [ x , F , gF , hF ] = optobjf( wObjF, ObjF, gObjF, hObjF , x0 , dat , tol , meth , mess, alp, amp)
%optobjf Optimize objective function
%   find minimal objective function value in global assemblage
%   unknown must be unnumbered, knowns must be numbered in global system
%
%   To be added: Global optimization strategies,
%       Learning [amp & tau] can be computed a diagonal weight matrix from
%           hessian connectivity improving the global solution
%       Ant colony optimization
%       Amoeboid optimiation
%       Optimization using a color function

% correct input
if nargin <7
    tol=10^(-3);
end
if  nargin <8
    meth='newton'; % 
end
if nargin<9
    mess='disp';
end
if nargin<10
    alp=1; % initial condition learning variable
end
if nargin<11
    amp=2; % "Mohr's Law" learning amplification / powerlaw learning vs power law convergence
end

% initiate
x=x0; % solution vector
switch meth % optimization problem
    case 'gradient'
        
        [ F, gF ] = assoptp( x , dat , wObjF, ObjF , gObjF );
        hF=sparse(length(gF),length(gF));
        
    otherwise % newton[-based] BFGS ... congrad matr.
        
        if ~isempty(hObjF) % full initial problem
            
            [ F, gF , hF] = assoptp( x , dat , wObjF, ObjF , gObjF , hObjF );
            
        else % BFGS ... without initial hessian
            
            [ F, gF ] = assoptp( x , dat , wObjF, ObjF , gObjF );
            
            hF=speye(length(gF));
            
        end
        
end
tau=1/sqrt(dot(gF,gF)); % sensitivity time scale
err=inf; % error norm
iter=0; % iteration count

% start message
if strcmp(mess,'disp')
    disp(['iter. = ', num2str(iter) ,' [#]',', res. = ',num2str(F),' [R^2], err. = ',num2str(sqrt(F/(size(dat,2)-2))),' <R>'])
end

% start procedure
while err>tol && F>eps % iter nothing less than machine precision
    
    % initiate linearized path search
    dt=alp*tau; % current step
    if alp~=inf % 'line/curve search'
        xup=x-(speye(length(gF))+dt*hF)\(dt*gF);
    else % direct newton avoid inf/inf above
        xup=x-hF\gF;
    end
    Fup = assoptp( xup , dat , wObjF, ObjF );
    
    % linearized path search
    if amp~=1 || alp~=inf % then learning is able
        
        % first force Fup below current objective
        while Fup>F 
            alp=alp/amp; % decrease
            dt=alp*tau; % current step
            xup=x-(speye(length(gF))+dt*hF)\(dt*gF);
            Fup = assoptp( xup , dat , wObjF, ObjF );
        end
        
        dt=alp/amp*tau; % decreased step
        xdec=x-(speye(length(gF))+dt*hF)\(dt*gF);
        Fdec = assoptp( xdec , dat , wObjF, ObjF );
        
        dt=alp*amp*tau; % increased step
        xinc=x-(speye(length(gF))+dt*hF)\(dt*gF);
        Finc = assoptp( xinc , dat , wObjF, ObjF );
        
        while ( (Fup>Finc) || (Fup>Fdec) ) %... % descend conditions
            
            % && alp<=1/tol % learning upper bound && alp>eps
            
            % in sharp elongated shapes, dont force along gradient zigzag
            
            if Fup>Fdec % && alp/amp>eps % decrease step size
                
                alp=alp/amp;
                xup=xdec;
                Fup=Fdec;
                
                dt=alp/amp*tau;
                xdec=x-(speye(length(gF))+dt*hF)\(dt*gF);
                Fdec = assoptp( xdec , dat , wObjF, ObjF ); % info at that point
                %         else
                %             Fdec=Fup; % force execute
                
            end
            
            if Fup>Finc % && alp*amp<=1/tol % increase step size
                
                alp=alp*amp;
                xup=xinc;
                Fup=Finc;
                
                dt=alp*amp*tau;
                xinc=x-(speye(length(gF))+dt*hF)\(dt*gF);
                Finc = assoptp( xinc , dat , wObjF, ObjF ); % info at that point
                %         else
                %             Finc=Fup; % force execute
                
            end
            
        end
        
    end
    
    % select method for updates
    switch meth
        case 'gradient'
            
            % evaluate
            [ ~, gFup ] = assoptp( xup , dat , wObjF, [] , gObjF , hObjF ); % info at that point
            
            % empty hessian
            hFup=hF;
            
        case 'newton' % first ord stepping
            
            % evaluate
            [ ~, gFup, hFup ] = assoptp( xup , dat , wObjF, [] , gObjF , hObjF ); % info at that point
            
        case 'BFGS'
            
            % evaluate
            [ ~, gFup ] = assoptp( xup , dat , wObjF, [] , gObjF ); % info at that point
            
            % solution change
            s=xup-x;
            
            % gradient change
            y=gFup-gF;
            
            % update hessian
            hFup=hF+y*y'/(y'*s)-hF*s*(hF*s)'/(s'*hF*s);
            
    end
        
    % error criterium
    err=(F-Fup)/F;%norm(v,inf); % vanishing gradient: can we sign update while conv?
    if err>=0%-tol
        
        % iterate
        iter=iter+1; % iterations
        
        % message
        if strcmp(mess,'disp')
            disp(['iter. = ', num2str(iter) ,' [#], alp. = ',num2str(alp),', decr. = ',num2str(err),' / ',num2str(tol),' [%], res. = ',num2str(Fup),' [R^2], err. = ',num2str(sqrt(Fup/(size(dat,2)-2))),' <R>'])
        end
        
        % updates
        x=xup; % solution vector
        F=Fup; % ojective value
        gF=gFup; % gradient vector
        hF=hFup; % hessian matrix
        tau=1/sqrt(dot(gF,gF)); % time scale
        
    else
        
        % end message
        if strcmp(mess,'disp')
            disp(['terminated: blow-up @',' alp. = ',num2str(alp),', res. = ',num2str(Fup),' [R^2], err. = ',num2str(sqrt(Fup/(size(dat,2)-2))),' <R>'])
        end
        
    end
end

end