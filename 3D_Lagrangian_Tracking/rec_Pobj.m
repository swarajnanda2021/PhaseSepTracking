function [Plink,Pobj] = rec_Pobj(Fsol,frm_proc)
%rec_Pobj Reconstruct physical object data.
%
%   Input:
%       Cpln    Camera plane ellipse idenitifications
%       Fsol    Feasible solution.
%   
%   Output:
%       Pobj    Newly computed physical object reconstruction. (Here we
%               choose to retriangulate from Fsol to prevent overprocessing
%               from smoothing and differentiation / kalman filtering in
%               feasible solutions.)
%
%   Note: better spherical regularization procedure for the quadrics
%       introduce weighted vector and perhaps use cameras rotated from
%       previous frames.
%   

%% Get global variables
global prop ctrl Pmat Kmat

%% initiate data

% Physical object data
Plink=zeros(19,0); % ( track index | indentification index | camera index | frame number | camera conic | camera velocity | object position | object velocity | res )
Pobj=zeros(17,0); % ( track index | frame number | object quadric | object velocity | peak intensity | res )

%% Reconstruct
disp('Reconstruct matched identifications')

% loop over the matched and indexed Plinks
for n=frm_proc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% get data from matching and files %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % optimized feasible solutions
    fsol=Fsol(:, Fsol(5,:)==n ); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get adjacencie mappings %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % define feasible solution adjacency
    [ufl,~,fl]=unique(fsol(2, : ));%ismember(fsol(1,:),ufl)
    G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
    
    %figure; spy(G_sol)
    
    f_sum=full(sum(G_sol,2))';
    
    %figure; spy(F_sol)
    %figure; plot(f_sum,'.')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% solve midpoint triangulation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initiate triangulation data
    xdat=cell(1,prop.res(4));
    
    % Assemble triangulation data
    for c=1:prop.res(4)
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % Get midpoints for triangulation
        x=cvec2eshp(fsol(6:11,g_cam ));
        
        % adjacencie camera data and solution
        A=G_sol(:,g_cam );
        
        % find
        [ai,aj]=find(A);
        
        % initiate and write
        xdat{c}=nan(3,length(ufl));
        xdat{c}(:,ai)=inhc2homc(x(:,aj));
        
    end
    
    % Triangulate
    X=objtriang(xdat,Pmat);
    
    %figure; plot3(X(1,:),X(2,:),X(3,:),'b.'); axis equal; view(2); drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% displaced midpoint triangulation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initiate triangulation data
    xdat=cell(1,prop.res(4));
    
    % Assemble triangulation data
    for c=1:prop.res(4)
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % Get midpoints for triangulation
        x=cvec2eshp(fsol(6:11,g_cam))+fsol(12:13,g_cam);
        
        % adjacencie camera data and solution
        A=G_sol(:,g_cam);
        
        % find
        [ai,aj]=find(A);
        
        % initiate and write
        xdat{c}=nan(3,length(ufl));
        xdat{c}(:,ai)=inhc2homc(x(:,aj));
        
    end
    
    % Triangulate
    Xacc=objtriang(xdat,Pmat); % [,~,r]
    
    %figure; plot3(X(1,:),X(2,:),X(3,:),'b.'); axis equal; view(2); drawnow
    %hold on; plot3(Xacc(1,:),Xacc(2,:),Xacc(3,:),'r.'); axis equal; view(2); drawnow
    
    % complute velocity vector
    V=Xacc-X; 
    
    %figure; plot3(X(1,:),X(2,:),X(3,:),'b.'); axis equal; view(2); drawnow
    %hold on; quiver3(X(1,:),X(2,:),X(3,:),V(1,:),V(2,:),V(3,:),0,'r')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% quadric reconstruction %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initiate conic data
    cdat=cell(1,prop.res(4));
    
    % loop cameras
    for c=1:prop.res(4)
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % Get midpoints for triangulation
        con=fsol(6:11,g_cam);
        
        %Xe=reshape(cvec2pnts(con),2,[]);
        %figure; plot(Xe(1,:),Xe(2,:),'.r')
        
        % adjacencie camera data and solution
        A=G_sol(:,g_cam);
        
        % find
        [ai,aj]=find(A);
        
        % initiate conic data
        cdat{c}=nan(6,length(ufl));
        
        % porject midpoint triangulation
        xp=homc2inhc(Pmat{c}*inhc2homc(X));
        
        %figure; plot(xp(1,:),xp(2,:),'.')
        
        % shift conics
        con=contrans(con(:,aj),xp(:,ai)-cvec2eshp(con(:,aj)),'displace');
        
        %Xe=reshape(cvec2pnts(con),2,[]);
        %figure; plot(Xe(1,:),Xe(2,:),'.r')
        
        % write conic data
        cdat{c}(:,ai)=con;
        
    end
    
    % constrained triangulation
    Qt=quadrec(cdat,Pmat,X,'spher');
	
	% normalize quadric
    Qt=quadtrans(Qt,[],'normalize');
    
    %Xq=reshape(qvec2pnts(Qt),3,[]);
    %figure; plot3(Xq(1,:),Xq(2,:),Xq(3,:),'r.'); axis tight; axis equal; view(2); drawnow
    %xlim([-3 3]); ylim([11 17]); zlim([-5 1])
    
    % unconstrained triangulation
    Qf=quadrec(cdat,Pmat);
    
    % positiveness 
    [~,g_val]=qvec2pnts(Qf);
    Qf(:,~g_val)=nan;
    
	% normalize quadric
	Qf=quadtrans(Qf,[],'normalize');
    
    %Xq=reshape(qvec2pnts(Qf),3,[]);
    %figure; plot3(Xq(1,:),Xq(2,:),Xq(3,:),'r.'); axis tight; axis equal; view(2); drawnow
    %xlim([-3 3]); ylim([11 17]); zlim([-5 1])
    
    % average level sets
    Q=nansum(cat(3,f_sum.^(-1).*Qt,(1-f_sum.^(-1)).*Qf),3)/2; % for now
    
	% normalize quadric
	Q=quadtrans(Q,[],'normalize'); % consistency
    
    %Xq=reshape(qvec2pnts(Q),3,[]);
    %figure; plot3(Xq(1,:),Xq(2,:),Xq(3,:),'r.'); axis tight; axis equal; view(2); drawnow
    %xlim([-3 3]); ylim([11 17]); zlim([-5 1])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% average residual rationalize by body length %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % midpoint from quadric
    X=qvec2eshp(Q);
    
    % initiate disparity
    D=zeros(1,length(ufl));
    
    % loop cameras
    for c=1:prop.res(4)
        
        % regularization for camera accuracy
        Dreg=prop.ecal(c)/sqrt(abs(det(Kmat{c}/Kmat{c}(end))));
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % expand displacement allowance
        cdat=contrans(fsol(6:11,g_cam),ctrl.dspl,'resize'); %+ctrl.dspl/(1+2*ctrl.tres) the higher the time resolution the higher the certainty of the position
        
        % expand disparity allowance
        cdat=contrans(cdat,ctrl.dspr,'resize'); % the higher the time resolution the higher the certainty of the position
        
        % camera uncertainty
        cdat=contrans(cdat,Dreg+sqrt(sum(fsol(12:13,g_cam).^2,1))/(1+2*ctrl.tres),'expand');
        
        % adjacencie camera data and solution
        A=G_sol(:,g_cam);
        
        % find
        [ai,aj]=find(A);
        
        % loop image displacement
        for t=0:ctrl.tres
            
            % camera conic along camera trajectory
            con=contrans(cdat,fsol(12:13,g_cam )*t/(2*ctrl.tres+1),'displace'); %fsol(12:13,g_cam & g_frm)*t/(2*ctrl.tres+1) fsol(6:11,g_cam & g_frm);%
            
            % object track
            pos=X+V*t/(2*ctrl.tres+1);
            
            %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
            
            % camera reference frame
            pos=Pmat{c}*inhc2homc(pos) ;
            
            %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
            
            % project coordinates
            xp=homc2inhc( pos ); % empty cell can comprimise
            
            %figure; plot(x(1,:),x(2,:),'.')
            
            % compute ellipse distance
            d=ellipsedist(con(:,aj),xp(:,ai),'vector');
            
            %figure; plot(d,'.')
            
            % write
            D(ai)=D(ai)+d/(2*ctrl.tres+1);
            
        end % t
        
    end % c
    
    % average disparity
    D=D./sum(G_sol,2)';
    
    %figure; plot(D,'.')
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% peak intensity %%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    % adjacencie camera data and solution
    [~,~,~,I]=cvec2lset(fsol(6:11,:),cvec2eshp(fsol(6:11,:)));
    
    % average intensity data for the conic from different views
    I=I*G_sol'./sum(G_sol,2)'; % sum by matrix
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% write data to file %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plink object reconstruction
    plink=fsol(2:end,:);
    
    % new object data ( index | frame number | quadric | parmetrization | peak intensity | res )
    pobj=[ufl
        n*ones(size(ufl))
        Q
        V
        I
        D];
    
    % write solution
    Plink=cat(2,Plink,plink); 
    Pobj=cat(2,Pobj,pobj); 
    
    % message
    disp(['Reconstructed ',num2str(length(ufl)),...
        ' physical objects at average disparity ',num2str(mean(D))...
        ' in frame ',num2str(n)])
    
end

end

