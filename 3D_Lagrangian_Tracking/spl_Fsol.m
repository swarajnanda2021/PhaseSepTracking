function Fsol = spl_Fsol(Fsol,n)
%spl_Fsol Split feasible solutions by triangulation
%
%   Input:
%       Fsol    Feasible solutions
%       n       Current frame.
%
%   Output:
%       Fsol 	Feasible solutions with new triangulations and updates
%               trajectory cost.
%
%   working principle:
%       Split solution that violate the reprojection error in their
%		(n-1) focal subsets to extend along the best direction
%

%% Get globals
global prop ctrl plotting Pmat Kmat

%% Initiate indexing feasible tracks

% maximum feasioble index
if isempty(Fsol)
    imaxfl=size(Fsol,2);
else
    imaxfl=max(Fsol(1,:));
end

% maximum track index
if isempty(Fsol)
    imaxtl=size(Fsol,2);
else
    imaxtl=max(Fsol(2,:));
end

%% Match nfocal sets

%-- #) End message

% keep sets complete in previous focality
disp(['Split feasible solutions frame ',num2str(n)])

%-- #) cluster focality

% related feasible solutions
fsol=Fsol(:, ismember( Fsol(1,:) , Fsol(1, Fsol(5,:)==n ) ) ) ; % & isnan(Fsol(20,:)) 

%-- #) Open for reassessment

% feasible solution (for reindixing when tracks split and merge)
Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));

%-- #) Split timeseries and reopen existing triangulation by nan

% if exist before and after then split
g_spl=ismember(fsol(1,:),fsol(1,fsol(5,:)<n)) ... 
    & ismember(fsol(1,:),fsol(1,fsol(5,:)>n));

% Past feasible solutions
fsolm=fsol(:, fsol(5,:)<=n & g_spl );

% Future and current feasible solutions
fsolp=fsol(:, fsol(5,:)>=n & g_spl );

% reindex future branches
[ufl,~,fl]=unique(fsolp(1,:));
ind_spl=imaxfl+fl;
imaxfl=imaxfl+length(ufl);
fsolp(1,:)=ind_spl; % overwrite feasible solution, track index will be kept in optimization

% reindex future branches
[~,~,ft]=unique(fsolp(2,:));
ind_spl=imaxtl+ft;
% imaxtl=imaxtl+length(uft);
fsolp(2,:)=ind_spl; % overwrite feasible solution, track index will be kept in optimization

% currently branching
fsol=fsol(:,~g_spl);

% append
fsol=cat(2,fsol,fsolm,fsolp);

% reopen to triangulate all tracks
fsol(14:20,fsol(5,:)==n)=nan;

%-- #) Split solution if triangulation does not fall in reprojection

% initiate removal
nsol=zeros(20,0);

% while loop when nans remain
while nnz(isnan(fsol))~=0 %for looop=1:4 % quick and dirty
    
    %-- #) Redefine adjacencies and sets
    
    % track solution mapping
    [ufl,~,fl]=unique(fsol(1, : ));
    G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
    
    %figure; spy(G_sol)
    
    % triangulation set
    g_tri=isnan(sum(fsol,1));
    
    % triang set
    g_ufl=sum(G_sol(:,g_tri),2)' >0.5;
    
    %-- #) Triangulation
    
    % Assemble triangulation data
    xdat=cell(1,prop.res(4));
    for c=1:prop.res(4)
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % Get midpoints for triangulation
        x=cvec2eshp(fsol(6:11,g_cam & g_tri));
        
        % adjacencie camera data and solution
        A=G_sol(g_ufl,g_cam & g_tri);
        
        % find
        [ai,aj]=find(A);
        
        % initiate and write
        xdat{c}=nan(3,nnz(g_ufl));
        xdat{c}(:,ai)=inhc2homc(x(:,aj));
        
    end
    
    % Triangulate
    X=objtriang(xdat,Pmat);
    
    %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2); colorbar; drawnow
    
    % find indices
    dum=find( g_tri );
    [gi,gj]=find(G_sol(g_ufl, g_tri ));
    gj=dum(gj);
    
    % write position data
    fsol(14:16,gj)=X(:,gi);
    
    %figure; plot3(fsol(14,:),fsol(15,:),fsol(16,:),'.'); axis equal; view(2); drawnow
    
    %-- #) Displacement
    
    % Assemble triangulation data
    xdat=cell(1,prop.res(4));
    for c=1:prop.res(4)
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % Get midpoints for triangulation
        x=cvec2eshp(fsol(6:11,g_cam & g_tri))+fsol(12:13,g_cam & g_tri);
        
        % adjacencie camera data and solution
        A=G_sol(g_ufl,g_cam & g_tri);
        
        % find
        [ai,aj]=find(A);
        
        % initiate and write
        xdat{c}=nan(3,nnz(g_ufl));
        xdat{c}(:,ai)=inhc2homc(x(:,aj));
        
    end
    
    % Triangulate
    Xacc=objtriang(xdat,Pmat); % [,~,r]
    
    %figure; plot3(Xacc(1,:),Xacc(2,:),Xacc(3,:),'.'); axis equal; view(2); colorbar; drawnow
    
    % velocity
    V=Xacc-X;
    
    % find indices
    dum=find(g_tri);
    [gi,gj]=find(G_sol(g_ufl,g_tri));
    gj=dum(gj);
    
    % write position data
    fsol(17:19,gj)=V(:,gi);
    
    %--# Exclusion by nan
    
    % exclusion by nan
    g_excl=isnan(sum(fsol(14:19,:)));
    
    % remove exclusion from trajectory data
    nsol=cat(2,nsol,fsol(:,g_excl));
    fsol=fsol(:, ~g_excl );
    
    %-- #) Redefine adjacencies
    
    % % track solution mapping
    % [ufl,~,fl]=unique(fsol(1, : ));
    % G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
    %
    % %figure; spy(G_sol)
    
    % smooth new triangulations
    g_smt=ismember(fsol(1,:),fsol(1,isnan(fsol(20,:))));
    
    % track fit stencil
    g_fit=fsol(5,:)>=( n+1-ctrl.tfit ) & fsol(5,:)<=( n-1+ctrl.tfit ) & g_smt;%fsol(5,:)==n;fsol(5,:)>=( n - ceil((ctrl.tfit-1)/2) ) & fsol(5,:)<=( n + ceil((ctrl.tfit-1)/2) );%fsol(5,:)==n;
    
    % frame set
    g_frm=fsol(5,:)==n & g_smt;%cpln(2,:)==c;
    
    %-- Filter camera and object trajectories
    
    % fit data
    fdat=zeros(5,0);
    for t=0:ctrl.tres
        fdat=cat(2,fdat,[fsol([1 5],g_fit)+[0 t/(2*ctrl.tres+1)]'
            fsol(14:16,g_fit)+fsol(17:19,g_fit)*t/(2*ctrl.tres+1)]);
    end
    
    % Define trajectories
    ftraj=polytraj(fdat,ctrl.ord); % inherently weighted by camera occurence
    
    fsol(14:16,g_frm)=polyeval(ftraj,fsol([1 5],g_frm),0);
    
    fsol(17:19,g_frm)=polyeval(ftraj,fsol([1 5],g_frm),1);
    
    %-- Redefine adjacencies
    
    % triangulation set
    g_rep=isnan(fsol(20,:));
    
    % triangulation set
    g_frm=fsol(5,:)==n;
    
    %-- #) Reprojection
    
    % define invalid part trajectories
    g_excl=false(1,size(fsol,2));
    
    % initiate cost
    fsol(20,g_rep)=0; % to recompute and update by traj fit
    
    % loop cameras
    for c=1:prop.res(4)
        
        % regularization for camera accuracy
        Dreg=prop.ecal(c)/sqrt(abs(det(Kmat{c}/Kmat{c}(end))));
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % expand displacement allowance
        cdat=contrans(fsol(6:11,g_cam & g_frm & g_rep),ctrl.dspl,'resize'); %+ctrl.dspl/(1+2*ctrl.tres) the higher the time resolution the higher the certainty of the position
        
        % expand disparity allowance
        cdat=contrans(cdat,ctrl.dspr,'resize'); % the higher the time resolution the higher the certainty of the position
        
        % camera uncertainty
        cdat=contrans(cdat,Dreg+sqrt(sum(fsol(12:13,g_cam & g_frm & g_rep).^2,1))/(1+2*ctrl.tres),'expand');
        
        % loop image displacement
        for t=0:ctrl.tres
            
            % camera conic along camera trajectory
            con=contrans(cdat,fsol(12:13,g_cam & g_frm & g_rep)*t/(1+2*ctrl.tres),'displace');
            
            % object track
            X=fsol(14:16, g_cam & g_frm & g_rep)+fsol(17:19, g_cam & g_frm & g_rep)*t/(1+2*ctrl.tres);
            
            %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
            
            % camera reference frame
            X=Pmat{c}*inhc2homc(X) ;
            
            %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
            
            % segment dof
            g_excl( g_cam & g_frm & g_rep)=g_excl( g_cam & g_frm & g_rep) | X(3,:)<ctrl.ldof(1) | X(3,:)>ctrl.ldof(2);
            
            % project coordinates
            x=homc2inhc( X ); % empty cell can comprimise
            
            %figure; plot(x(1,:),x(2,:),'.')
            
            % compute ellipse distance
            d=ellipsedist(con,x,'vector');
            
            %figure; plot(d,'.')
            
            %     % reprojection
            %     x=cvec2eshp(con);
            %     r=sqrt(sum((x(:,ai)-xp).^2,1));
            
            %     %Xe=reshape(cvec2pnts(con),2,[]);
            %     %hold on; plot(Xe(1,:),Xe(2,:),'.r')
            
            % update exclusion
            g_excl( g_cam & g_frm & g_rep)=g_excl( g_cam & g_frm & g_rep) | ( d > 1 );
            
            % write
            fsol(20, g_cam & g_frm & g_rep)=fsol(20, g_cam & g_frm & g_rep)+d/(2*ctrl.tres+1);
            
        end % t
        
    end % c
    
    %figure; plot(fsol(20,:),'.')
    
    % exclusion by feasible solution
    g_excl=ismember(fsol(1,:),fsol(1,g_excl));
    
    %-- split invalid triangulation in subtriangulation
    
    % reopen invalid solution
    fsol(14:20,g_frm & g_excl)=nan;
    
    % initiate split
    f_sol=fsol(:, ~g_excl );
    
    % loop cameras
    pf_sol=zeros(20,0);
    for c=1:prop.res(4)
        
        % camera set
        g_cam=fsol(4,:)==c;
        
        % split present in camera
        pf_split=fsol(:, g_frm & ismember(fsol(1,:),fsol(1,g_cam & g_frm & g_excl)) );
        
        % remove that camera from the split
        pf_split=pf_split(:,~(pf_split(4,:)==c));
        
        % correcsponding previous to split
        pf_prev=fsol(:,~g_frm & ismember(fsol(1,:),pf_split(1,:)));
        
        % new
        pf_new=cat(2,pf_prev,pf_split);
        
        % unique indexing in split
        [up,~,ip]=unique(pf_new(1,:));
        
        % define re-indexing branching tracks
        ifl_brn=imaxfl+ip; % can be done reusing a subset of present indices
        imaxfl=imaxfl+length(up);
        
        % overwrite indexing to split
        pf_new(1,:)=ifl_brn; %
        
        % split solution removing one camera
        pf_sol=cat(2,pf_sol,pf_new);
        
    end
    
    % remove exclusion from trajectory data
    nsol=cat(2,nsol,fsol(:,g_frm & g_excl));
    fsol=cat(2,f_sol,pf_sol);
    
    %-- #) Redefine adjacencies
    
    % track solution mapping
    [ufl,~,fl]=unique(fsol(1, : ));
    G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
    
    %figure; spy(G_sol)
    
    % frame set
    g_frm=fsol(5,:)==n;%cpln(2,:)==c;
    
    %-- #) Minimum focality
    
    % current camera cover
    cf=sum( G_sol(:,g_frm) ,2)'; % previous tracking
    
    %figure; plot(cf,'.')
    
    % define exclusion
    g_excl=g_frm & ( ismember(fsol(1,:),...
        ufl( cf<ctrl.lfoc(1) ) ) );
    
    % remove wrong triangulation as whole & fsol(5,:)==n
    % g_excl=g_excl | ( fsol(5,:)>n & ismember(fsol(1,:),fsol(1,g_excl)) ); % fsol(5,:)>=n & ismember(fsol(1,:),fsol(1,g_excl)) ;
    
    % remove exclusion from trajectory data
    nsol=cat(2,nsol,fsol(:,g_excl));
    fsol=fsol(:, ~g_excl );
    
    %-- #) Redefine adjacencies
    
%     % track solution mapping
%     [ufl,~,fl]=unique(fsol(1, : ));
%     G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
%     
%     %figure; spy(G_sol)
    
    % track fit stencil
    g_prd=fsol(5,:)>=(n-ctrl.tfit) & fsol(5,:)<=(n+ctrl.tfit) ;%fsol(5,:)==n;fsol(5,:)>=( n - ceil((ctrl.tfit-1)/2) ) & fsol(5,:)<=( n + ceil((ctrl.tfit-1)/2) );%fsol(5,:)==n;
    
    % feasible track index mapping
    [ufl,~,fl]=unique(fsol(1, : ));
    [utl,~,tl]=unique(fsol(2, : ));
    T_sol=double(sparse(tl,fl,ones(size(tl)),length(utl),length(ufl))>0.5); % else >0.5 !
    
    %figure; spy(T_sol
    
    % new camera coordinates mapping g_prd
    [ufl,~,fl]=unique(fsol(1, g_prd ));
    [ucl,~,cl]=unique(fsol([3 5], g_prd )','rows');
    F_sol=double(sparse(cl,fl,ones(size(cl)),size(ucl,1),length(ufl))>0.5);%
    
    %figure; spy(F_sol)
    
    f_sum=sum(F_sol,1);
    
    %figure; plot(f_sum,'.')
    
    %-- Remove identical track from splitting sequence
    
    % adjacent tracks
    A=T_sol'*T_sol;
    
    %figure; spy(A)
    
    % find
    [ai,aj]=find(tril(A,-1));
    ai=reshape(ai,1,[]);
    aj=reshape(aj,1,[]);
    
    a_sum=sum(F_sol(:,ai).*F_sol(:,aj),1);
    
    %figure; plot(a_sum,'.')
    %figure; plot(f_sum(ai),'.')
    %figure; plot(f_sum(aj),'.')
    
    % exclusion
    g_excl=ismember(fsol(1,:),ufl( unique(ai(a_sum==f_sum(ai) & a_sum==f_sum(aj))) )); % ismember(fsol(1,:),ufl( unique(ai(ak==f_sum(ai) & ak==f_sum(aj) ))));
    
    % remove exclusion from trajectory data
    nsol=cat(2,nsol,fsol(:,g_excl));
    fsol=fsol(:, ~g_excl );
    
end % ~isnan

%-- #) write

Fsol=cat(2,Fsol,fsol);

%-- #) End message

% keep sets complete in previous focality
disp(['Splitted ',num2str(length(unique(nsol(1,:)))),...
    ' from remaining ',num2str(length(unique(fsol(1,:)))),' feasible solutions'])

%% plot results
if strcmp(plotting,'on')
    
    % frame
    N=Fsol(5,:)==n;
    
    % figure cost function
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(20,N),'.')
    axis tight
    axis equal
    view(2)
    title(['Evaluated solution, cost value frame ',num2str(n)])
    colorbar
    drawnow
    
    % figure indexing
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(1,N),'.')
    colormap lines
    axis tight
    axis equal
    view(2)
    title(['Evaluated solution, indexing frame ',num2str(n)])
    drawnow
    
end

end

