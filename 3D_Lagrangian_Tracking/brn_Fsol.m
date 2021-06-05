function Fsol = brn_Fsol(Cpln,Tlink,Fsol,n)
%brn_Fsol Branch feasible [camera] tracks, connect track segemnt in image 
%   plane to impose a overall time scale.
%   
%   Input:
%       Cpln    Camera plane ellipse idenitifications
%       Tlink   Temporal links between ellipse idenitifications
%       Fsol	Feasible solution
%       n       Current frame.
%   
%   Output:
%       Fsol    Feasible solution with new branch appendages.
%   
%   Working principle:
%       we branch tracks by 
%           1) temporal links, these act as a forcing to current feasible
%               tracks
%           2) fitting a polynomial curve to camera data and extrapolate
%           3) fitting a polynomial curve to object data and extrapolate
%               (2) and (3) act as a best prediction to also jump over 
%               a number of (ctrl.prd) empty frames.
%   
%   Note: Here kalman filtering or a RNN could be used to further 
%         improve the prediction step.
%   

%% Get globals
global prop ctrl plotting Kmat Pmat % folder date rec Dmap Dpar 

%% Initiate indexing feasible tracks

% maximum feasioble index
if isempty(Fsol)
    imaxfl=size(Fsol,2);
else
    imaxfl=max(Fsol(1,:));
end

%% Branch new feasible tracks

%-- #) Loop cameras

% loop cameras
for c=1:prop.res(4)
    
    %-- #) Begin message
    
    % keep sets complete in previous focality
    disp(['Branch new feasible tracks camera ',num2str(c),' frame ',num2str(n)])

    %-- #) Load camera data
    
    % get camera feasible tracks
    fsol=Fsol( : , ismember( Fsol(1,:) , Fsol(1, ...
        Fsol(5,:)>=( n - ctrl.prd ) & Fsol(5,:)<n ) ) ... % part of from where we want to extend
        & ~ismember( Fsol(1,:),Fsol(1, Fsol(5,:)>=n & ~isnan(sum(Fsol,1)) ) ) ); % not part of current or greater frame if not as nan currently (by code)
%     fsol=Fsol( : , ismember( Fsol(1,:) , Fsol(1, ...
%         Fsol(4,:)==c & Fsol(5,:)>=( n - ctrl.prd ) & Fsol(5,:)<n ) ) ... % part of from where we want to extend
%         & ~ismember( Fsol(1,:),Fsol(1, Fsol(5,:)>=n & ~isnan(sum(Fsol,1)) ) ) ); % not part of current or greater frame if not as nan currently (by code)
    
    % get camera identifications in current frame
    cpln=Cpln( : , Cpln(2,:)==c &  Cpln(3,:)==n ); %
    
    % get temporal links
    tlink=Tlink(:, Tlink(3,:)>=( n-ctrl.prd ) & Tlink(5,:)==n ...
        & ismember(Tlink(4,:),fsol(3,:)) & ismember(Tlink(6,:),cpln(1,:)) ); % camera
    
    %-- #) Open for reassessment
    
    % feasible solution (for reindixing when tracks split and merge)
    Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));
    
    %-- #) Define sets and adjacency matrices
    
    % camera set
    g_cam=fsol(4,:)==c;
    
    % prediction set
    g_prd=fsol(5,:) >= (n - ctrl.prd); % - ctrl.prd
    
    % track solution mapping
    [ufl,~,fl]=unique(fsol(1, : ));
    G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
    
    %figure; spy(G_sol)
    %figure; plot(sum(G_sol,2),'.')
    
    % track solution mapping
    [ufl,~,fl]=unique(fsol(1, g_prd ));
    [ucl,~,cl]=unique(fsol(3, g_prd )); 
    F_sol=sparse(cl,fl,ones(size(fl)),length(ucl),length(ufl));
    
    %figure; spy(F_sol)
    
    % define camera sorting
    [ucp,~,icp]=unique(cpln(1,:));
    C_pln=sparse(icp,1:length(icp),ones(size(icp)),length(ucp),length(icp));
    
    %figure; spy(C_pln)
    
    % define temporal connection
    [ult,~,ilt]=unique(tlink(4,:)); % unique links in tracked part
    iult=find(ismember(ucl,ult));
    ilt=iult(ilt);
    [uln,~,iln]=unique(tlink(6,:)); % unique links in new part
    iuln=find(ismember(ucp,uln));
    iln=iuln(iln);
    T_lnk=sparse(iln,ilt,ones(size(ilt)),length(ucp),length(ucl)); % temporal asymetry allowed
    
    %figure; spy(T_lnk)
    
    %-- #) Temporal link adjaceny
    T_lnk = zeros(size(T_lnk));
    % Temporal link adjacency
    Lij=C_pln'*T_lnk*F_sol;%(g_cpln,:);
    
    %figure; spy(Lij)
    
    %-- #) Get camera trajectory data
    
    % find max frame camera
    dum=find(g_cam);
    [n_max,j_max]=max(G_sol(:,g_cam)*spdiags(fsol(5,g_cam)',0,nnz(g_cam),nnz(g_cam)),[],2); % maximum frame in camera selected
    g_max=n_max>0.5;
    j_max=dum(j_max(n_max>0.5));
    
    % fit stencil camera
    n_max=n_max'*G_sol; % map back to camera data
    g_fit=fsol(5,:)>= n_max-ctrl.tfit & fsol(5,:)<n; % select frames to fit
    
    % fit trajectory to camera data
    cdat=zeros(4,0);
    xdat=cvec2eshp(fsol(6:11,g_fit & g_cam));
    for t=0:ctrl.tres%-ctrl.tres:ctrl.tres
        cdat=cat(2,cdat,[fsol([1 5],g_fit & g_cam)+[0 t/(2*ctrl.tres+1)]'
            xdat+fsol(12:13,g_fit & g_cam)*t/(2*ctrl.tres+1)]);
    end
    
    % Define trajectories(1:2,ctrk(2,:)<=n)
    ctraj=polytraj(cdat,ctrl.ord); % incl vel, fit full stencil
    
    %-- #) Find track branching adjacency
    
    % camera conic in new camera data
    cdat_i=contrans(cpln(4:9,:),'resize',ctrl.dspl);
    
    % camera conics track position data
    cdat_j=contrans(fsol(6:11,j_max),'resize',ctrl.dspl);
    
    % positions in new camera data
    xdat_i=cvec2eshp(cdat_i); %
    
    % positions in track position data
    xdat_j=cvec2eshp(cdat_j); %
    
    % predicted position data
    ydat_j=polyeval(ctraj,[fsol(1,j_max) ; n*ones(1,length(j_max))] ,0); %
    
    % displace the tracked conic on a predicted path
    cdat_j=contrans(cdat_j,ydat_j-xdat_j,'displace'); %
    
    % compute distance
    Dji=ellipsedist(cdat_j,xdat_i,'matrix'); % compute distance on expanded conic and new position only
    
    % compute distance
    Dij=ellipsedist(cdat_i,ydat_j,'matrix'); % compute distance on expanded conic and new position only
    
    % Ellipse adjacency
    dum=find(g_max);
    [ti,tj]=find(C_pln'*double(Dij<=1 | Dji'<=1));
    tj=dum(tj);
    Tij=sparse(ti,tj,ones(size(tj)),length(ucp),length(ufl));
    
    %figure; spy(Tij)
    
    %-- #) get object trajectories
    
    % complete stencil camera
    [n_max,j_max]=max(G_sol*spdiags(fsol(5,:)',0,size(fsol,2),size(fsol,2)),[],2); % get maximum frame in object space (including nans appended)
    n_max=n_max'*G_sol; % map to data
    g_fit=fsol(5,:)>= n_max-ctrl.tfit & fsol(5,:)<n; % select fit data (and don't include possible nans)
    
    % fit data
    fdat=zeros(5,0);
    for t=0:ctrl.tres%-ctrl.tres:ctrl.tres
        fdat=cat(2,fdat,[fsol([1 5],g_fit)+[0 t/(2*ctrl.tres+1)]'
            fsol(14:16,g_fit)+fsol(17:19,g_fit)*t/(2*ctrl.tres+1)]);
    end
    
    % Define trajectories
    ftraj=polytraj(fdat,ctrl.ord); % 3d
    
    %-- #) branch via object trajectory
    
    % expand ellipse
    Dreg=prop.ecal(c)/sqrt(abs(det(Kmat{c}/Kmat{c}(end))));
    
    % resize ellipse 
    cdat_i=contrans(cpln(4:9,:),ctrl.dspr,'resize');
    cdat_i=contrans(cdat_i,Dreg,'expand');
    
    % camera conics track position data
    cdat_j=contrans(fsol(6:11,j_max),'resize',ctrl.dspr);
    cdat_j=contrans(cdat_j,Dreg,'expand');
    
    % positions in track position data
    xdat_j=cvec2eshp(cdat_j); %
    
    % predicted position data
    Xobj_j=polyeval(ftraj,[ufl ; n*ones(1,length(ufl))] ,0); %fobj(1,sj)
    ydat_j=homc2inhc(Pmat{c}*inhc2homc(Xobj_j));
    
    % displace the tracked conic on a predicted path
    cdat_j=contrans(cdat_j,ydat_j-xdat_j,'displace'); %
    
    % compute distance
    Rij=ellipsedist(cdat_i,ydat_j,'matrix'); % compute distance on expanded conic and new position only
    
    Rji=ellipsedist(cdat_j,xdat_i,'matrix');
    
    % segment
    Pij=C_pln'*double(Rij<=1 | Rji'<=1 );%(:,dum)
    
    %figure; spy(Pij)
    
    %-- #) Find track branch adjacency
    
    % segment valid distance
    A=( Lij | Tij | Pij ); % 
    
    %figure; spy( A )
    %figure; spy( Lij )
    %figure; spy( Tij )
    %figure; spy( Pij )
        
    % define pairing
    [ai,aj]=find( A );
    
	% define unconnected and new feasible solutions
    g_unc=sum(A,1)==0;%true(1,size(A,2));%
    
    %-- #) keep track
    
    % keep tracks
    f_sol=fsol(:,ismember(fsol(1,:),ufl(g_unc))); 
    
    %-- #) Write splitting tracks
    
    % define re-indexing branching tracks
    ifl_brn=imaxfl+(1:length(ai)); % can be done reusing a subset of present indices
    imaxfl=imaxfl+length(ai);
    
    % find mapping to copy data
    [gi,gj]=find(G_sol(aj,:));
    
    %figure; spy(G_sol(aj,:))
    
    % pairing temporal index corresponces
    pf_sol=zeros( size(fsol,1) , length(gj) );
    pf_sol(1,:)=ifl_brn(gi); % 
    for i=2:size(pf_sol,1)
        pf_sol(i,:)=fsol(i,gj); 
    end
    
    %figure; plot(pf_sol(3,:),'.')
    
    % initiate branche positions
    [~,gs,~]=unique(gi);
    pf_new=[ifl_brn(gi(gs))
        fsol(2,gj(gs))
        cpln(1:3,ai(gi(gs)))
        cpln(4:9,ai(gi(gs)))
        cpln(10:11,ai(gi(gs)))
        nan(7,nnz(gs))]; % intiate nan is important read property to identify track extension
    
    %-- #) Append
    
    % write dataf_sol,
    fsol=cat(2,f_sol,pf_sol,pf_new); 
    
    %-- #) write new tracking data

    % write data
    Fsol=cat(2,Fsol,fsol); 
    
    %-- #) End message
    
    % message
    disp(['Branched ',num2str(nnz(pf_new(4,:)==c)),' out of ',...
        num2str(length(unique(fsol(1,fsol(4,:)==c))))])
    
end

%% plotting
if strcmp(plotting,'on')
    
    % frame
    N=Fsol(5,:)==n;
    
    % figure cost function
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(20,N),'.')
    axis tight
    axis equal
    view(2)
    title(['Branched solution, cost value frame ',num2str(n)])
    colorbar
    drawnow
    
    % figure indexing
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(1,N),'.')
    colormap lines
    axis tight
    axis equal
    view(2)
    title(['Branched solution, indexing frame ',num2str(n)])
    drawnow
    
end

end

