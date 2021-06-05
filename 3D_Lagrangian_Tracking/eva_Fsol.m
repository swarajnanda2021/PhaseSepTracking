function Fsol = eva_Fsol(Fsol,n)
%eva_Fsol Evaluate feasible solutions.
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
%       1) Filter temporal coherence camera plane
%       2) Triangulate solution [deviates from matched sets triang]
%       3) Filter and correct trajectories, and compute trajectory error
%       4) Filter spurious triangulation by reprojection.
%   
%   Note: Here kalman filtering would improve the trajectorie corrections
%   

%% Get globals
global prop ctrl plotting Pmat Kmat

%% Match nfocal sets

%-- #) End message

% keep sets complete in previous focality
disp(['Evaluate feasible solutions frame ',num2str(n)])

%-- #) cluster focality

% related feasible solutions
fsol=Fsol(:, ismember( Fsol(1,:) , Fsol(1, Fsol(5,:)==n ) ) ) ;

%-- #) Open for reassessment

% feasible solution (for reindixing when tracks split and merge)
Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));

%-- #) Define adjacencies and sets

% track fit stencil
g_fit=fsol(5,:)>=( n+1-ctrl.tfit ) & fsol(5,:)<=( n );%fsol(5,:)==n;fsol(5,:)>=( n - ceil((ctrl.tfit-1)/2) ) & fsol(5,:)<=( n + ceil((ctrl.tfit-1)/2) );%fsol(5,:)==n;

% frame set
g_frm=fsol(5,:)==n ;%cpln(2,:)==c;

%-- #) Filter camera trajectories

% define invalid part trajectories
g_excl=false(1,size(fsol,2));

% loop cameras
for c=1:prop.res(4)
    
    % camera set
    g_cam=fsol(4,:)==c ;
    
    % fit data
    cdat=zeros(4,0);
    xdat=cvec2eshp(fsol(6:11, g_cam & g_fit ));
    for t=0:ctrl.tres%-ctrl.tres:ctrl.tres
        cdat=cat(2,cdat,[fsol([1 5], g_cam & g_fit )+[0 t/(2*ctrl.tres+1)]'
            xdat+fsol([12 13], g_cam & g_fit )*t/(2*ctrl.tres+1)]); 
    end
    
    %figure; plot(cdat(3,:),cdat(4,:),'.')
    
    % fit trajectories
    ctraj=polytraj(cdat,ctrl.ord);
    
    % resize to displacement
    cdat=contrans(fsol(6:11, g_cam & g_frm),ctrl.dspl,'resize');%+ctrl.dspr
    
    % loop between frame
    for t=0:ctrl.tres%-ctrl.tres:ctrl.tres
        
        % displace conic
        con=contrans(cdat, fsol(12:13, g_cam & g_frm)*t/(2*ctrl.tres+1) ,'displace');
        
        % evaluate trajectory
        ydat=polyeval(ctraj,fsol([1 5], g_cam & g_frm)+[0 t/(2*ctrl.tres+1)]',0); % null
        
        %     xdat=cvec2eshp(cdat);
        %     figure; plot(ydat(1,:),ydat(2,:),'.')
        %     hold on; plot(xdat(1,:),xdat(2,:),'.'); drawnow
        
        % compute vector ellipse distance
        d=ellipsedist(con,ydat,'vector'); % coherence
        
        %figure; plot(d,'.'); drawnow
        
        % segment disparity
        g_excl( g_cam & g_frm )=d > 1;
        
    end % t
    
end % c

% also select remaining appendage of the track to remove
g_excl=g_excl | ( fsol(5,:)>n & ismember(fsol(1,:),fsol(1,g_excl)) ); %   

% remove exclusion from trajectory data
nsol=fsol(:,g_excl);
fsol=fsol(:, ~g_excl ); 

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

% remove nans
nsol=cat(2,nsol,fsol(:,isnan(sum(fsol(1:19,:)))));%fsol(:,isnan(sum(fsol(1:19,:))));%
fsol=fsol(:,~isnan(sum(fsol(1:19,:))));

%-- #) Redefine adjacencies

% % track solution mapping
% [ufl,~,fl]=unique(fsol(1, : ));
% G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
% 
% %figure; spy(G_sol)

% track fit stencil
g_fit=fsol(5,:)>=( n+1-ctrl.tfit ) & fsol(5,:)<=( n );%fsol(5,:)==n;fsol(5,:)>=( n - ceil((ctrl.tfit-1)/2) ) & fsol(5,:)<=( n + ceil((ctrl.tfit-1)/2) );%fsol(5,:)==n;

% frame set
g_frm=fsol(5,:)==n ;%cpln(2,:)==c;

%-- Filter camera and object trajectories

% fit data
fdat=zeros(5,0);
for t=-ctrl.tres:ctrl.tres
    fdat=cat(2,fdat,[fsol([1 5],g_fit)+[0 t/(2*ctrl.tres+1)]'
        fsol(14:16,g_fit)+fsol(17:19,g_fit)*t/(2*ctrl.tres+1)]);
end

% Define trajectories
ftraj=polytraj(fdat,ctrl.ord); % inherently weighted by camera occurence

% initiate cost
fsol(20,g_frm)=0; % to recompute and update by traj fit

% define invalid part trajectories
g_excl=false(1,size(fsol,2));

% loop cameras
for c=1:prop.res(4)
    
    % regularization for camera accuracy
    Dreg=prop.ecal(c)/sqrt(abs(det(Kmat{c}/Kmat{c}(end))));
    
    % camera set
    g_cam=fsol(4,:)==c;
    
    % fit camera data
    cdat=zeros(4,0);
    xdat=cvec2eshp(fsol(6:11, g_cam & g_fit ));
    for t=0:ctrl.tres%-ctrl.tres:ctrl.tres
        cdat=cat(2,cdat,[fsol([1 5], g_cam & g_fit )+[0 t/(2*ctrl.tres+1)]'
            xdat+fsol([12 13], g_cam & g_fit )*t/(2*ctrl.tres+1)]); 
    end
    
    %figure; plot(cdat(3,:),cdat(4,:),'.')
    
    % fit camera trajectories
    ctraj=polytraj(cdat,ctrl.ord);

    % expand disparity allowance
    cdat=contrans(fsol(6:11,g_cam & g_frm),ctrl.dspr+ctrl.dspl/(1+2*ctrl.tres),'resize'); % the higher the time resolution the higher the certainty of the position
    
    % camera uncertainty
    cdat=contrans(cdat,Dreg,'expand');

    % loop image displacement
    for t=0:ctrl.tres%-ctrl.tres:ctrl.tres
        
        % camera conic shift
        dshft=polyeval(ctraj,fsol([1 5], g_cam & g_frm )+[0 t/(2*ctrl.tres+1)]') - cvec2eshp(fsol(6:11,g_cam & g_frm));
        
        % camera conic along camera trajectory
        con=contrans(cdat,dshft,'displace'); %fsol(12:13,g_cam & g_frm)*t/(2*ctrl.tres+1) fsol(6:11,g_cam & g_frm);%
        
        % object track
        X=polyeval(ftraj,fsol([1 5], g_cam & g_frm )+[0 t/(2*ctrl.tres+1)]',0);%fsol(14:16, g_cam & g_frm )+fsol(17:19, g_cam & g_frm )*t/(2*ctrl.tres+1);%
        
        %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
        
        % camera reference frame
        X=Pmat{c}*inhc2homc(X) ;
        
        %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
        
        % segment dof
        g_excl( g_cam & g_frm )=g_excl( g_cam & g_frm ) | X(3,:)<ctrl.ldof(1) | X(3,:)>ctrl.ldof(2);
        
        % project coordinates
        x=homc2inhc( X ); % empty cell can comprimise
        
        %figure; plot(x(1,:),x(2,:),'.')
        
        % compute ellipse distance
        d=ellipsedist(con,x,'vector');
        
        %figure; plot(d,'.')
        
        %     % reprojection
        %     x=cvec2eshp(con);
        %     r=sqrt(sum((x(:,ai)-xp).^2,1));
        
        % update exclusion
        g_excl( g_cam & g_frm )=g_excl( g_cam & g_frm ) | d > 1;
        
        % write
        fsol(20, g_cam & g_frm )=fsol(20, g_cam & g_frm )+d/(2*ctrl.tres+1);
        
    end % t
    
    % position
    xdat=cvec2eshp(fsol(6:11, g_cam & g_frm ));
    
    %hold on; plot(xdat(1,:),xdat(2,:),'.')
    
    % evaluate
    ydat=polyeval(ctraj,fsol([1 5], g_cam & g_frm ),0);
    
    % evaluate
    vdat=polyeval(ctraj,fsol([1 5], g_cam & g_frm ),1);
    
    %hold on; plot(ydat(1,:),ydat(2,:),'.')
    
    % correct coordinates
    fsol(6:11, g_cam & g_frm )=contrans(fsol(6:11, g_cam & g_frm ),ydat-xdat,'displace');
    
    % correct vel
    fsol(12:13, g_cam & g_frm )=vdat;
    
end % c

% correct coordinates
fsol(14:16, g_frm )=polyeval(ftraj,fsol([1 5], g_frm ),0);

% correct vel
fsol(17:19, g_frm )=polyeval(ftraj,fsol([1 5], g_frm ),1);

%figure; plot(fsol(20,:),'.')

% remove wrong triangulation as whole & fsol(5,:)==n
g_excl=g_excl | ( fsol(5,:)>n & ismember(fsol(1,:),fsol(1,g_excl)) ); % fsol(5,:)>=n & ismember(fsol(1,:),fsol(1,g_excl)) ;

% remove exclusion from trajectory data
nsol=cat(2,nsol,fsol(:,g_excl));
fsol=fsol(:, ~g_excl ); 

%-- #) Minimum focality

% frame set
g_frm=fsol(5,:)==n ;%cpln(2,:)==c;

% track solution mapping
[ufl,~,fl]=unique(fsol(1, g_frm )); %ismember(fsol(5,:),uc) 
[ucl,~,cl]=unique(fsol(3, g_frm )); % ismember(fsol(5,:),uc)
F_sol=sparse(cl,fl,ones(size(fl)),length(ucl),length(ufl));

f_sum=sum(F_sol,1);

% exclude minimal focality
g_excl= g_frm & ismember(fsol(1,:),ufl(f_sum < 3));

% remove wrong triangulation as whole & fsol(5,:)==n
g_excl= fsol(5,:)>=n & ismember(fsol(1,:),fsol(1,g_excl)) ;

% remove exclusion from trajectory data
nsol=cat(2,nsol,fsol(:,g_excl));
fsol=fsol(:, ~g_excl ); 

%-- #) write

Fsol=cat(2,Fsol,fsol);

%-- #) End message

% keep sets complete in previous focality
disp(['Filtered ',num2str(length(unique(nsol(1,:)))),...
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

