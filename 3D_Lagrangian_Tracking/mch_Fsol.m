function Fsol = mch_Fsol(Cpln,Mset,Fsol,n)
%mch_Fsol Branch feasible solution over matched sets and filter trajectory
%   extrapolation/continuation.
%
%   Input:
%       Camera plane ellipse idenitifications
%       Matched sets between ellipse idenitifications
%       Feasible solution
%       Frame to process
%   
%   Output:
%       Fsol
%   
%   working principle
%       1) branch geometry by matched sets
%       2) filter temporal coherence camera plane
%       3) triangulate solution [allowed to deviate from matched sets triang]
%       4) filter spurious triangulation by reprojection
%       5) filter object trajectorie coherence and compute trajectory error
%
%   TODO
%       There are some obvious speed ups e.g. double operations inside
%       instead of single operations outside loop, which can be optimized
%       for runtime.

%% Get globals
global prop ctrl plotting Pmat Kmat

%% Initiate data

% maximum index
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

%-- #) cluster focality

% Select time frame
cpln=Cpln(:, Cpln(3,:)>=( n + 1 - ctrl.tfit ) & Cpln(3,:)<=n );%  - 1

% corresponding matched sets
mset=Mset(:, ~ismember(Mset(1,:),Mset(1, ~ismember(Mset(2,:),cpln(1, cpln(3,:)==n )) )) ); % complete in focality

% related feasible solutions
fsol=Fsol(:, ismember(Fsol(1,:),Fsol(1,ismember(Fsol(5,:),cpln(1, cpln(3,:)==n )) )) );

%-- #) Open for reassessment

% feasible solution (for reindixing when tracks split and merge)
Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));

%-- 2) define adjacencies

% unique indices feasible track indices
uc=unique( cpln(1,:) ); % mset(2,:)fsol(4,:)

% adjacency camera data
cpi=(1:size(cpln,2))';
[ucp,ai,cp]=unique(cpln(1,:));
iucp=find(ismember(uc,ucp));
cp=iucp(cp);
C_pln=sparse(cp,cpi,ones(size(cpi)),length(uc),length(cpi));

c_cam=cpln(2,ai);

%figure; spy(C_pln)

% Define existing matched set
[ums,~,ms]=unique(mset(1, : ));
[umsc,~,msc]=unique(mset(2, : ));
iumsc=find(ismember(uc,umsc));
msc=iumsc(msc);
M_set=sparse(msc,ms,ones(size(ms)),length(uc),length(ums));

%figure; spy(M_set)
%figure; plot(sum(M_set,1),'.')

% define feasible solution adjacency
[ufl,~,fl]=unique(fsol(1,:));
G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));

%figure; spy(G_sol)
%figure; plot(sum(G_sol & fsol(3,:)==n ,2),'.')

% track solution mapping
[ufl,~,fl]=unique(fsol(1, fsol(4,:)==n )); %ismember(fsol(5,:),uc) 
[ucl,~,cl]=unique(fsol(5, fsol(4,:)==n )); % ismember(fsol(5,:),uc)
iucl=find(ismember(uc,ucl));
cl=iucl(cl);
F_sol=sparse(cl,fl,ones(size(fl)),length(uc),length(ufl));

f_sum=sum(F_sol,1);

%figure; spy(F_sol)
%figure; plot(sum(F_sol,1),'.')

%-- #) Branch geometry and split solutions

% adjacency
A=M_set'*F_sol;%( F_sol(:,dum) ); % 

%figure; spy(A)

% find
[ai,aj,ak]=find(A);

% validity
val=ak==f_sum(aj)';
ai=ai(val);
aj=aj(val);
ak=ak(val);

% redefine A
A=sparse(ai,aj,ak,size(A,1),size(A,2));

% find mapping to copy data
[gi,gj]=find( G_sol(aj,:) );

%unique track
[~,gs,~]=unique(gi);
ift_brn=fsol(2,gj(gs));

%figure; plot(ift_brn,'.')

% define unconnected and new feasible solutions
g_unc=sum(A,1)==0;%true(1,size(A,2));%
g_new=sum(A,2)==0;%true(1,size(A,1));%

%-- #) Write unbranched track data

% unmatched solution
f_sol=fsol(:,ismember(fsol(1,:),ufl(g_unc)));

%figure; scatter3(f_sol(6,:),f_sol(7,:),f_sol(8,:),[],f_sol(9,:),'.'); axis equal; view(2); drawnow;

%-- #) Write splitting tracks

% define re-indexing branching tracks
ifl_brn=imaxfl+(1:length(aj));
imaxfl=imaxfl+length(aj);

dum= find(~ismember(fsol(5,:),cpln(1, cpln(3,:)==n ))) ;
% find mapping to copy data
[gi,gj]=find( G_sol(aj,~ismember(fsol(5,:),cpln(1, cpln(3,:)==n ))) );
gj=dum(gj);
%figure; spy(G_sol(aj,:))

% pairing temporal index corresponces
pf_sol=zeros( size(fsol,1) , length(gj) );
pf_sol(1,:)=ifl_brn(gi); %
for i=2:size(pf_sol,1)
    pf_sol(i,:)=fsol(i,gj); %
end

%figure; scatter3(pf_sol(5,:),pf_sol(6,:),pf_sol(7,:),[],pf_sol(8,:),'.'); axis equal; view(2); drawnow;

% branch
B= M_set(:,ai) ;

% find
[bi,bj]=find(B);

%figure; spy(B)

% initiate branche positions
pf_new=[ifl_brn(bj)
    ift_brn(bj)
    c_cam(bi)
    n*ones(1,nnz(bi))
    uc(bi)
    zeros(8,nnz(bi))];

%figure; scatter3(pf_new(6,:),pf_new(7,:),pf_new(8,:),[],pf_new(9,:),'.'); axis equal; view(2); drawnow;

%-- #) Newly seeding matched sets that are not included in branches

% enummeration new solutions
ifl_new=imaxfl+(1:nnz(g_new));
% imaxfl=imaxfl+nnz(g_new);

% enummeration new solutions
ift_new=imaxtl+(1:nnz(g_new));
% imaxtl=imaxtl+nnz(g_new);

% new
B= M_set(:,g_new)  ;

% find
[bi,bj]=find(B);

%figure; spy(B)

% initiate new feasible tracks
f_new=[ifl_new(bj)
    ift_new(bj)
    c_cam(bi)
    n*ones(1,nnz(bi))
    uc(bi)
    zeros(8,nnz(bi))];

%-- #) Append

% write new tracking data
fsol=cat(2,f_sol,pf_sol,pf_new,f_new);%

%figure; plot3(fsol(6,:),fsol(7,:),fsol(8,:),'.'); axis equal; view(2); drawnow

%-- redefine adjacencies

% % track solution mapping
% [ufl,~,fl]=unique(fsol(1, : ));
% G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
% 
% %figure; spy(G_sol)

% track solution mapping
[ufl,~,fl]=unique(fsol(1, ismember(fsol(5,:),uc) ));
[ucl,~,cl]=unique(fsol(5, ismember(fsol(5,:),uc) )); %ismember(ftrk(3,:),cpln(1,:))
iucl=find(ismember(uc,ucl));
cl=iucl(cl);
F_sol=sparse(cl,fl,ones(size(fl)),length(uc),length(ufl));

%figure; spy(F_sol)
%figure; plot(sum(F_sol,1),'.')

%-- #) Filter connectivity trajectories

% define invalid part trajectories
g_excl=false(1,size(fsol,2));

% loop cameras
for c=1:prop.res(4)
    
    % camera set
    g_cam=cpln(2,:)==c;
    
    % data
    cdat=cpln(:,g_cam);
    
    % adjacencie camera and solution 
    C_trk=double(C_pln(:,g_cam)'*F_sol>0.5);%g_cam' & 
    
    %figure; spy(C_trk); drawnow
    %figure; plot(sum(C_trk,2),'.')
    
    % select frames in data
    N_trk=spdiags( cdat(3,:)' , 0 , size(cdat,2) , size(cdat,2) )*C_trk; % only exploration
    [n_max,ind]=max(N_trk,[],1); % 
    
    ind=cdat(1,ind);
    
    % find
    dum=find( n_max==n);
    [ni,nj]=find(N_trk(:, n_max==n)); %  ,nk
    nj=dum(nj);
    
    ind=ind(nj);%cdat(1,ni);
    
    % define trajectory data
    ctrk=zeros( 2 + 6 + 2 , length(ni) ); % [t ind n con vel]
    ctrk(1,:)=ufl(nj);
    for i=2:size(ctrk,1)
        ctrk(i,:)=cdat(i+1,ni);
    end
    
    % fit data
    cdat=zeros(4,0);
    xdat=cvec2eshp(ctrk(3:8,:));
    for t=-ctrl.tres:ctrl.tres
        cdat=cat(2,cdat,[ctrk(1:2,:)+[0 t/(2*ctrl.tres+1)]'
            xdat+ctrk(9:10,:)*t/(2*ctrl.tres+1)]); 
    end
    
    % fit trajectories
    ctraj=polytraj(cdat,ctrl.ord);
    
    % loop between frame
    for t=-ctrl.tres:ctrl.tres
        
        % resize to displacement
        cdat=contrans(ctrk(3:8,:),ctrl.dspl+ctrl.dspr,'resize');
        
        % displace conic
        cdat=contrans(cdat, ctrk(9:10,:)*t/(2*ctrl.tres+1) ,'displace');
        
        % evaluate trajectory
        ydat=polyeval(ctraj,ctrk(1:2, : )+[0 t/(2*ctrl.tres+1)]',0); % null
        
        %     xdat=cvec2eshp(ctrk(3:8, : ));
        %     figure; plot(ydat(1,:),ydat(2,:),'.')
        %     hold on; plot(xdat(1,:),xdat(2,:),'.'); drawnow
        
        % compute vector ellipse distance
        d=ellipsedist(cdat,ydat,'vector'); % coherence
        
        %figure; plot(d,'.'); drawnow
        
        % segment disparity
        g_rem=d > 1;
        
        % update removal
        g_excl=g_excl | ismember(fsol([1 5],:)',[ctrk(1,g_rem)' ind(g_rem)'],'rows')'; %
        
    end % t
    
end % c

% remove exclusion from trajectory data
fsol=fsol(:, ~g_excl );

%-- redefine adjacencies

% track solution mapping
[ufl,~,fl]=unique(fsol(1, : ));
G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));

%figure; spy(G_sol)

% track solution mapping
[ufll,~,fl]=unique(fsol(1, ismember(fsol(5,:),uc) ));
iufl=find(ismember(ufl,ufll));
fl=iufl(fl);
[ucl,~,cl]=unique(fsol(5, ismember(fsol(5,:),uc) )); %ismember(ftrk(3,:),cpln(1,:))
iucl=find(ismember(uc,ucl));
cl=iucl(cl);
F_sol=sparse(cl,fl,ones(size(fl)),length(uc),length(ufl));

%figure; spy(F_sol)
%figure; plot(sum(F_sol,1),'.')

%-- #) Triangulation

% Assemble triangulation data
xdat=cell(1,prop.res(4));
for c=1:prop.res(4)
    
    g_cam=cpln(2,:)==c & cpln(3,:)==n;
    
    % Get midpoints for triangulation
    x=cvec2eshp(cpln(4:9,g_cam));
    
    % adjacencie camera data and solution
    A=C_pln(:,g_cam)'*F_sol;
    
    % find
    [ai,aj]=find(A);
    
    xdat{c}=nan*zeros(3,size(F_sol,2)); %0
    
    xdat{c}(:,aj)=inhc2homc(x(:,ai));
    
end

% Triangulate 
X=objtriang(xdat,Pmat);%[,~,r] 

% find indices
dum=find(fsol(4,:)==n);
[gi,gj]=find(G_sol(:,fsol(4,:)==n));

% write position data
fsol(6:8,dum(gj))=X(:,gi);

%-- #) Displacement

% Assemble triangulation data
xdat=cell(1,prop.res(4));
for c=1:prop.res(4)
    
    g_cam=cpln(2,:)==c & cpln(3,:)==n;
    
    % Get midpoints for triangulation
    x=cvec2eshp(cpln(4:9,g_cam))+cpln(10:11,g_cam);
    
    % adjacencie camera data and solution
    A=C_pln(:,g_cam)'*F_sol;
    
    % find
    [ai,aj]=find(A);
    
    xdat{c}=nan*zeros(3,size(F_sol,2)); %0
    
    xdat{c}(:,aj)=inhc2homc(x(:,ai));
    
end

% Triangulate 
Xacc=objtriang(xdat,Pmat); % [,~,r] 

% figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2); colorbar; drawnow

% velocity
V=Xacc-X;

% find indices
dum=find(fsol(4,:)==n);
[gi,gj]=find(G_sol(:,fsol(4,:)==n));

% write position data
fsol(9:11,dum(gj))=V(:,gi);

% figure; plot3(fsol(6,:),fsol(7,:),fsol(8,:),'.'); axis equal; view(2); drawnow

% exclude nan values
fsol = fsol(:, ~isnan(sum(fsol,1))); 

%-- Exclusion and redefine adjacencies

% track solution mapping
[ufl,~,fl]=unique(fsol(1, : ));
G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));

%figure; spy(G_sol)

% track solution mapping
[ufll,~,fl]=unique(fsol(1, ismember(fsol(5,:),uc) ));
iufl=find(ismember(ufl,ufll));
fl=iufl(fl);
[ucl,~,cl]=unique(fsol(5, ismember(fsol(5,:),uc) )); %ismember(ftrk(3,:),cpln(1,:))
iucl=find(ismember(uc,ucl));
cl=iucl(cl);
F_sol=sparse(cl,fl,ones(size(fl)),length(uc),length(ufl));

%figure; spy(F_sol)
%figure; plot(sum(F_sol,1),'.')

%-- #) Compute reprojection error and update statistics

% change body size by disparity allowance
cpln(4:9, : )=contrans(cpln(4:9, : ),ctrl.dspl+ctrl.dspr,'resize');

% Regularize conics at camera precision for evaluation
for c=1:prop.res(4)
    
    g_cam=cpln(2,:)==c;
    
    Dreg=prop.ecal(c)/sqrt(abs(det(Kmat{c}/Kmat{c}(end))));
    
    cpln(4:9,g_cam)=contrans(cpln(4:9,g_cam),Dreg,'expand');
    
end

% define invalid part trajectories
g_excl=false(1,size(fsol,2));

% loop cameras
for c=1:prop.res(4)
    
    % loop image displacement
    for t=-ctrl.tres:ctrl.tres
        
        % camera set
        g_cam=cpln(2,:)==c & cpln(3,:)==n;
        
        % solution set
        g_sol=fsol(3,:)==c & fsol(4,:)==n;
        dum=find(g_sol);
        
        % camera conic data
        con=contrans(cpln(4:9,g_cam),t*cpln(10:11,g_cam),'displace');
        
        % object triangulation
        X=fsol(6:8,g_sol)+t*fsol(9:11,g_sol);
        
        %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
        
        % find adjacency
        [ai,aj]=find(C_pln(:,g_cam)'*F_sol*G_sol(:,g_sol) );
        
        % camera reference frame
        X=Pmat{c}*inhc2homc(X) ;
        
        %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
        
        % segment dof
        g_excl(dum(aj))=g_excl(dum(aj)) | X(3,aj)<ctrl.ldof(1) | X(3,aj)>ctrl.ldof(2);
        
        % project coordinates
        x=homc2inhc( X ); % empty cell can comprimise
        
        %figure; plot(xp(1,:),xp(2,:),'.')
        
        % compute ellipse distance
        d=ellipsedist(con(:,ai),x(:,aj),'vector');
        
        %figure; plot(d,'.')
        
        %     % reprojection
        %     x=cvec2eshp(con);
        %     r=sqrt(sum((x(:,ai)-xp).^2,1));
        
        % update exclusion
        g_excl(dum(aj))=g_excl(dum(aj)) | d > 1;
        
        % write
        fsol(12,dum(aj))=fsol(12,dum(aj))+d/(2*ctrl.tres+1);
        
    end % t
    
end % c

%figure; plot(fsol(12,:),'.')

% remove exclusion from trajectory data
fsol=fsol(:, ~(ismember(fsol(1,:),fsol(1,g_excl)) & fsol(4,:)==n) );

%-- redefine adjacencies

% track solution mapping
[ufl,~,fl]=unique(fsol(1, : ));
G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));

%figure; spy(G_sol)

% track solution mapping
[ufll,~,fl]=unique(fsol(1, ismember(fsol(5,:),uc) ));
iufl=find(ismember(ufl,ufll));
fl=iufl(fl);
[ucl,~,cl]=unique(fsol(5, ismember(fsol(5,:),uc) )); %ismember(ftrk(3,:),cpln(1,:))
iucl=find(ismember(uc,ucl));
cl=iucl(cl);
F_sol=sparse(cl,fl,ones(size(fl)),length(uc),length(ufl));

%figure; spy(F_sol)
%figure; plot(sum(F_sol,1),'.')

%-- #) Filter object trajectories

% select frames in data
dum=find( ismember(fsol(5,:),cpln(1,:)) );
N_sol=G_sol(:, ismember(fsol(5,:),cpln(1,:)) ) ;

%figure; spy(N_sol)

N_sol=N_sol*spdiags(fsol(4,ismember(fsol(5,:),cpln(1,:)))',0,length(dum),length(dum));
n_max=max(N_sol,[],2);

% find
[~,nj]=find(N_sol(n_max==n,:) );
nj=dum(nj);

% data in track
fobj=zeros( 2+3+3, length(nj) ); % [t n ind pos vel]
fobj(1,:)=fsol(1,nj);
fobj(2,:)=fsol(4,nj);
for i=3:size(fobj,1)
    fobj(i,:)=fsol(i+3,nj);
end

%figure; plot3(fobj(3,:),fobj(4,:),fobj(5,:),'.'); axis equal; view(2); drawnow

% fit data
fdat=zeros(5,0);
for t=-ctrl.tres:ctrl.tres
    fdat=cat(2,fdat,[fobj(1:2,:)+[0 t/(2*ctrl.tres+1)]'
        fobj(3:5,:)+fobj(6:8,:)*t/(2*ctrl.tres+1)]);
end

% Define trajectories
ftraj=polytraj(fdat,ctrl.ord); % inherently weighted by camera occurence

% unique track frames
[~,ia,~]=unique(fobj([1 2],:)','rows'); % removing the weight above
fobj=fobj(:,ia);

% track solution mapping
[ufo,~,fo]=unique(fobj(1, : ));
iufo=find(ismember(ufl,ufo));
fo=iufo(fo);
G_obj=sparse(fo,1:length(fo),ones(size(fo)),length(ufl),length(fo));

%figure; spy(G_obj)

% % initiate exclusion set
g_excl=false(1,size(fobj,2));

% loop cameras
for c=1:prop.res(4)
    
    % total disparity
    D=zeros(1,size(fobj,2)); % per view
    
    % loop trajectory frames
    for np= max(ctrl.tproc(1),( n + 1 - ctrl.tfit )) : n %
        
        % camera set
        g_cam=cpln(2,:)==c & cpln(3,:)==np;
        
        % object set
        g_obj=fobj(2,:)==np;
        dum=find(g_obj);
        
        % camera conic data
        cdat=cpln(:,g_cam);
        
        % trajectory data
        fdat=fobj(:,g_obj);
        
        % find adjacency
        [ai,aj]=find( C_pln(:,g_cam)'*F_sol*G_obj(:,g_obj) ); %
        
        % loop frame interpolant time stencil
        for t=-ctrl.tres:ctrl.tres
            
            % change body size to displacement
            con=contrans(cdat(4:9,:),ctrl.dspl,'resize');
            
            % change body size by [halve] step in velocity
            con=contrans(con, cdat(10:11,:)*t/(2*ctrl.tres+1) ,'displace'); % undirected motion + uncertainty 1/2+1/2
            
            % evaluate track position
            X=polyeval(ftraj,fdat(1:2,:)+[0 t/(2*ctrl.tres+1)]' ,0); 
            
            %figure; plot3(X(3,:),X(4,:),X(5,:),'.'); axis equal; view(2); drawnow
            
%             % evaluate track velocity
%             V=polyeval(ftraj,fdat(1:2,:)+[0 t/(2*ctrl.tres+1)]' ,1); 
            
            %figure; plot3(X(3,:),X(4,:),X(5,:),'.'); axis equal; view(2); drawnow
            
            % camera refernce frame
            X=Pmat{c}*inhc2homc(X) ;
            
            %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
            
            % segment dof
            g_excl(dum(aj))=g_excl(dum(aj)) | X(3,aj)<ctrl.ldof(1) | X(3,aj)>ctrl.ldof(2);
            
            % project track
            x=homc2inhc(X);
            
            % figure; plot(x(1,:),x(2,:),'.'); drawnow
            
            % compute disparity on track
            d=ellipsedist(con(:,ai),x(:,aj),'vector');
            
            %figure; plot(d,'.'); drawnow;
            
            % segment disparity
            g_rem=d>1;
            
            % update exclusion
            g_excl(dum(aj))=g_excl(dum(aj)) | g_rem;
            
            %         % reprojection on dewarped image (allows better crossing ellipses, not preferring big ellipses)
            %         x=cvec2eshp(con);
            %         R0=sqrt(sum((x(:,ai)-xp0).^2,1)); % more convex, segmentation already done
            
%             % Swept reprojection error in object space
%             d=d.*sqrt(sum(V(:,aj).^2,1));
            
            % write track disparity
            D(dum(aj))=D(dum(aj))+d/(2*ctrl.tres+1);
            
        end % t
        
    end % np
    
    %figure; plot(D,'.'); drawnow
    %figure; histogram(D,50); drawnow
    
    % mean track error
    D=(D*G_obj')./sum(G_obj,2)'; 
    
    %figure; plot(D,'.'); drawnow
    
    % write track error
    dum=find(fsol(3,:)==c & fsol(4,:)==n);
    [gi,gj]=find(G_sol(:,fsol(3,:)==c & fsol(4,:)==n));
    fsol(13,dum(gj))=D(gi);
    
end % c 

%figure; plot(fsol(13,:),'.')

% exclude spurious solutions
fsol=fsol(:, ~(ismember(fsol(1,:),fobj(1,g_excl)) & fsol(4,:)==n));

%figure; plot3(fsol(5,:),fsol(6,:),fsol(7,:),'.'); axis equal; view(2); drawnow

%-- #) write

Fsol=cat(2,Fsol,fsol);

%-- #) End message

% keep sets complete in previous focality
disp(['Matched ',num2str(length(unique(fsol(1,fsol(4,:)==n)))),...
    ' new out of ',num2str(length(unique(fsol(1,:)))),' sets'])

%% plot results
if strcmp(plotting,'on')
    
    figure(prop.res(4)+1)
    cla
    
    plot3(X(1,:),X(2,:),X(3,:),'o');
    hold on
    scatter3(fsol(6,fsol(4,:)==n),fsol(7,fsol(4,:)==n),fsol(8,fsol(4,:)==n),...
        [],fsol(9,fsol(4,:)==n),'.')
    hold off
    axis equal
    view(2)
    colorbar

    title(['frame ',num2str(n)])
    
    xlim([-3 3]); ylim([11 17])
    
    drawnow
end

end

%     mu1=(D*G_obj')./sum(G_obj,2)'; % 1-norm mean
%     mu2=sqrt((D.^2*G_obj')./sum(G_obj,2)'); % 2-norm fit
%     sig=sqrt( abs( mu2.^2 - mu1.^2 ) ); % variance (avoid imagionary val comp acc)
%     D=mu1+sqrt(3)*sig; % error to minimize based on a uniform distrobution
    