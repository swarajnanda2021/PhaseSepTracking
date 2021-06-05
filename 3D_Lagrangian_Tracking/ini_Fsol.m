function Fsol = ini_Fsol(Cpln,Mset,Fsol,n)
%ini_Fsol Initiate and seed feasible solution over matched sets
%   extrapolation/continuation.
%
%   Input:
%       Cpln    Camera plane ellipse idenitifications
%       Mset    Matched sets between ellipse idenitifications
%       Fsol    Feasible solution
%       n       Current frame.
%   
%   Output:
%       Fsol    Feasible solution with new seeding.
%   
%  Working principle:
%       Seed new feasible solution from unused matched sets.
%   

%% Get globals
global plotting % ctrl Pmat prop 

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

%-- #) Begin message

% keep sets complete in previous focality
disp(['Seed new feasible solutions frame ',num2str(n)])

%-- #) cluster focality

% Current feasible solutions
fsol=Fsol(:, ismember( Fsol(1,:), Fsol(1,...
    Fsol(5,:) >= n-1 & Fsol(5,:) <= n+1 ) ) ); % commented in for Junaid 

% Select time frame
cpln=Cpln(:, Cpln(3,:) >= n-1 & Cpln(3,:) <= n+1 ); % commented in for Junaid 

% corresponding matched sets
mset=Mset(:, ~ismember(Mset(1,:),Mset(1, ~ismember(Mset(4,:),cpln(1, : )))));

%-- #) Open for reassessment

% feasible solution (for reindixing when tracks split and merge)
Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));

%-- 2) define adjacencies

% adjacency camera data
[uc,~,cp]=unique(cpln(1,:)); %ai
C_pln=sparse(cp,1:length(cp),ones(size(cp)),length(uc),length(cp));

%figure; spy(G_pln)

% Define existing matched set
[ums,mu,ms]=unique(mset(1, : ));
[umsc,~,msc]=unique(mset(4, : ));
iumsc=find(ismember(uc,umsc));
msc=iumsc(msc);
M_set=sparse(msc,ms,ones(size(ms)),length(uc),length(ums));

m_sum=sum(M_set,1);

%figure; spy(M_set)
%figure; plot(sum(M_set,1),'.')

% track solution mapping
[ufl,~,fl]=unique(fsol(1, ismember(fsol(3,:),uc)  )); %ismember(fsol(5,:),uc) 
[ucl,~,cl]=unique(fsol(3, ismember(fsol(3,:),uc)  )); % ismember(fsol(5,:),uc)
iucl=find(ismember(uc,ucl));
cl=iucl(cl);
F_sol=sparse(cl,fl,ones(size(fl)),length(uc),length(ufl));

%figure; spy(F_sol)
%figure; plot(f_sum,'.')

%-- #) Branch geometry and split solutions

% adjacency
A=M_set'*F_sol; % double(sum(F_sol,2)>0.5);%

%figure; spy(A)

% find
[ai,aj,ak]=find(A);

% validity
val=ak==m_sum(ai)';
ai=ai(val);
aj=aj(val);
ak=ak(val);

% redefine A
A=sparse(ai,aj,ak,size(A,1),size(A,2));

% define new seeding
g_new=sum(A,2)==0; %

%-- #) Newly seeding matched sets that are not included in branches

% enummeration new solutions
ifl_new=imaxfl+(1:nnz(g_new));
% imaxfl=imaxfl+nnz(g_new);

% enummeration new solutions
ift_new=imaxtl+(1:nnz(g_new));
% imaxtl=imaxtl+nnz(g_new);

% find
dum=find(g_new);
[bi,bj]=find( C_pln'*M_set(:,g_new) );
bi=reshape(bi,1,[]);
bj=reshape(bj,1,[]); % because when appending 1 -> [1 1]' but should be [1 1]

% initiate new feasible tracks
f_new=[ifl_new(bj)
    ift_new(bj)
    cpln(1:3,bi)
    cpln(4:9,bi)
    cpln(10:11,bi)
    mset(5:7,mu(dum(bj)))
    mset(8:10,mu(dum(bj)))
    mset(11,mu(dum(bj)))];

% %-- #) Redefine adjacencies and sets
% 
% % track solution mapping
% [ufl,~,fl]=unique(f_new(1, : ));
% G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));
% 
% %figure; spy(G_sol)
% 
% %--#) Shift ellipses from object space
% 
% % loop cameras
% cdat=cell(1,prop.res(4));
% % Xdat=nan(3,length(ufl));
% for c=1:prop.res(4)
%     
%     % camera set
%     g_cam=f_new(4,:)==c ;
%     
%     % camera conics
%     con=f_new(6:11, g_cam );
%     
%     %Xe=reshape(cvec2pnts(con),2,[]);
%     %figure; plot(Xe(1,:),Xe(2,:),'.r')
%     
%     % postion
%     x=cvec2eshp(con);
%     
%     % object position
%     X=f_new(14:16, g_cam  );%
% %     V=f_new(17:19, g_cam);%
%     
%     %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
%     
%     y=homc2inhc( Pmat{c}*inhc2homc(X) );
%     
%     con=contrans(con,y-x,'displace');
%     
% %     % project coordinates
% %     v=homc2inhc( Pmat{c}*inhc2homc(X+V) )-y; % empty cell can comprimise
% % 
% %     f_new(6:11, g_cam )=con;
% %     f_new(12:13, g_cam)=v;
%     
%     %Xe=reshape(cvec2pnts(con),2,[]);
%     %figure; plot(Xe(1,:),Xe(2,:),'.r')
%     
%     % adjacencie camera data and solution
%     A=G_sol(:,g_cam);
%     
%     % find
%     [ai,aj]=find(A);
%     
%     % conic data
%     cdat{c}=nan(6,length(ufl));
%     cdat{c}(:,ai)=con(:,aj);
% %     Xdat(:,ai)=X(:,aj);%quick fix overwrite
%     
% end
% 
% % constrained triangulation
% Q=quadrec(cdat,Pmat);%,Xdat);
% 
% % positiveness
% [~,g_val]=qvec2pnts(Q);
% Q=Q(:,g_val);%(:,~g_val)=nan;
%     
% %Xq=reshape(qvec2pnts(Q),3,[]);
% %figure; plot3(Xq(1,:),Xq(2,:),Xq(3,:),'r.'); axis tight; axis equal; view(2); drawnow
% %xlim([-3 3]); ylim([11 17]); zlim([-5 1])
% 
% % project quadric
% for c=1:prop.res(4)
%     
%     % camera set
%     g_cam=f_new(4,:)==c; 
%     
%     % Get midpoints for triangulation
%     conim=f_new(6:11,g_cam );
%     
%     %Xe=reshape(cvec2pnts(conim),2,[]);
%     %figure; plot(Xe(1,:),Xe(2,:),'.r')
%     
%     % postion
%     x=cvec2eshp(conim);
%     
% %     % object position
% %     X=f_new(14:16, g_cam  );%
% %     
% %     %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
% %     
% %     y=homc2inhc( Pmat{c}*inhc2homc(X) );
% %     
% %     conim=contrans(conim,y-x,'displace');
%     
%     % adjacencie camera data and solution
%     A=G_sol(g_val,g_cam);
%     
%     % find
%     [ai,aj]=find(A);
%     
%     con=conim;
%     for loop=1:length(aj)
%         % conic vector
%         con(:,aj(loop))=cmat2cvec(inv(Pmat{c}*inv(qvec2qmat(Q(:,ai(loop))))*Pmat{c}'));%*(0.5-rand(1)); % randly scale
%     end
%     
%     %Xe=reshape(cvec2pnts(con),2,[]);
%     %figure; plot(Xe(1,:),Xe(2,:),'.r')
%     
%     % postion
%     y=cvec2eshp(con);
%     
%     con=contrans(con,x-y,'displace');
%     
%     %Xe=reshape(cvec2pnts(con),2,[]);
%     %figure; plot(Xe(1,:),Xe(2,:),'.r')
%     
%     [~,~,~,peak1]=cvec2lset(con,cvec2eshp(con));
%     [~,~,~,peak2]=cvec2lset(conim,cvec2eshp(conim));
%     
%     cdat=((peak2./peak1).*con + conim)/2; % nanmean(cat(3,contrans(con,[],'normalize'),contrans(conim,[],'normalize')),3);
%     
%     %Xe=reshape(cvec2pnts(cdat),2,[]);
%     %figure; plot(Xe(1,:),Xe(2,:),'.r')
%     
%     f_new(6:11, g_cam )=cdat;
%     
% end

%-- #) Append

% write new tracking data
fsol=cat(2,fsol,f_new);

%figure; plot3(fsol(6,:),fsol(7,:),fsol(8,:),'.'); axis equal; view(2); drawnow

%-- #) write

Fsol=cat(2,Fsol,fsol);

%-- #) End message

% keep sets complete in previous focality
disp(['Seeded ',num2str(length(unique(f_new(1,:)))),...
    ' new out of ',num2str(length(unique(fsol(1,:))))])

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
    title(['New seeded solution, cost value frame ',num2str(n)])
    colorbar
    drawnow
    
    % figure indexing
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(1,N),'.')
    colormap lines
    axis tight
    axis equal
    view(2)
    title(['New seeded solution, indexing frame ',num2str(n)])
    drawnow
    
end

end

