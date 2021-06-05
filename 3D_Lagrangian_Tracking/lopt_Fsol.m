function Fsol=lopt_Fsol(Fsol,n)
%opt_Fsol Locally optimize feasible solutions.
%
%   Input:
%       Fsol    Feasible solutions
%       n       Current frame.
%   
%   Output:
%       Fsol    Feasible solutions.
%   
%   working principle:
%       Iterative greedy solution by cost function sortation. Compute
%       optimal solution in new frame appendages, then compute overlap
%       between solution and penalize overlapping coordinates.
%   

%% Get globals
global plotting ctrl % prop

%% Branch and optimize feasible solutions

%-- 0) Begin message

disp(['Optimize feasible solutions frame ',num2str(n)])

%-- #) Initiate data

% get feasible solution
fsol=Fsol( : , ismember( Fsol(1,:) , Fsol(1,Fsol(5,:)>=(n-ctrl.prd) & Fsol(5,:)<=n ) ) );% Fsol(5,:)==n

%--#) Open for reaccessment

% feasible solutions
Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));

%-- #) Define adjacencies

% define feasible solution adjacency
[ufl,~,fl]=unique(fsol(1, : ));
G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));

% track sets
g_min= fsol(5,:)< n;
g_plus=fsol(5,:)>=n; 

%figure; spy(G_sol)

% feasible track index mapping
[ufl,~,fl]=unique(fsol(1, : )); 
[utl,~,tl]=unique(fsol(2, : )); 
T_sol=double(sparse(tl,fl,ones(size(tl)),length(utl),length(ufl))>0.5); % else >0.5 !

%figure; spy(T_sol)
%figure; plot(sum(T_sol>0.5,1),'.')
%figure; plot(sum(T_sol>0.5,2),'.') 

% new camera coordinates mapping
[ufp,~,fp]=unique(fsol(1, g_plus )); 
iufp=find(ismember(ufl,ufp));
fp=iufp(fp);
[ucp,~,cp]=unique(fsol(3, g_plus )); 
P_sol=sparse(cp,fp,ones(size(cp)),length(ucp),length(ufl));%

%figure; spy(P_sol)

%-- #) Initiate optimization

% % track start
% [~,dum]=max(G_sol*spdiags(max(fsol(5,:))+1-fsol(5,:)',0,size(fsol,2),size(fsol,2)),[],2);
% n_min=fsol(5,dum);
% 
% %figure; plot(n_min,'.'); drawnow
% 
% % track end
% n_max=max(G_sol*spdiags(fsol(5,:)',0,size(fsol,2),size(fsol,2)),[],2)';
% 
% %figure; plot(n_max,'.'); drawnow

% camera weight track
wfm=sum( G_sol(:, g_min  ) ,2)'; % previous tracking
wfp=sum( G_sol(:, g_plus ) ,2)'; % new segments

%figure; plot(wfm,'.')
%figure; plot(wfp,'.')
%figure; plot(wfm+wfp,'.')

% average new disparity/dissimilarity
dfm=((fsol(20, g_min  ))*G_sol(:, g_min  )')...
    ./sum(G_sol(:, g_min  ),2)'; % previous tracking
dfp=((fsol(20, g_plus ))*G_sol(:, g_plus )')...
    ./sum(G_sol(:, g_plus ),2)'; % new segments 

%figure; plot(dfm,'.')
%figure; plot(dfp,'.')
%figure; plot(nanmean([dfm ; dfp],1),'.')

%-- 4) Iterative greedy optimization using cost function sortation

% initiate open set
g_ope= true(1,length(ufl)); 

% initiate solution set
g_sol=false(1,length(ufl)); 

% iterate cost function
while nnz(g_ope)~=0
    
    % open solution
    ope_sol=find(g_ope);
    
    % open solution overlap penalty total overlap
    ope_pfp=sum(P_sol( sum( P_sol(:,g_sol) ,2)>0.5 ,ope_sol),1);
    
    %figure; plot(ope_pf,'.'); drawnow
    
    % sort cost function -n_min(ope_sol)' n_max(ope_sol)' wfp(ope_sol)' 
    [~,ind]=sortrows([wfm(ope_sol)' -dfm(ope_sol)' ( wfp(ope_sol)'-ope_pfp' ) -dfp(ope_sol)' fliplr(1:length(ope_sol))'],'descend');
    
    % sort open solution to cost
    ope_sol=ope_sol(ind);
    
    %figure; plot(nanmean([dfm(ope_sol) ; dfp(ope_sol)],1),'.'); drawnow
    
    % solve optimal cost by maximize the sortation (T_sol can be important here)
    [val,subsol]=max( double( T_sol(:,ope_sol) >0.5 )*spdiags(fliplr(1:length(ope_sol))',0,length(ope_sol),length(ope_sol)) ,[],2); % first unique el
    
    %figure; plot(sum(T_sol(:,ope_sol),2),'.') 

    % current optimal solution
    opt_sol=ope_sol(unique(subsol(val>0)));
    
    % adjancency new camera coordinates in new solution
    A=tril(P_sol(:,opt_sol)'*P_sol(:,opt_sol),-1); % symmetric, keep lower part
    
    %figure; spy(A)
    
    % current free solution (found first in overlap)
    free_sol=opt_sol( sum(A,2)==0  );
    
    % write free solution
    g_sol(free_sol)=1; 
    
    % update open set to picked current solution
    g_ope(sum(T_sol(sum(T_sol(:,free_sol),2)>0.5,:),1)>0.5)=0;
    
    % find open solution
    ope_sol=find(g_ope);
    
    % adjancency new camera coordinates to free solution
    A=P_sol(:,ope_sol)'*P_sol(:,free_sol); % asymmetric
    
    %figure; spy(A)
    
    % find nonzero adjacency
    [ai,~,ak]=find(A);
    
    % full overlapping solution
    ovlp_sol=ope_sol(unique(ai(ak==wfp(ope_sol(ai))')));
    
%     % update wf here
%     wfp(ovlp_sol)=0; % penalize wfp only in complete overlap
    
    % remove overlapping solution when being only new appendage
    rem_sol=ovlp_sol(wfm(ovlp_sol)==0);
    
    % update opensolution
    g_ope(rem_sol)=0;
    
end

% keep valid solutions
fsol=fsol(:,ismember(fsol(1,:),ufl(g_sol)));

%-- #) Write data

% feasible solutions
Fsol=cat(2,Fsol,fsol);

%-- #) End message

disp(['Optimized ',num2str(nnz(g_sol)),' out of ',num2str(length(ufl)),' feasible solutions'])

%% plot results

%-- #) Make figure

if strcmp(plotting,'on')
   
    % frame
    N=Fsol(5,:)==n;
    
    % figure cost function
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(20,N),'.')
    axis tight
    axis equal
    view(2)
    title(['Optimal solution, cost value frame ',num2str(n)])
    colorbar
    drawnow
    
    % figure indexing
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(1,N),'.')
    colormap lines
    axis tight
    axis equal
    view(2)
    title(['Optimal solution, indexing frame ',num2str(n)])
    drawnow
    
end

end

