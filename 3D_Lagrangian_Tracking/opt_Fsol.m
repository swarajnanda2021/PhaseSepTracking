function Fsol=opt_Fsol(Cpln,Fsol,n)
%opt_Fsol Optimize feasible solution by unique selecting best trajectory
%   with least camera overlap.
%
%   Input:
%       Camera plane ellipse idenitifications
%       Feasible solution
%       Frame to process
%   
%   Output:
%       Fsol
%
%   working principle
%       Optimize by iteratively (D&C) selecting best unique trajectories by 
%       computing a penalty for overlap and resorting the cost function
%       (greedy solution). Somewhere in between Divide and Conquer and
%       Greedy solution, is that the best language to classify this?

%% Get globals
global plotting ctrl prop

%% Branch and optimize feasible solutions

%-- 0) Begin message

disp(['Optimize feasible solutions frame ',num2str(n)])

%-- #) Initiate data

% get track segments in frames floor(ctrl.tfit*1/2)
cpln=Cpln( : , Cpln(3,:)>=( n - ctrl.prd ) & Cpln(3,:)<=n );

% get feasible solution
fsol=Fsol(:,ismember(Fsol(1,:),...
    Fsol(1,ismember(Fsol(5,:),cpln(1,:))))); 
    
%--#) Open for reaccessment

% feasible solutions
Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));

%-- #) Define adjacencies

% define feasible solution adjacency
[ufl,~,fl]=unique(fsol(1, : ));%ismember(fsol(1,:),ufl)
G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));

%figure; spy(G_sol)

% track solution mapping
[ufl,~,fl]=unique(fsol(1, : )); % ismember(fsol(5,:),uc)
[utl,~,tl]=unique(fsol(2, : )); %ismember(fsol(5,:),uc)double(~isnan(sum(fsol,1)))
T_sol=sparse(tl,fl,ones(size(tl)),length(utl),length(ufl));%

%figure; spy(T_sol)
%figure; plot(sum(T_sol>0.5,1),'.')
%figure; plot(sum(T_sol>0.5,2),'.') % >0.5 !

% track solution mapping
[ufl,~,fl]=unique(fsol(1, : )); % 
[usl,~,sl]=unique(fsol(5, : )); %double(~isnan(sum(fsol,1)))
S_sol=sparse(sl,fl,ones(size(sl)),length(usl),length(ufl));%

%figure; spy(S_sol)
%figure; plot(sum(S_sol>0.5,1),'.')
%figure; plot(sum(S_sol>0.5,2),'.')

%-- #) Initiate optimization

% track start
[~,dum]=max(G_sol*spdiags(n+1-fsol(4,:)',0,size(fsol,2),size(fsol,2)),[],2);
n_min=fsol(4,dum);

% figure; plot(n_min,'.'); drawnow

% track end
n_max=max(G_sol*spdiags(fsol(4,:)',0,size(fsol,2),size(fsol,2)),[],2)';

% figure; plot(n_max,'.'); drawnow

% track set
g_trk=fsol(4,:)~=n;%~ismember(fsol(5,:),uc); % not within new exploration

% camera weight track
wf=sum( G_sol(:,g_trk) ,2)'; %(:,g_trk)C_pln'*C_cam*average

% figure; plot(wf,'.'); drawnow

% average track error (sortation first, non-unique)
ef=(fsol(13,g_trk)*G_sol(:,g_trk)')...
    ./sum(G_sol(:,g_trk),2)'; % mean

% figure; plot(ef,'.'); drawnow

% camera cover new
cf=sum( G_sol(:,~g_trk) ,2)'; %(:,g_trk)C_pln'*C_cam*average

% figure; plot(cf,'.'); drawnow

% average new disparity/dissimilarity
df=(fsol(13,~g_trk)*G_sol(:,~g_trk)')...
    ./sum(G_sol(:,~g_trk),2)'; % mean

% figure; plot(df,'.'); drawnow

%-- 4) Greedy optimizization to keep valid subset of feasible solutions

% initiate open set
g_ope=true(1,size(S_sol,2)); % all

% initiate solution set
g_sol=false(1,size(S_sol,2)); % none

% % iterate cost function
% while nnz(g_ope)~=0
%     
%     % reward function free coordinates
%     rf=sum(S_sol(sum(S_sol(:,g_sol),2)==0,:),1);
%     
%     % figure; plot(rf,'.'); drawnow
    
    % initiate indices open solution set
    sol=find(g_ope);
    
    % sort cost function to priotitize in descending order: -n_min(sol)' n_max(sol)' rf(sol)' 
    [~,ind]=sortrows([wf(sol)' -ef(sol)' cf(sol)' -df(sol)' fliplr(1:length(sol))'],'descend');
    
    % figure; plot(cf(ind),'.'); drawnow
    
    % sort solution to cost
    sol=sol(ind);
    
    % solve optimal cost by maximize the sortation (T_sol is important here)
    [val,subsol]=max( double( T_sol(:,sol) >0.5 )*spdiags(fliplr(1:length(sol))',0,length(sol),length(sol)) ,[],2); % first unique el
    
    % keep unique nonzero solution
    sol=sol(unique(subsol(val>0.5)));
    
    % sort cost function to priotitize in descending order:  -n_min(sol)' n_max(sol)' 
    [~,ind]=sortrows([wf(sol)' -ef(sol)' cf(sol)' -df(sol)' fliplr(1:length(sol))'],'descend');
    
    % figure; plot(cf(ind),'.'); drawnow
    
    % sort solution to cost
    sol=sol(ind);
    
    % solve optimal cost by maximize the sortationS_sol(:,sol)'*
    [val,subsol]=max( double( S_sol(:,sol) >0.5 )*spdiags(fliplr(1:length(sol))',0,length(sol),length(sol)) ,[],2); % first unique el
    
    % keep unique nonzero solution
    sol=sol(unique(subsol(val>0.5)));
    
    % update solution
    g_sol(sol)=1; % =g_sol | ismember(ufl,ufl_sort(sol));
    
%     % get total solution
%     sol=find(g_sol);
%     
%     % sort solution to optimal subtriangulation
%     [~,ind]=sortrows([ -n_min(sol)' n_max(sol)' -ef(sol)' -df(sol)' fliplr(1:length(sol))' ],'descend'); %rf(sol)' ./cf(sol)'
%     
%     % sort solution to cost
%     sol=sol(ind);
%     
%     % solve optimal cost by maximize the sortation
%     [val,subsol]=max( double( S_sol(:,sol) >0.5 )*spdiags(fliplr(1:length(sol))',0,length(sol),length(sol)) ,[],2); % first unique el
%     
%     % keep unique nonzero solution
%     sol=ismember(find(g_sol),sol(unique(subsol(val>0))));
%     
%     % solution
%     g_sol(g_sol)=sol;
    
%     % find g_ope & complete overlapping matched setssum(M_sort'*M_sort(:,g_sol)==sum(M_sort,1)',2)'>0.5;%
%     g_ope=sum(T_sol(sum(T_sol(:,g_sol),2)==0,:),1)>0.5 ... this conditions is enough
%         & sum(S_sol(sum(S_sol(:,g_sol),2)==0,:),1)>0.5; % sum( (S_sol(:,g_sol)'*S_sol)==sum(S_sol,1) ,1)>0.5;%
    
% end

% keep valid solutions
fsol=fsol(:,ismember(fsol(1,:),ufl(g_sol)));

%-- #) Write data

% feasible solutions
Fsol=cat(2,Fsol,fsol);

%-- #) End message

disp(['Optimized ',num2str(nnz(g_sol)),' new out of ',num2str(length(ufl)),' feasible solutions'])

%% plot results

%-- #) Make figure

if strcmp(plotting,'on')
    
    figure(prop.res(4)+2)
    
    cla
%     hold on
    scatter3(Fsol(6,:),Fsol(7,:),Fsol(8,:),[],Fsol(4,:),'.')
    axis equal
    view(2)
    hold off
    title(['frame ',num2str(n)])
    colorbar
    xlim([-3 3]); ylim([11 17])
    
    drawnow
    
    figure(prop.res(4)+3)
    
    cla
%     hold on
    [~,~,col]=unique(Fsol(1,:));
    scatter3(Fsol(6,:),Fsol(7,:),Fsol(8,:),[],col,'.')
    colormap lines
    axis equal
    view(2)
    hold off
    title(['frame ',num2str(n)])
    xlim([-3 3]); ylim([11 17])
    
    drawnow
    
end

end


%     % sort variables according to cost
%     S_sort=S_sol(:,ind);
%     ufl_sort=ufl(ind);
%     g_sort=g_ope(ind); % open sortation
%     
%     %figure; spy(T_sort)
%     %figure; plot(sum(T_sort>0.5,1),'.'); drawnow
%     
%     % initiate solution
%     sol=find(g_sort);
    
% % unique tracks
% uc=unique(cpln(1,:)); % [fsol(5,:)]

% % track solution mapping
% [ucp,~,cp]=unique(cpln(1, : ));
% iucp=find(ismember(uc,ucp));
% cp=iucp(cp);
% C_pln=sparse(cp,1:length(cp),ones(size(cp)),length(uc),length(cp));
% 
% %figure; spy(C_pln)

% % track solution mapping
% [ufl,~,fl]=unique(fsol(1, ismember(fsol(5,:),uc) )); % 
% [ucl,~,cl]=unique(fsol(5, ismember(fsol(5,:),uc) )); %double(~isnan(sum(fsol(:,ismember(fsol(5,:),uc) ),1)))
% iucl=find(ismember(uc,ucl));
% cl=iucl(cl);
% F_sol=sparse(cl,fl,ones(size(fl)),length(uc),length(ufl));%
% 
% %figure; spy(F_sol)
% %figure; plot(sum(F_sol,1),'.')


% % sort cost function to priotitize in descending order: -df'
% [~,ind]=sortrows([ -n_min' n_max' cf' -ef' -df' fliplr(1:length(df))' ],'descend');
% 
% % figure; plot(cf(ind),'.'); drawnow
% 
% % sort variables according to cost
% S_sort=S_sol(:,ind);
% T_sort=T_sol(:,ind);
% ufl_sort=ufl(ind);
% 
% %figure; spy(T_sort)
% %figure; plot(sum(T_sort>0.5,1),'.'); drawnow
% 
% % initiate open set
% g_ope=true(1,size(S_sort,2)); % all
% 
% % initiate solution set
% g_sol=false(1,size(S_sort,2)); % none
% 
% % iterate cost function
% while nnz(g_ope)~=0
%     
%     % initiate solution
%     sol=find(g_ope);
%     
%     % solve optimal cost by maximize the sortation
%     [val,subsol]=max( double( S_sort(:,sol) >0.5 )*spdiags(fliplr(1:length(sol))',0,length(sol),length(sol)) ,[],2); % first unique el
%     
%     % keep unique nonzero solution
%     sol=sol(unique(subsol(val>0.5)));
%     
%     % solve optimal cost by maximize the sortation
%     [val,subsol]=max( double( S_sort(:,sol)'*S_sort(:,sol) >0.5 )*spdiags(fliplr(1:length(sol))',0,length(sol),length(sol)) ,[],2); % first unique el
%     
%     % keep unique nonzero solution
%     sol=sol(unique(subsol(val>0.5)));
%     
%     % update solution
%     g_sol=g_sol | ismember(ufl_sort,ufl_sort(sol));
%     
%     % find g_ope & complete overlapping matched setssum(M_sort'*M_sort(:,g_sol)==sum(M_sort,1)',2)'>0.5;%
%     g_ope=sum(T_sort(sum(T_sort(:,g_sol),2)==0,:),1)>0.5 ... this conditions is enough
%         & sum(S_sort(sum(S_sort(:,g_sol),2)==0,:),1)>0.5; %sum( (S_sol(:,g_sol)'*S_sol)==sum(S_sol,1) ,1)>0.5;%
%     
% end
% 
% % keep valid solutions
% fsol=fsol(:,ismember(fsol(1,:),ufl_sort(g_sol)));

