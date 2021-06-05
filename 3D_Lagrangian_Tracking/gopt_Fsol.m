function Fsol=gopt_Fsol(Fsol,frames)
%gopt_Fsol Global optimization feasible solutions.
%
%   Input:
%       Fsol		Feasible solutions
%       frames 		Frames set
%   
%   Output:
%       Fsol 		Feasible solutions.
%   
%   working principle:
%       Greedy solution by cost function sortation.
%   	Fitting best tracks to image identifications.
%	


%% Get globals
global plotting % ctrl prop

%% Branch and optimize feasible solutions

%-- 0) Begin message

disp(['Optimize global feasible solutions frame ',num2str(min(frames)),' to ',num2str(max(frames))])

%-- #) Initiate data

% get feasible solution
fsol=Fsol( : , ismember( Fsol(1,:) , Fsol(1,ismember(Fsol(5,:),frames))));

%--#) Open for reaccessment

% feasible solutions
Fsol=Fsol( : , ~ismember( Fsol(1,:) , fsol(1,:) ));

%-- #) Define adjacencies

% define feasible solution adjacency
[ufl,~,fl]=unique(fsol(1, : ));%ismember(fsol(1,:),ufl)
G_sol=sparse(fl,1:length(fl),ones(size(fl)),length(ufl),length(fl));

% track set

%figure; spy(G_sol)

% track solution mapping
[ufl,~,fl]=unique(fsol(1, : )); % g_prd
[usl,~,sl]=unique(fsol(3, : )); %double(~isnan(sum(fsol,1)))
S_sol=sparse(sl,fl,ones(size(sl)),length(usl),length(ufl));%

%figure; spy(S_sol)
%figure; plot(sum(S_sol>0.5,1),'.')
%figure; plot(sum(S_sol>0.5,2),'.')

%-- #) Initiate optimization

% track start
[~,dum]=max(G_sol*spdiags(max(fsol(5,:))+1-fsol(5,:)',0,size(fsol,2),size(fsol,2)),[],2);
n_min=fsol(5,dum);

%figure; plot(n_min,'.'); drawnow

% track end
n_max=max(G_sol*spdiags(fsol(5,:)',0,size(fsol,2),size(fsol,2)),[],2)';

%figure; plot(n_max,'.'); drawnow

% camera weight track
wf=sum( G_sol ,2)';

% figure; plot(wf,'.'); drawnow

% % average track error (sortation first, non-unique)g_fit & g_fit & g_fit &
% ef=((fsol(20,~g_frm))*G_sol(:,~g_frm)')...
%     ./sum(G_sol(:, ~g_frm),2)'; % mean

% figure; plot(ef,'.'); drawnow

% % camera cover new
% cf=sum( G_sol(:, g_frm ) ,2)'; %(:,g_trk)C_pln'*C_cam*average

% figure; plot(cf,'.'); drawnow

% average new disparity/dissimilarity
df=((fsol(20, : ))*G_sol')...
    ./sum(G_sol,2)'; % mean(:, g_frm )(:, g_frm ) g_frm

% figure; plot(df,'.'); drawnow

%-- 4) Greedy optimizization to keep valid subset of feasible solutions

% initiate open set
g_ope=true(1,size(S_sol,2)); %cf>=3;% minimum focality here all

% initiate solution set
g_sol=false(1,size(S_sol,2)); % none

% % iterate cost function
% while nnz(g_ope)~=0
%     
%     % reward function free coordinates
%     rf=sum(S_sol( sum( S_sol( : ,g_sol),2)==0 ,:),1);%ismember(usl,fsol(5,~g_trk))' & 
%     
%     % figure; plot(rf,'.'); drawnow
    
    % initiate indices open solution set
    sol=find(g_ope);
    
    % sort cost function to priotitize in descending order: rf(sol)'  cf(sol)' -nf(sol)' rf(sol)' -ef(sol)'   
    [~,ind]=sortrows([-n_min(sol)' n_max(sol)' wf(sol)' -df(sol)' fliplr(1:length(sol))'],'descend'); % make sure to select first starting last ending, to continue written tracks first
    
    % figure; plot(cf(ind),'.'); drawnow
    
    % sort solution to cost
    sol=sol(ind);
    
    % solve optimal cost by maximize the sortation (T_sol can be important here)
    [val,subsol]=max( double( S_sol(:,sol) >0.5 )*spdiags(fliplr(1:length(sol))',0,length(sol),length(sol)) ,[],2); % first unique el
    
    % solution
    sol=sol(unique(subsol(val>0)));%(sel));
    
    % update solution
    g_sol(sol)=1; % =g_sol | ismember(ufl,ufl_sort(sol));
    
%     % find g_ope & complete overlapping matched setssum(M_sort'*M_sort(:,g_sol)==sum(M_sort,1)',2)'>0.5;%
%     g_ope=g_ope & sum(T_sol(sum(T_sol(:,g_sol),2)==0,:),1)>0.5 ... by def this conditions is enough
%         & sum(S_sol(sum(S_sol(:,g_sol),2)==0,:),1)>0.5; % sum( (S_sol(:,g_sol)'*S_sol)==sum(S_sol,1) ,1)>0.5;%
%     
% end

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
    N=ismember(Fsol(5,:),frames);
    
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

