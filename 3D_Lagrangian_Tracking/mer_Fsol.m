function Fsol = mer_Fsol(Fsol,n)
%mer_Fsol Merge feasible solutions.
%
%   Input:
%       Fsol 	Feasible solutions
%       n 		Current frame.
%   
%   Output:
%       Fsol 	Feasible solutions that has merge sevaral branches into
%               one.
%   
%   working principle:
%       To track that have complete overlap in the camera plane merge, this
%       means either over segment full overlaps with the other.
%   

%% Get globals
global plotting % prop ctrl 

%% Initiate data

% maximum index
if isempty(Fsol)
    imaxfl=size(Fsol,2);
else
    imaxfl=max(Fsol(1,:));
end

%% Match nfocal sets

%-- #) Begin message

% keep sets complete in previous focality
disp(['Merge feasible tracks frame ',num2str(n)])

%-- #) Load data

% Past feasible solutions
fsolm=Fsol(:, ismember( Fsol(1,:) , Fsol(1, Fsol(5,:)==n ) ) ...
    &  ismember( Fsol(1,:) , Fsol(1, Fsol(5,:) < n ) ) ... % to remove frame merge
    & ~ismember( Fsol(1,:) , Fsol(1, Fsol(5,:) > n ) ) );

% Future feasible solutions
fsolp=Fsol(:, ismember( Fsol(1,:) , Fsol(1, Fsol(5,:)==n ) ) ...
    & ~ismember( Fsol(1,:) , Fsol(1, Fsol(5,:) < n ) ) );

%-- #) Open for reassessment

% feasible solution (for reindexing when tracks split and merge)
Fsol=Fsol( : , ~ismember( Fsol(1,:) , [fsolm(1,:) fsolp(1,:)] ));

%-- 2) define adjacencies

% camera indices
uc=unique([fsolm(3,:) fsolp(3,:)]);

% define feasible solution adjacency
[ufm,~,fm]=unique(fsolm(1,:));
G_solm=sparse(fm,1:length(fm),ones(size(fm)),length(ufm),length(fm));

g_solm=fsolm(5,:)==n;

%figure; spy(G_solm)

% track solution mapping
dum=find(g_solm);
[ufm,fmu,fm]=unique(fsolm(1, g_solm )); %ismember(fsol(5,:),uc) 
[ucm,~,cm]=unique(fsolm(3, g_solm )); % ismember(fsol(5,:),uc)
iucm=find(ismember(uc,ucm));
cm=iucm(cm);
F_solm=sparse(cm,fm,ones(size(fm)),length(uc),length(ufm));
fmu=dum(fmu);

f_summ=sum(F_solm,1);

%figure; spy(F_solm)
%figure; plot(f_solm,'.')

% define feasible solution adjacency
[ufp,~,fp]=unique(fsolp(1,:));
G_solp=sparse(fp,1:length(fp),ones(size(fp)),length(ufp),length(fp));

g_solp=fsolp(5,:)==n;

%figure; spy(G_solp)

% track solution mapping
[ufp,~,fp]=unique(fsolp(1, g_solp )); %ismember(fsol(5,:),uc) 
[ucp,~,cp]=unique(fsolp(3, g_solp )); % ismember(fsol(5,:),uc)
iucp=find(ismember(uc,ucp));
cp=iucp(cp);
F_solp=sparse(cp,fp,ones(size(fp)),length(uc),length(ufp));

f_sump=sum(F_solp,1);

%figure; spy(F_solp)
%figure; plot(f_solp,'.')

%-- #) Branch geometry and split solutions

% adjacency
A=F_solm'*F_solp;

%figure; spy(A)

% find
[ai,aj,ak]=find(A);

%figure; plot(ak,'.' )

% validity
val=ak'==f_sump(aj) | ak'==f_summ(ai) ; % & not! | no full overlap
ai=ai(val);
aj=aj(val);
ak=ak(val);

% redefine A
A=sparse(ai,aj,ak,size(A,1),size(A,2));

%figure; spy(A)

% define unconnected and new feasible solutions
g_unm=sum(A,2)==0;%true(1,size(A,2)); %
% g_unp=sum(A,1)==0;%true(1,size(A,1)); !% keep all in case merge is removed

%-- #) Write unbranched track data

% unmatched solution  
f_solm=fsolm(:,ismember(fsolm(1,:),ufm(g_unm))); % partial overlap but referenced base track
f_solp=fsolp;%(:,ismember(fsolp(1,:),ufp(g_unp))); % in full overlap but another track

%figure; scatter3(f_sol(14,:),f_sol(15,:),f_sol(16,:),[],f_sol(20,:),'.'); axis equal; view(2); drawnow;

%-- #) Write splitting tracks

% define re-indexing branching tracks
ifl_brn=imaxfl+(1:length(ak));
% imaxfl=imaxfl+length(ak);
itl_brn=fsolm(2,fmu(ai));

% find mapping to copy data
dum=find(~g_solm);
[gi,gj]=find( G_solm(ai,~g_solm) ); % always contain minus 1
gj=dum(gj);

%figure; spy(G_solm(aj,:))

% pairing temporal index corresponces
pf_solm=zeros( size(fsolm,1) , length(gj) );
pf_solm(1,:)=ifl_brn(gi); %
pf_solm(2,:)=itl_brn(gi); %
for i=3:size(pf_solm,1) % 3
    pf_solm(i,:)=fsolm(i,gj); %
end

%figure; scatter3(pf_solm(14,:),pf_solm(15,:),pf_solm(16,:),[],pf_solm(20,:),'.'); axis equal; view(2); drawnow;

% find current frame
dum1=find(g_solm);
dum2=find(f_summ(ai)'==ak & f_sump(aj)'~=ak);
[gi,gj]=find( G_solm( ai(f_summ(ai)'==ak & f_sump(aj)'~=ak) ,g_solm ) );
gj=dum1(gj);
gi=dum2(gi);

%figure; spy(G_solm(aj,:))

% pairing temporal index corresponces
nf_solm=zeros( size(fsolm,1) , length(gj) );
nf_solm(1,:)=ifl_brn(gi); %
nf_solm(2,:)=itl_brn(gi); %
for i=3:size(nf_solm,1) % 3
    nf_solm(i,:)=fsolm(i,gj); %
end

%figure; scatter3(pf_sol(5,:),pf_sol(6,:),pf_sol(7,:),[],pf_sol(8,:),'.'); axis equal; view(2); drawnow;

% find mapping to copy data
dum=find(~g_solp);
[gi,gj]=find( G_solp(aj,~g_solp) );
gj=dum(gj);

%figure; spy(G_sol(aj,:))

% pairing temporal index corresponces
pf_solp=zeros( size(fsolp,1) , length(gj) );
pf_solp(1,:)=ifl_brn(gi); %
pf_solp(2,:)=itl_brn(gi); %
for i=3:size(pf_solp,1)
    pf_solp(i,:)=fsolp(i,gj); %
end

%figure; scatter3(pf_solp(14,:),pf_solp(15,:),pf_solp(16,:),[],pf_solp(20,:),'.'); axis equal; view(2); drawnow;

% find current frame
dum1=find(g_solp);
dum2=find(f_sump(aj)'==ak);
[gi,gj]=find( G_solp( aj(f_sump(aj)'==ak) ,g_solp ) );
gj=dum1(gj);
gi=dum2(gi);

%figure; spy(G_solm(aj,:))

% pairing temporal index corresponces
nf_solp=zeros( size(fsolp,1) , length(gj) );
nf_solp(1,:)=ifl_brn(gi); %
nf_solp(2,:)=itl_brn(gi); %
for i=3:size(nf_solp,1) % 3
    nf_solp(i,:)=fsolp(i,gj); %
end

%figure; scatter3(pf_sol(5,:),pf_sol(6,:),pf_sol(7,:),[],pf_sol(8,:),'.'); axis equal; view(2); drawnow;

%-- #) Append

% write new tracking data
fsol=cat(2,f_solm,f_solp,pf_solm,nf_solm,pf_solp,nf_solp);%

%figure; plot3(fsol(14,:),fsol(15,:),fsol(16,:),'.'); axis equal; view(2); drawnow

%-- #) write

Fsol=cat(2,Fsol,fsol);

%-- #) End message

% keep sets complete in previous focality
disp(['Merged ',num2str(length(ai)),...
    ' new out of ',num2str(length(unique(fsol(1,:)))),...
    ' feasible solution frame ',num2str(n)])

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
    title(['Merged solution, cost value frame ',num2str(n)])
    colorbar
    drawnow
    
    % figure indexing
    figure
    scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(1,N),'.')
    colormap lines
    axis tight
    axis equal
    view(2)
    title(['Merged solution, indexing frame ',num2str(n)])
    drawnow
    
end

end

