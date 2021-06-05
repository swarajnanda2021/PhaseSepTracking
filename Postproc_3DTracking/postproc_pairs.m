function postproc_pairs
%postproc_pairs Post process pairs for pair analysis
%   
%   Define [spatio-temporal] pairs for
%       lyapunov analysis
%       pair corrolation
%       information tansfer
%   
%   Note: Omitted a pair frame, which is ill-define; even when contraining
%       gravity vector.

%% get globals
global folder date rec post pair % prop

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);
Curve=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Curve.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs'])
end

%% Initiate files to write
PairIndex=zeros(1+1+2,0); % pairindices [p-index n-frame i-index j-index]
PairPosition=zeros(3,0); % pair mid position 
PairVelocity=zeros(3,0); % pair mid velocity 
PairAccelaration=zeros(3,0); % pair mid accelaration 
PairCurve=zeros(size(Curve,1),0); % pair mid curve
PairDistance=zeros(3,0); % pair separation
PairVelocityDifference=zeros(3,0); % pair separation velocity
PairAccelarationDifference=zeros(3,0); % pair separation accelaration
PairFluctuation=zeros(size(Curve,1),0); % pair separation curve

%% post processing
disp('post proc')
tic

% intiate pair indexing
pmax=0;

% loop over the matched and indexed Plinks
for n=post.tproc(1):post.tproc(2) % post.tproc(1)+(cwin):post.tproc(2)-(cwin+mran)
    
    %-- #) Get data
    
    % get data frame range n
    N= Index(2,:)==n;
    
    % get the indexing
    ind=Index(1:2,N);
    
    % get the positions
    pos=Position(:,N);
    
    % get the velocity
    vel=Velocity(:,N);
    
    % get the positions
    acc=Accelaration(:,N);
    
    % get the positions
    curve=Curve(:,N);
    
    %figure; plot3(pos(1,:),pos(2,:),pos(3,:),'.'); axis equal;   xlim([-2 2]); ylim([10 16]); zlim([-4 0.5]); view(2)
    
    % get previous pair indices
    pairprev=PairIndex(:,PairIndex(2,:)==(n-1) ...
        & ismember(PairIndex(3,:),ind(1,:))...
        & ismember(PairIndex(4,:),ind(1,:)));
    
    %-- #) Test length non continuing pair indexing
    
    % get not selected previous
    pairtest=PairIndex(:,ismember(PairIndex(1,:),...
        PairIndex(1, PairIndex(2,:)==(n-2))) ...
        & ~ismember(PairIndex(1,:),pairprev(1,:)));
    
    % test length
    if ~isempty(pairtest)
        up=unique(pairtest(1,:));
        cnt=histcounts(pairtest(1,:),[up-1/2 max(up)+1/2]);
        rem=cnt<pair.pfrm;
        G_keep=~ismember(PairIndex(1,:),up(rem));
    else
        G_keep=true(1,size(PairIndex,2));
    end
    
    %-- #) Pair adjacency
    
    % define unique track indices
    uc=unique(ind(1,:));
    
    % define unique pairs posI posJ
    [pairI,pairJ]=ndgrid(ind(1,:),ind(1,:)); % unique(indI),unique(indJ));
    K=pairI>pairJ;
    pairIJ=[reshape(pairI(K),1,[])
        reshape(pairJ(K),1,[])]; % keep unique upper traingulator part (also removing self pair)
%     pairIJ=unique([pairI(:) pairJ(:)],'rows')';
    
    % compute metric distance midpoint tracks on pairs
    [uip,~,ip]=unique(pairIJ(1,:)); %
    iuc=find(ismember(uc,uip));
    ip=iuc(ip);
    AdjI= sparse(ip,(1:length(ip)),ones(1,length(ip)),length(uc),length(ip)) ;
    [ujp,~,jp]=unique(pairIJ(2,:)); %
    juc=find(ismember(uc,ujp));
    jp=juc(jp);
    AdjJ= sparse(jp,(1:length(jp)),ones(1,length(jp)),length(uc),length(jp)) ;
    
    %figure; spy(AdjI)
    %figure; spy(AdjJ)
    
    % segment pairs in metric range @ time consuption
    segIJ=sqrt(sum((pos*AdjJ-pos*AdjI).^2,1))<=pair.pdis; % segmentation selection
    pairIJ=pairIJ(:,segIJ); % pairs in range
    
    % update pair adjacency
    AdjI=AdjI(:,segIJ);
    AdjJ=AdjJ(:,segIJ);
    
    %figure; spy(AdjI)
    %figure; spy(AdjJ)
    
    % define previous pairs adjacency
    [uip,~,ip]=unique(pairprev(3,:));
    iuc=find(ismember(uc,uip));
    ip=iuc(ip);
    [ujp,~,jp]=unique(pairprev(4,:));
    juc=find(ismember(uc,ujp));
    jp=juc(jp);
    AdjP=sparse(ip,1:length(ip),ones(size(ip)),length(uc),length(ip)) ...
        + sparse(jp,1:length(jp),ones(size(jp)),length(uc),length(jp)) ;
    
    %figure; spy(AdjP)
    
    % relate pairs to previous frame else [re-]define
    B=( AdjP'*(AdjI+AdjJ) )==2;
    
    % find
    [bi,bj]=find(B);
    
    % pair indexing
    pairind=[zeros(1,size(pairIJ,2))
                n*ones(1,size(pairIJ,2))
                pairIJ];
    
    % copy indexing
    pairind(1,bj)=pairprev(1,bi);
    
    % initiate new indexing
    i_new=1:nnz(pairind(1,:)==0);
    pairind(1,pairind(1,:)==0)=pmax+i_new;
    pmax=pmax+length(i_new);
    
    %-- #) compute usefull pair quantities
    
    % compute pair (mean) position
    pairpos=(pos*AdjJ+pos*AdjI)/2;
    
    %figure; plot3(pairpos(1,:),pairpos(2,:),pairpos(3,:),'.'); axis equal;   xlim([-2 2]); ylim([10 16]); zlim([-4 0.5]); view(2)
    
    % compute pair velocity
    pairvel=(vel*AdjJ+vel*AdjI)/2;
    
    % compute pair accelaration
    pairacc=(acc*AdjJ+acc*AdjI)/2;
    
    % compute pair curve
    paircurve=(curve*AdjJ+curve*AdjI)/2;
    
    % compute pair distance
    pairdis=(pos*AdjJ-pos*AdjI);
    
    % compute pair Velocity difference
    pairveldiff=(vel*AdjJ-vel*AdjI);
    
    % compute pair accelaration difference
    pairaccdiff=(acc*AdjJ-acc*AdjI);
    
    % compute pair fluctuation curve
    pairfluc=(curve*AdjJ-curve*AdjI);
    
    %dum1=sqrt(sum(pairdis.^2,1));
    %dum2=sqrt(sum((pairveldiff./pairvel).^2,1));
    %figure; plot(dum1,dum2,'.'); 
    
    %-- #) write pairs
    
    % write pairs
    PairIndex=cat(2,PairIndex(:,G_keep),pairind);
    
    % write pair positions
    PairPosition=cat(2,PairPosition(:,G_keep),pairpos);
    
    % write pair curve
    PairCurve=cat(2,PairCurve(:,G_keep),paircurve);
    
    % write pair velocity
    PairVelocity=cat(2,PairVelocity(:,G_keep),pairvel);
    
    % write pair accelaration
    PairAccelaration=cat(2,PairAccelaration(:,G_keep),pairacc);
    
    % write pair distance
    PairDistance=cat(2,PairDistance(:,G_keep),pairdis);
    
    % write pair velocity difference
    PairVelocityDifference=cat(2,PairVelocityDifference(:,G_keep),pairveldiff);
    
    % write pair velocity difference
    PairAccelarationDifference=cat(2,PairAccelarationDifference(:,G_keep),pairaccdiff);
    
    % write pair velocity difference
    PairFluctuation=cat(2,PairFluctuation(:,G_keep),pairfluc);
    
    % message
    disp(['Post processed pairs in frame ',num2str(n)])
    
end

toc

%% Save pairs

% Pairindexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'PairIndex.dat'],PairIndex,'w');

% Pair position
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'PairPosition.dat'],PairPosition,'w');

% Pair Curve
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'pairCurve.dat'],PairCurve,'w');

% Pair velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'pairVelocity.dat'],PairVelocity,'w');

% Pair accelaration
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'pairAccelaration.dat'],PairAccelaration,'w');

% Pair distance
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'pairDistance.dat'],PairDistance,'w');

% Pair velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'PairVelocityDifference.dat'],PairVelocityDifference,'w');

% Pair Accelaration
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'PairAccelarationDifference.dat'],PairAccelarationDifference,'w');

% Pair Fluctuation
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'PairFluctuation.dat'],PairFluctuation,'w');

end

