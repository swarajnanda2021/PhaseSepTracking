function postproc_lyapanalysis
%postproc_lyapanalysis Proces Lyapunov finite time separation analysis
%       
%   See how pairs in promity come together and separte when been in
%   proximity
%
%   Note that simply computing all pairs not neccesary in proximity blow up
%   in numbers

%% get globals
global folder date rec post % prop traj pair

%% Create memory maps
PairIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Pairs' vsl 'PairIndex.dat']);
TrajIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
TrajPosition=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
TrajVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
TrajAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);
TrajCurve=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Curve.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis'])
end

%% Initiate files to write
LyapIndex=zeros(4,0); % pairindices [p-index n-frame i-index j-index]
LyapPosition=zeros(3,0); % pair mid position 
LyapVelocity=zeros(3,0); % pair mid velocity 
LyapAccelaration=zeros(3,0); % pair mid accelaration 
LyapCurve=zeros(size(TrajCurve,1),0); % pair mid curve
LyapDistance=zeros(3,0); % pair separation
LyapVelocityDifference=zeros(3,0); % pair separation velocity
LyapAccelarationDifference=zeros(3,0); % pair separation accelaration
LyapFluctuation=zeros(size(TrajCurve,1),0); % pair separation curve

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for n=post.tproc(1):post.tproc(2) % post.tproc(1)+(cwin):post.tproc(2)-(cwin+mran)
    
    % get data frame range n
    N= TrajIndex(2,:)==n;
    
    % get the indexing
    trajindex=TrajIndex(1,N);
    
    % get the positions
    trajpos=TrajPosition(:,N);
    
    % get the velocity
    trajvel=TrajVelocity(:,N);
    
    % get the positions
    trajacc=TrajAccelaration(:,N);
    
    % get the positions
    trajcurve=TrajCurve(:,N);
    
    %figure; plot3(trajpos(1,:),trajpos(2,:),trajpos(3,:),'.'); axis equal;   xlim([-2 2]); ylim([10 16]); zlim([-4 0.5]); view(2)
    
    % get previous pair indices
    pairindex=PairIndex(:, ismember(PairIndex(3,:),trajindex(1,:))...
        & ismember(PairIndex(4,:),trajindex(1,:)));
    
    % keep unique part
    pairindex=unique(pairindex([1 3 4],:)','rows')';
    
    %-- #) Pair adjacency
    
    % find unique
    uc=unique([trajindex pairindex(2,:) pairindex(3,:)]);
    
    % compute metric distance midpoint tracks on pairs
    [ucp,~,ip]=unique(pairindex(2,:)); %
    iucp=find(ismember(uc,ucp));
    ip=iucp(ip);
    AdjI= sparse(ip,(1:length(ip)),ones(1,length(ip)),length(uc),length(ip)) ;
    [ucp,~,jp]=unique(pairindex(3,:)); %
    iucp=find(ismember(uc,ucp));
    jp=iucp(jp);
    AdjJ= sparse(jp,(1:length(jp)),ones(1,length(jp)),length(uc),length(jp)) ;
    
    %figure; spy(AdjI)
    %figure; spy(AdjJ)
    
    %-- #) compute usefull pair quantities
    
    % pair indexing
    lyapindex=[pairindex(1,:)
                n*ones(1,size(pairindex,2))
                pairindex([2 3],:)];
    
    % compute pair (mean) position
    lyappos=(trajpos*AdjJ+trajpos*AdjI)/2;
    
    %figure; plot3(lyappos(1,:),lyappos(2,:),lyappos(3,:),'.'); axis equal;   xlim([-2 2]); ylim([10 16]); zlim([-4 0.5]); view(2)
    
    % compute Lyap velocity
    lyapvel=(trajvel*AdjJ+trajvel*AdjI)/2;
    
    % compute Lyap accelaration
    lyapacc=(trajacc*AdjJ+trajacc*AdjI)/2;
    
    % compute Lyap curve
    lyapcurve=(trajcurve*AdjJ+trajcurve*AdjI)/2;
    
    % compute Lyap distance
    lyapdis=(trajpos*AdjJ-trajpos*AdjI);
    
    % compute Lyap Velocity difference
    lyapveldiff=(trajvel*AdjJ-trajvel*AdjI);
    
    % compute Lyap accelaration difference
    lyapaccdiff=(trajacc*AdjJ-trajacc*AdjI);
    
    % compute Lyap fluctuation curve
    lyapfluc=(trajcurve*AdjJ-trajcurve*AdjI);
    
    %dum1=sqrt(sum(lyapdis.^2,1));
    %dum2=sqrt(sum((lyapveldiff./lyapvel).^2,1));
    %figure; plot(dum1,dum2,'.'); 
    
    %-- #) write Lyaps
    
    % write pairs
    LyapIndex=cat(2,LyapIndex,lyapindex);
    
    % write Lyap positions
    LyapPosition=cat(2,LyapPosition,lyappos);
    
    % write Lyap curve
    LyapCurve=cat(2,LyapCurve,lyapcurve);
    
    % write Lyap velocity
    LyapVelocity=cat(2,LyapVelocity,lyapvel);
    
    % write Lyap accelaration
    LyapAccelaration=cat(2,LyapAccelaration,lyapacc);
    
    % write Lyap distance
    LyapDistance=cat(2,LyapDistance,lyapdis);
    
    % write Lyap velocity difference
    LyapVelocityDifference=cat(2,LyapVelocityDifference,lyapveldiff);
    
    % write Lyap velocity difference
    LyapAccelarationDifference=cat(2,LyapAccelarationDifference,lyapaccdiff);
    
    % write Lyap velocity difference
    LyapFluctuation=cat(2,LyapFluctuation,lyapfluc);
    
    % message
    disp(['Post processed lyapunov analysis in frame ',num2str(n)])
    
end

toc

%% save

% Pairindexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapIndex.dat'],LyapIndex,'w');

% Pair position
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapPosition.dat'],LyapPosition,'w');

% Pair Curve
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapCurve.dat'],LyapCurve,'w');

% Pair velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapVelocity.dat'],LyapVelocity,'w');

% Pair accelaration
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapAccelaration.dat'],LyapAccelaration,'w');

% Pair distance
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapDistance.dat'],LyapDistance,'w');

% Pair velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapVelocityDifference.dat'],LyapVelocityDifference,'w');

% Pair Accelaration
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapAccelarationDifference.dat'],LyapAccelarationDifference,'w');

% Pair Fluctuation
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapFluctuation.dat'],LyapFluctuation,'w');

end
