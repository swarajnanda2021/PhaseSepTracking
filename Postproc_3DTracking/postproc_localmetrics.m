function postproc_localmetrics
%postproc_localmetrics Summary of this function goes here
%   Detailed explanation goes here

global folder date rec post locm

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'Localmetrics'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'Localmetrics'])
end

%% Initiate files to write
Ranking=zeros(locm.nngb,0);
Distance=zeros(locm.nngb,0);
VelocityDifference=zeros(locm.nngb,0);
Polarization=zeros(locm.nngb,0);
% Milling=zeros(locm.nngb,0);

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for n=post.tproc(1):post.tproc(2) %unique(Index(2,:))
    % get trajectory of interest
    N=Index(2,:)==n ; % rem nans 
    
    % indexing
    index=Index(1,N);
    
    % get position in time
    pos=Position(:,N);
    
    % get velocity in time
    vel=Velocity(:,N);
    
    % initiate index relations
    I=repmat(index',1,length(index));
    
    % compute euclidean distances
    D=cellfun(@(x,y)norm(x-y,2),repmat(num2cell(pos,1)',1,size(pos,2)),...
        repmat(num2cell(pos,1),size(pos,2),1));
    
    % compute velocity difference
    V=cellfun(@(x,y)norm(x-y,2),repmat(num2cell(vel,1)',1,size(vel,2)),...
        repmat(num2cell(vel,1),size(vel,2),1));
    
    % compute polarization
    P=cellfun(@(x,y)x'*y/norm(x,2)/norm(y,2),repmat(num2cell(vel,1)',1,size(vel,2)),...
        repmat(num2cell(vel,1),size(vel,2),1));%1-acos(/pi)
    
    % treshold
    [D,i]=sort(D,1);
    I=I(i);
    V=V(i);
    P=P(i);
    
    % remove self 
    D(1,:)=[];
    I(1,:)=[];
    V(1,:)=[];
    P(1,:)=[];
    
    % remove outside number neighbors
    D=D(1:locm.nngb,:);
    I=I(1:locm.nngb,:);
    V=V(1:locm.nngb,:);
    P=P(1:locm.nngb,:);
    
    %%% write variables beneath
    
    Ranking=cat(2,Ranking,I);
    Distance=cat(2,Distance,D);
    VelocityDifference=cat(2,VelocityDifference,V);
    Polarization=cat(2,Polarization,P);

    % message
    disp(['post processed local metrics frame ',num2str(n)])
end

toc

%% save

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Localmetrics' vsl 'Index.dat'],Index,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Localmetrics' vsl 'Ranking.dat'],Ranking,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Localmetrics' vsl 'VelocityDifference.dat'],VelocityDifference,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Localmetrics' vsl 'Distance.dat'],Distance,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Localmetrics' vsl 'Polarization.dat'],Polarization,'w');

end