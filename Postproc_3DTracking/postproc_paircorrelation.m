function postproc_paircorrelation
%postproc_paircorrolatiom Proces pair correlation see how pairs from lyapunow
%   analysis corrolate on distance and approach/separation
%       
%   See how pairs in promity come together and separte when been in
%   proximity.
%
%   type of correlation functions
%       plain correlation
%       binned correlation
%       see Cavagna 2018 paper for more interesting correlations
%           (for sure the fluctuation pair correlation! e.g. using the eulerian reference vectors)
%
%   Note: Binning the data is commented out for now, because: Does it bring
%       new insight w.r.t. spatial velocity vector correlation(?), is it
%       not just a visualization step?

%% get globals
global folder date rec post eulr pair % prop 

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapIndex.dat']);% Pairindexing
% Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapPosition.dat']);% Pair position
%PairFrame=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapFrame.dat']);% Pair frame
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapVelocity.dat']);% Pair velocity
%Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapAccelaration.dat']);% Pair accelaration
%Curve=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapCurve.dat']);% Pair Curve
Distance=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapDistance.dat']);% Pair distance
VelocityDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapVelocityDifference.dat']);% Pair velocity
%AccelarationDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapAccelarationDifference.dat']);% Pair Accelaration
%Fluctuation=fload([folder date rec vsl 'Postproc3DTracking' vsl 'LyapunovAnalysis' vsl 'LyapFluctuation.dat']);% Pair Fluctuation

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'PairCorrelation'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'PairCorrelation'])
end

%% Define grids to bin data

% % edges for binning the data
% cedges=-1-1/(2*pair.sampl):1/pair.sampl:1+1/(2*pair.sampl); % Correlation edges
% dedges=-1/(2*pair.sampl):1/pair.sampl:max(sqrt(sum(Distance(1,:).^2,1)))+1/(2*pair.sampl); % distance edges
% 
% % meshed grid discetization variable
% [D,C]=ndgrid(dedges(1:end-1)+1/(2*pair.sampl),cedges(1:end-1)+1/(2*pair.sampl)); % grid

% % edges for spatial discretization
% xedges=round(min(Distance(1,:)))-eulr.gres/2:eulr.gres:round(max(Distance(1,:)))+eulr.gres/2; 
% yedges=round(min(Distance(2,:)))-eulr.gres/2:eulr.gres:round(max(Distance(2,:)))+eulr.gres/2; 
% zedges=round(min(Distance(3,:)))-eulr.gres/2:eulr.gres:round(max(Distance(3,:)))+eulr.gres/2; 
% 
% % meshed grid discetization variable
% [dX,dY,dZ]=ndgrid(xedges(1:end-1)+eulr.gres/2,...
%                     yedges(1:end-1)+eulr.gres/2,...
%                     zedges(1:end-1)+eulr.gres/2); % grid

% % wave number discretization
% kx=(-(length(xedges(1:end-1))-1)/2:(length(xedges(1:end-1))-1)/2)*(length(xedges(1:end-1))/eulr.gres); 
% ky=(-(length(yedges(1:end-1))-1)/2:(length(yedges(1:end-1))-1)/2)*(length(yedges(1:end-1))/eulr.gres);
% kz=(-(length(zedges(1:end-1))-1)/2:(length(zedges(1:end-1))-1)/2)*(length(zedges(1:end-1))/eulr.gres);
% 
% % meshed grid discetization variable
% [kX,kY,kZ]=ndgrid(kx,ky,kz); % grid

%% Initiate files to write
CorIndex=zeros(1+1+2,0);
CorDistance=zeros(3,0);
CorVelocity=zeros(6,0);

%SpeVelocity=zeros(1+1+2,0); %
%CorPolarization=zeros(1+1+2,0); %
%CorAccelaration=zeros(1+1+2,0); %

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for n=post.tproc(1):post.tproc(2) 
    
    %-- #) Get data
    
    % get previous pairs
    N=Index(2,:)==n;
    
    % get previous pair indices
    ind=Index(:,N);
    
    % get the positions
    dis=Distance(:,N);
    
    % get the velocity
    vel=Velocity(:,N);
    
    % get the velocity
    veldif=VelocityDifference(:,N);
    
    % vel index i
    veli=vel+veldif/2;
    
    % vel index j
    velj=vel-veldif/2; % check with pair definition
    
    %-- #) Compute Correlation
    
    % velocity Correlation
    velcor=[veli(1,:).*velj(1,:)
        veli(1,:).*velj(2,:)
        veli(1,:).*velj(3,:)
        veli(2,:).*velj(2,:)
        veli(2,:).*velj(3,:)
        veli(3,:).*velj(3,:)];
    
    % normalize Correlation
    velcor=velcor./(sqrt(sum(veli.^2,1)).*sqrt(sum(velj.^2,1)));
    
    %figure; plot(sqrt(sum(dis.^2,1)),velcor(1,:)+velcor(4,:)+velcor(6,:),'.')
    
    %%% accelaration Correlation etc..
    
    %-- #) Spatial binning the Correlation data
    
%     % intiate
%     spatcor_dens=zeros(size(dX,1),size(dY,2),size(dZ,3),size(velcor,1));
%     spatcor_mean=zeros(size(dX,1),size(dY,2),size(dZ,3),size(velcor,1));
%     spatcor_std=zeros(size(dX,1),size(dY,2),size(dZ,3),size(velcor,1));
%     
%     % spatially grid data
%     for i=1:size(velcor,1) % loop Correlation variables
%         for k=1:size(dZ,3)% loop z slices (since histcounts3 or n is not existing)
%             
%             % select data in slice k
%             sel=dis(3,:)>zedges(k) & dis(3,:)<=zedges(k+1);
%             
%             % make Correlation map
%             [map,~,~,binX,binY]=histcounts2(dis(1,sel),dis(2,sel),xedges,yedges);
%             
%             % make nan map if not data present
%             map(map==0)=nan;
%             
%             % remove data outside domain
%             rem=binX==0 | binY==0;
%             binX=binX(:,~rem);
%             binY=binY(:,~rem);
%             sel(sel)=~rem;
%             
%             % write point density map
%             spatcor_dens(:,:,k,i)=map;
%             
%             % get index
%             binI=sub2ind(size(map),binX,binY);
%             
%             % adjacencie
%             A=sparse(binI,1:length(binI),ones(size(binI)),numel(map),length(binI));
%             
%             %figure; spy(A)
%             
%             % binned Correlation
%             cormean=velcor(i,sel)*A'./sum(A,2)';
%             corstd=sqrt(velcor(i,sel).^2*A'./sum(A,2)'-cormean.^2);
%             
%             % reshape data
%             spatcor_mean(:,:,k,i)=reshape(cormean',size(map,1),size(map,2),size(cormean,1));
%             spatcor_std(:,:,k,i)=reshape(corstd',size(map,1),size(map,2),size(corstd,1));
%             
%         end % k
%     end % i
%         
%     %figure; scatter3(dX(:),dY(:),dZ(:),100,reshape(sum(spatcor_dens(:,:,:,[1 4 6]),4),[],1),'.'); axis equal; colorbar;
%     %figure; scatter3(dX(:),dY(:),dZ(:),100,reshape(sum(spatcor_mean(:,:,:,[1 4 6]),4),[],1),'.'); axis equal; colorbar;
%     %figure; errorbar(sqrt(dX(:).^2+dY(:).^2+dZ(:).^2),reshape(sum(spatcor_mean(:,:,:,[1 4 6]),4),[],1),reshape(sum(spatcor_std(:,:,:,[1 4 6]),4),[],1),'.');
    
%     % supress nans
%     spatcor_mean(isnan(spatcor_mean))=0;
%     spatcor_std(isnan(spatcor_mean))=0;
%     
%     % compute spectrum
%     corspec=fftn(spatcor_mean);
%     corspec=fftshift(corspec) % to be checked if usefull..
%     
%     %figure; plot(corspec(1,:)+corspec(4,:)+corspec(6,:),'.')
%     
%     % normalize spectrum
%     corspec=corspec./..;
    
    %-- #) Binning the Correlation value over distance 
    
%     % distance
%     disij=sqrt(sum(dis.^2,1));
%     
%     % trace of the Correlation -1 to +1
%     trcor=corij(1,:)+corij(4,:)+corij(6,:);
%     
%     % Q: What would be the eigenvalues? Or the Correlation matrix invariants?
%     
%     % make Correlation map
%     [map,~,~,binX,binY]=histcounts2(disij,trcor,dedges,cedges);
%     
%     %figure; surf(D',C',map'); view(2); shading interp
%     
%     % bin index
%     binI=1:length(binX);
%     
%     % remove not binned data
%     rem=binX==0 | binY==0;
%     binX=binX(~rem);
%     binY=binY(~rem);
%     binI=binI(~rem);
    
    %-- #) Distance Correlation
    
%     % adjacencie
%     A=sparse(binX,binI,ones(size(binX)),length(dedges)-1,length(disij));
%     
%     %figure; spy(A)
%     
%     % binned Correlation
%     bincor_mean=corij*A'./sum(A,2)';
%     bincor_std=sqrt(corij.^2*A'./sum(A,2)'-cormean.^2);
%     
%     %figure; plot(unique(D(:)),cormean(1,:)+cormean(4,:)+cormean(6,:),'.')
%     %figure; plot(unique(D(:)),corstd(1,:)+corstd(4,:)+corstd(6,:),'.')
%     %figure; errorbar(unique(D(:)),cormean(1,:)+cormean(4,:)+cormean(6,:),corstd(1,:)+corstd(4,:)+corstd(6,:),'.')
    
    %-- #) Write Correlation
    
    % Index Correlation
    CorIndex=cat(2,CorIndex,ind); % same as lyap index
    
    % write pair Correlation distance
    CorDistance=cat(2,CorDistance,dis);
    
    % write velocity pair Correlation
    CorVelocity=cat(2,CorVelocity,velcor);
    
%     % Spatial indexing
%     
%     % write binned spatial Correlation
%     SpatCorVelocity=cat(2,SpatCorVelocity,cormap); 
%     
%     % write binned spatial Correlation
%     SpeCorVelocity=cat(2,SpatCorVelocity,cormap); 
%     
%     % write pair positions
%     BinCorVelocity=cat(2,CorMapVelocity,[cormean
%                                             corstd]);
    
    % message
    disp(['Post processed pair correlation in frame ',num2str(n)])
    
end

toc

%% Save pair correlations

% Pairindexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'PairCorrelation' vsl 'CorIndex.dat'],CorIndex,'w');

% Pairindexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'PairCorrelation' vsl 'CorDistance.dat'],CorDistance,'w');

% Pair position
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'PairCorrelation' vsl 'CorVelocity.dat'],CorVelocity,'w');

end

