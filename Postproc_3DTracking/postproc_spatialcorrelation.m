function postproc_spatialcorrelation
%postproc_spatialcorrelation Post process the spatial 2-point correlation 
%   function.
%
%   Different things can be computed
%       Spatial correlation function
%       Spatial wavespectrum
%
%   Normalizations is choosen s.t.
%       The correlation is weighted by the volume of data intersection
%       The trace of the correlation and spectral matrix are at unity
%
%   Note, for now only velocity is implemented, but also think about the
%       polarization vector and for example the acelaration vector.

%% get globals
global folder date rec post prop eulr

%% load grid
vecgrid=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'vecgrid.mat']);

%% Get gridded data

% load general data
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Time.dat']);
% GridPosition=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Position.dat']);

% load different cases
% BinnedVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedVelocity.dat']);
% BinnedAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedAccelaration.dat']);
% Polarization=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Polarization.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Velocity.dat']);
% Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Accelaration.dat']);

%% make folder for data
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation'])
end

%% initiate data sets

% griddata
CorIndex=zeros(2,0); % index
CorTime=zeros(1,0); % time
CorDis=zeros(3,0); % distance
CorFun=zeros(6,0); % correlation function
CorWav=zeros(3,0); % wavenumber
CorSpe=zeros(12,0); % Spectrum (real and complex)
CorNor=zeros(1,0); % Normalization

%% loop frames in data, preserve time resolution

% loop frames
for n=post.tproc(1):post.tproc(2)%unique(Index(2,:))
    
    %-- Message
    disp(['Process spatial correlation function frame ',num2str(n)])
    
    %-- Get data
    
    % get frame
    N=Index(2,:)==n;
    
    % Indexing
    ind=Index(:,N);
    time=Time(:,N);

    % Get velocity data frame
    vel=Velocity(:,N);
    
    % Get grid positions
    pos=vecgrid.X(:,ind(1,:));
    
    % remove nans if there
    ind=ind(:,~isnan(sum(vel,1)));
    pos=pos(:,~isnan(sum(vel,1)));
    vel=vel(:,~isnan(sum(vel,1)));
    
    %figure; quiver3c(pos(1,:),pos(2,:),pos(3,:),vel(1,:),vel(2,:),vel(3,:),3)
    
    %-- Compute correlation
    
    % assemble on grid data
    grd=cat(1,vecgrid.IND,vecgrid.X);
    dat1=[ind(1,:)
            vel];
    dat2=dat1;
    
    % compute spatial properly normalized 2point correlation function
    [dis,cor,wav,spe,nor] = spatial2pointcor(grd,dat1,dat2,'normalized');
    
    % remove empty correlation and find indexing
    sel=~isinf(nor);
    ind=[find(sel)
        n*ones(1,nnz(sel))];
    time=unique(time)*ones(1,nnz(sel));
    dis=dis(:,sel);
    cor=cor(:,sel);
    wav=wav(:,sel);
    spe=spe(:,sel);
    nor=nor(:,sel);
    
    % select zero displacement
    sel=sum(dis==0,1)==3;
    
    % scalar normalization correlation
    cor0=sqrt(cor([1 4 6],sel)'*cor([1 4 6],sel)); % make it a unit vector
    cor=cor/cor0; % normalize scalar s.t. magnitude componentscan be compared
    
    % select zero wavenumber
    sel=sum(wav==0,1)==3;
    
    % kinetic energy by wavenumber vector
    spe0=sqrt(spe([1 4 6],sel)'*conj(spe([1 4 6],sel))); % make it a unit vector
    spe=spe/spe0; % normalize scalar s.t. magnitude componentscan be compared
    
    %figure; scatter3(dis(1,:),dis(2,:),dis(3,:),100,squeeze(cor(1,:)),'.'); colorbar;
    %figure; scatter3(dis(1,:),dis(2,:),dis(3,:),100,squeeze(cor(2,:)),'.'); colorbar;
    %figure; scatter3(dis(1, dis(3,:)==0 ),dis(2, dis(3,:)==0 ),dis(3, dis(3,:)==0 ),100,squeeze(cor(1, dis(3,:)==0 ) + cor(4, dis(3,:)==0 ) + cor(6, dis(3,:)==0 )),'.'); view(2); colorbar; 
    %figure; plot(dis(1, dis(2,:)==0 & dis(3,:)==0),squeeze(cor(1, dis(2,:)==0 & dis(3,:)==0)),'.'); colorbar;
    
    %-- Write correlation
    
    % write index
    CorIndex=cat(2,CorIndex,ind);
    
    % write time
    CorTime=cat(2,CorTime,time);
    
    % write distance
    CorDis=cat(2,CorDis,dis);
    
    % write correlation function
    CorFun=cat(2,CorFun,cor);
    
    % write wavenumber
    CorWav=cat(2,CorWav,wav);
    
    % write Spectrum
    CorSpe=cat(2,CorSpe,[real(spe)
                            imag(spe)]);
    
    % write Normalization
    CorNor=cat(2,CorNor,nor);
    
end

%% Save

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Index.dat'],CorIndex,'w');

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Time.dat'],CorTime,'w');

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Distance.dat'],CorDis,'w');

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Correlation.dat'],CorFun,'w');

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'WaveNumber.dat'],CorWav,'w');

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Spectrum.dat'],CorSpe,'w');

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Normalization.dat'],CorNor,'w');

end

