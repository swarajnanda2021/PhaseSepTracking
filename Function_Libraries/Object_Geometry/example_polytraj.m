%% Notes
%
%   Author: Koen Muller, 06 June 2020
%
% example script for fitting a polynomial trajectories to position data
%
% this function allows testing 2D and 3D data
%
% To high order trajectories for example on a 15 stencil and 14 order,
%   oscilate to much and loose accuracy while fitted at low residuals.
%   This results in ill-conditioning warning aswell.

%% Start
close all
clear all
clc

%% Controls

% trajectories
I=100; % number of trajectories
N=100; % number frames [even number]
tini='random'; % random or fullspan start/end [in]complete data
smth=0.75;% smoothness in 0=random, 1=smooth

% fit
tfit=15; % trajectorie stencil size
tord=3; % 0 | 1 | 2 | 3 trajectory order

%% Define some synthetic trajectories

% define trajectory data
tdat=zeros(5,0);
for i=1:I
    
    % initiate frames
    switch tini
        case 'random'
            nspan=sort([-randi(N/2,1,1) randi(N/2,1,1)]); % random start and end
        case 'full'
            nspan=[-N/2 N/2]; % full span time series
    end
    
    % initial condition position
    x=(2*rand(3,1)-1);
    
    % initial direction
    d=(2*rand(3,1)-1); 
    
    % loop span frames
    for n=nspan(1):nspan(2)
        
        % new direction
        d=smth*d+(1-smth)*(2*rand(3,1)-1);
        
        % new position
        x=x+2/N*d; % displacing one unit over all frames
        
        % append track
        tdat=cat(2,tdat,[i;n;x]);
        
    end
    
end

% message
disp(['generated ',num2str(I),' trajectories spanning ',...
    num2str(N), ' frames with a ',tini,'span track'])

%% Plot 3D trajectories

% indexing
figure
scatter3(tdat(3,:),tdat(4,:),tdat(5,:),[],tdat(1,:),'.')
colormap lines
axis tight
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectories Indexing')
drawnow

% timeseries
figure
scatter3(tdat(3,:),tdat(4,:),tdat(5,:),[],tdat(2,:),'.')

% layout
colormap parula
col=colorbar;
col.Label.String = 'frame number';
col.FontSize=11;
axis tight
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectories Timeseries')
drawnow

%% evaluate fit trajecoties over time

% loop frames
idat=zeros(2,0); % index
xdat=zeros(size(tdat(3:end,:),1),0); % position
vdat=zeros(size(tdat(3:end,:),1),0); % velocity
adat=zeros(size(tdat(3:end,:),1),0); % accelaration
rdat=zeros(1,0); % residuals
for n=unique(tdat(2,:))
    
    % select time frames
    N=ismember(tdat(1,:),tdat(1,tdat(2,:)==n)) ...
        & tdat(2,:)>=n-floor(tfit/2) & tdat(2,:)<=n+floor(tfit/2);
    
    % fit data
    fdat=tdat(:,N);
    
    % get coeffients
    [coef ,res]=polytraj(fdat,tord); %
    
    %dum=polyeval(coef,fdat(1:2,:),0);
    %figure; scatter3(dum(1,:),dum(2,:),dum(3,:),[],res,'.')
    
    % compute total residual track
    [uf,~,fi]=unique(fdat(1,:));
    A=sparse(fi,1:length(fi),ones(size(fi)),length(uf),length(fi));
    
    % figure; spy(A)
    
    % mean residuals
    res=res*A'./sum(A,2)' ;
    
    % evaluate fittered positions
    x=polyeval(coef,tdat(1:2,tdat(2,:)==n),0); % position
    v=polyeval(coef,tdat(1:2,tdat(2,:)==n),1); % velocity
    a=polyeval(coef,tdat(1:2,tdat(2,:)==n),2); % accelaration
    
    % write
    idat=cat(2,idat,tdat(1:2,tdat(2,:)==n));
    xdat=cat(2,xdat,x);
    vdat=cat(2,vdat,v);
    adat=cat(2,adat,a);
    rdat=cat(2,rdat,res);
    
end


%% Plot 3D trajectories

% trajectorie and fit
figure
plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.')
hold on
scatter3(xdat(1,:),xdat(2,:),xdat(3,:),[],rdat(1,:),'.')
hold off
col=colorbar;
col.Label.String = 'residuals' ;
col.FontSize=11;
axis tight
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectories Fit')
drawnow

% velocity
figure
scatter3(xdat(1,:),xdat(2,:),xdat(3,:),[],sqrt(sum(vdat.^2,1)),'.')

% layout
colormap parula
col=colorbar;
col.Label.String = 'velocity' ;
col.FontSize=11;
axis tight
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectories Velocity')
drawnow

% velocity
figure
scatter3(xdat(1,:),xdat(2,:),xdat(3,:),[],sqrt(sum(adat.^2,1)),'.')

% layout
colormap parula
col=colorbar;
col.Label.String = 'accelaration';
col.FontSize=11;
axis tight
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectories Accelaration')
drawnow




