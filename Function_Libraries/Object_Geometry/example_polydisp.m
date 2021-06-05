%% Notes
%
%   Author: Koen Muller, 08 May 2020
%
% example script for fitting a polynomial displacement field to trajectorie
%   data
%
% this function allows testing 2D and 3D displacement field fit, run and 
%   see code for the examples, also look inside polydisp and dispeval
%
% linear deformation field allows generating exact trajectories, and thus
%   an exact fit
%
% it is interesting to see how a time varying linear deformation field can
%   capture the exact motion of complex motions with in it
%
% nonlinear deformation field requires time integration to generate correct
%   trajectories
%
% multiple type of motion can be selected, see controls

%% Start
close all
clear all
clc

%% Controls

% trajectories
I=100; % number of trajectories
N=100; % number frames [even number]
tord=2; % 0 | 1 | 2 trajectory order
tini='full'; % random or fullspan start/end [in]complete data
tdyn='rotat'; % displac[e/ing] | stretch[ing] | rotat[e/ing]

% order deformation field
xord=1; % 1 | 2 

% order fit displacementfield
ord=[1 1 1 2]; % [1 1 1 3] is a exact fit for 2-ord trajectory in linear deformation field

% samping indexed grid
Ngrid=10;

%% Define some synthetic trajectories

% trajectory dynamics
d0=double(tord>=0)*(2*rand(3,1)-1);
d1=double(tord>=1)*(2*rand(3,1)-1);
d2=double(tord>=2)*(2*rand(3,1)-1);

% define trajectory data
tdat=zeros(5,0);
for i=1:I
    
    % initial condition
    s0=(2*rand(3,1)-1);
    
    % initiate frames
    switch tini
        case 'random'
            nspan=sort([-randi(N/2,1,1) randi(N/2,1,1)]); % random start and end
        case 'full'
            nspan=[-N/2 N/2]; % full span time series
    end
    
    % loop span frames
    for n=nspan(1):nspan(2)
        
        % generate a trajectory in a deformation field accordingly
        switch tdyn
            case 'displac'
                x=s0.*d0+...
                    d1*((n)/N)+...
                    d2*((n)/N)^2; % expanding accelarating trajectories
            case 'stretch'
                x=s0.*(d0+...
                    d1*((n)/N)+...
                    d2*((n)/N)^2); % expanding accelarating trajectories
            case 'rotat'
                x=s0.*d0+...
                    [0 1 0 ; -1 0 0 ; 0 0 1]*s0.*(d1*((n)/N)+...
                    d2*((n)/N)^2); % expanding accelarating trajectories
        end
        
        % displacement field order
        x=double(xord>=1)*x+double(xord>=2)*x.^2;
        
        % append track
        tdat=cat(2,tdat,[i;n;x]);
        
    end
    
end

% message
disp(['generated ',num2str(I),' ',...
    tdyn,'ing ',...
    num2str(tord),'-order trajectories spanning ',...
    num2str(N), ' frames with a ',tini,'span track'])

%% Plot 2D trajectories

% indexing
figure
scatter(tdat(3,:),tdat(4,:),[],tdat(1,:),'.')
colormap lines
axis tight
axis equal
xlabel('x')
ylabel('y')
title('2D Trajectories Indexing')
drawnow

% timeseries
figure
scatter(tdat(3,:),tdat(4,:),[],tdat(2,:),'.')

% layout
colormap parula
col=colorbar;
col.Label.String = 'frame number';
col.FontSize=11;
axis tight
axis equal
xlabel('x')
ylabel('y')
title('2D Trajectories Timeseries')
drawnow

%% fit 2D displacement field

% use partial data to fit 
[coef,fdat]=polydisp(tdat(1:4,:),ord([1,2,4]),4);

%% define indexed mesh on intial points
% read span grid initial frame from generated trajectories
X=linspace(min(fdat(2,:)),max(fdat(2,:)),Ngrid);
Y=linspace(min(fdat(3,:)),max(fdat(3,:)),Ngrid);
Z=linspace(min(fdat(4,:)),max(fdat(4,:)),Ngrid);

% mesh regular grid
[X,Y,Z]=ndgrid(X,Y,Z);

% indexed grid data
xdat=[X(:)'
    Y(:)'
    Z(:)'];

%% plot 2D fitted data

% evaluate displacemen
dX= dispeval( coef, fdat(1:4,:) ,[0 0 0] );

% compute positions
X=fdat(2:4,:)+dX;
 
% figure
figure
plot(tdat(3,:),tdat(4,:),'.','Color',[0.75 0.75 0.75])
hold on
scatter(X(1,:),X(2,:),[],log(fdat(5,:)+1),'.')
hold off

% layout
colormap parula
col=colorbar;
col.Label.String = 'residual';
col.FontSize=11;
axis tight
axis equal
xlabel('x')
ylabel('y')
title('2D Fit Displacement Field to Trajectories')
drawnow

%% Animate 2D displacement field 

% figure
figure

% loop frames
for n=-N/2:N/2
    
    % original trajectory
    plot(tdat(3,tdat(2,:)<=n),tdat(4,tdat(2,:)<=n),'.','Color',[0.75 0.75 0.75])
    hold on
    
    % evaluate displacement fit trajectory
    dX= dispeval( coef, fdat(1:4,tdat(2,:)<=n) ,[0 0 0] );
    
    % compute fit positions
    X=fdat(2:4,tdat(2,:)<=n)+dX;
    
    % plot fit
    scatter(X(1,:),X(2,:),[],log(fdat(5,tdat(2,:)<=n)+1),'.')
    
    % evaluate displacement end point trajectory
    dX= dispeval( coef, fdat(1:4,tdat(2,:)==n) ,[0 0 0] );
    
    % compute position end point
    X=fdat(2:4,tdat(2,:)==n)+dX;
    
    % velocity vector end point fit
    plot(X(1,:),X(2,:),'.k')
    
    % evaluate displacement fit
    U= N/Ngrid*dispeval( coef, fdat(1:4,tdat(2,:)==n) ,[0 0 1] ); % scale visual
    
    % velocity vector fit
    quiver(X(1,:),X(2,:),U(1,:),U(2,:),0,'k')
    
    % evaluate displacement field
    dX= dispeval( coef, [n*ones(1,Ngrid^3); xdat(1:2,:) ] ,[0 0 0] );
    
    % evaluate displaced grid position 
    X=xdat+dX;
    
    % evaluate velocity field 
    U= N/Ngrid*dispeval( coef, [n*ones(1,Ngrid^3); xdat(1:2,:) ] ,[0 0 1] ); % scale visual
    
    % plot displacement vector
    plot(X(1,:),X(2,:),'.','Color',[0.9290 0.6940 0.1250])
    quiver(X(1,:),X(2,:),U(1,:),U(2,:),0,'Color',[0.9290 0.6940 0.1250])
    
    % layout
    hold off
    colormap parula
    col=colorbar;
    col.Label.String = 'trajectory fit residual';
    axis tight
    axis equal
    xlabel('x')
    ylabel('y')
    title('Trajectories in Fitted Displacement field')
    drawnow
    
end

%% plot 3D trajectories
if 0
% indexing
figure
scatter3(tdat(3,:),tdat(4,:),tdat(5,:),[],tdat(1,:),'.')
colormap lines
axis tight
axis equal
xlabel('x')
ylabel('y')
title('2D Trajectories Indexing')
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

%% fit 3D displacement field

% use partial data to fit 
[coef, fdat]=polydisp(tdat(1:5,:),ord,4);

%% define indexed mesh on intial points
% read span grid initial frame from generated trajectories
X=linspace(min(fdat(2,:)),max(fdat(2,:)),Ngrid);
Y=linspace(min(fdat(3,:)),max(fdat(3,:)),Ngrid);
Z=linspace(min(fdat(4,:)),max(fdat(4,:)),Ngrid);

% mesh regular grid
[X,Y,Z]=ndgrid(X,Y,Z);

% indexed grid data
xdat=[X(:)'
    Y(:)'
    Z(:)'];

%% plot fitted data

% evaluate displacemen
dX= dispeval( coef, fdat(1:4,:) ,[0 0 0 0] );

% compute positions
X=fdat(2:4,:)+dX;
 
% figure
figure
plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.','Color',[0.75 0.75 0.75])
hold on
scatter3(X(1,:),X(2,:),X(3,:),[],log(fdat(5,:)+1),'.')
hold off

% layout
colormap parula
col=colorbar;
col.Label.String = 'residual';
col.FontSize=11;
axis tight
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Fit Displacement Field to Trajectories')
drawnow

%% Animate 3D displacement field 

% figure
figure

% loop frames
for n=-N/2:N/2
    
    % original trajectory
    plot3(tdat(3,tdat(2,:)<=n),tdat(4,tdat(2,:)<=n),tdat(5,tdat(2,:)<=n),'.','Color',[0.75 0.75 0.75])
    hold on
    
    % evaluate displacement fit trajectory
    dX= dispeval( coef, fdat(1:4,tdat(2,:)<=n) ,[0 0 0 0] );
    
    % compute fit positions
    X=fdat(2:4,tdat(2,:)<=n)+dX;
    
    % plot fit
    scatter3(X(1,:),X(2,:),X(3,:),[],log(fdat(5,tdat(2,:)<=n)+1),'.')
    
    % evaluate displacement end point trajectory
    dX= dispeval( coef, fdat(1:4,tdat(2,:)==n) ,[0 0 0 0] );
    
    % compute position end point
    X=fdat(2:4,tdat(2,:)==n)+dX;
    
    % velocity vector end point fit
    plot3(X(1,:),X(2,:),X(3,:),'.k')
    
    % evaluate displacement fit
    U= N/Ngrid*dispeval( coef, fdat(1:4,tdat(2,:)==n) ,[0 0 0 1] ); % scale visual
    
    % velocity vector fit
    quiver3(X(1,:),X(2,:),X(3,:),U(1,:),U(2,:),U(3,:),0,'k')
    
    % evaluate displacement field
    dX= dispeval( coef, [n*ones(1,Ngrid^3); xdat ] ,[0 0 0 0] );
    
    % evaluate displaced grid position 
    X=xdat+dX;
    
    % evaluate velocity field 
    U= N/Ngrid*dispeval( coef, [n*ones(1,Ngrid^3); xdat ] ,[0 0 0 1] ); % scale visual
    
    % plot displacement vector
    plot3(X(1,:),X(2,:),X(3,:),'.','Color',[0.9290 0.6940 0.1250])
    quiver3(X(1,:),X(2,:),X(3,:),U(1,:),U(2,:),U(3,:),0,'Color',[0.9290 0.6940 0.1250])
    
    % layout
    hold off
    colormap parula
    col=colorbar;
    col.Label.String = 'trajectory fit residual';
    axis tight
    axis equal
    xlabel('x')
    ylabel('y')
    ylabel('z')
    title('3D Trajectories in Fitted Displacement field')
    drawnow
    
end
end