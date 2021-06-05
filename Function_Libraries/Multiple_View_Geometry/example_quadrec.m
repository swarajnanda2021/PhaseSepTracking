%% notes

% author koen muller

% try-out for quadric reconstruction with constrained translation

%% start
close all
clear all
clc

%% Additional Paths
addpath(genpath([cd '\..\..\Function_Libraries\']))

%%%%%%%%%%%%%%%%
%%% controls %%%
%%%%%%%%%%%%%%%%

%% number quadrics
N=1000;

%% define cameras
% define projection matrix
Pmat=cell(1,4);

% define matrices omit the camera matrix
Pmat{1}=[ eye(3) [-3 0 0]' ]; % (K)[R t]
Pmat{2}=[ eye(3) [ 3 0 0]' ];
Pmat{3}=[ eye(3) [-3 1 0]' ];
Pmat{4}=[ eye(3) [ 3 1 0]' ];
% Pmat{5}=[ eye(3) [ 0 0 0 ]' ];
% Pmat{6}=[ eye(3) [ 10 10 5]' ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define single quadric to reconstruct %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define quadric
% define shape
Lam=inv(diag([3 1 1/2].^2)); % semi axis sort(rand(1,3),'descend')
lam=-1; % ellipsoid
Q0=[Lam zeros(3,1) ; zeros(1,3) lam]; % principle axis

% apply some coordinate transformation
R=eula2rotm(pi*1*rand(1,3),'ZYX');%pi*[0 0 0]
t=[0 0 10]';
H=[R t ; zeros(1,3) 1];

% transform quadric and compute its inverse
q0=qmat2qvec(Q0);
q=quadtrans(q0,H,'transmat');%inv(H)'*Q0*inv(H);

disp([['Original      Q = ' ;'                  ';'                  ' ;'                  ' ],num2str(qvec2qmat(q))])

%% compute camera conic projections
% initiate
cdat=cell(size(Pmat)); % conics
xdat=cell(size(Pmat)); % midpoints

% loop views
for i=1:length(cdat)
    
    % conic vector
    cdat{i}=projquad(q,Pmat{i});%*(0.5-rand(1)); % randly scale
    
    % image midpoint
    xdat{i}=cvec2eshp(cdat{i});
    
end

%% plot original quadric
% make figure
figure
hold on;

% implicit plot
fimplicit3(@(x,y,z)(q(1).*(x).^2+...
    q(2).*(x).*(y)+...
    q(3).*(x).*(z)+...
    q(4).*(y).^2+...
    q(5).*(y).*(z)+...
    q(6).*(z).^2+...
    q(7).*(x)+...
    q(8).*(y)+...
    q(9).*(z)+...
    q(10)),[t(1)-5 t(1)+5 t(2)-5 t(2)+5 t(3)-5 t(3)+5],'FaceAlpha',0.5,'EdgeColor','none')

% explicit
X=qvec2pnts(q);
plot3(X(1,:),X(2,:),X(3,:),'b.')

% midpoint
plot3(t(1,:),t(2,:),t(3,:),'c*')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
xlim([-5 5])
ylim([-5 5])
zlim([0 20])
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Original Quadric','interpreter','latex')
drawnow

%% plot projections
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic data
    c=cdat{i};
    x=xdat{i};
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % implicit plot
    fimplicit(@(x,y)(c(1).*(x).^2+...
        c(2).*(x).*(y)+...
        c(3).*(y).^2+...
        c(4).*(x)+...
        c(5).*(y)+...
        c(6)),[-1 1 -1 1],'-b')%0 size(sub,1) 0 size(sub,2)
    
    % explicit plot
    X=cvec2pnts(c);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title(['Ellipse Camera ',num2str(i)],'interpreter','latex')
    
end

% figure
drawnow

%% object triangulation
tr=objtriang(xdat,Pmat);

%% Align reprojection camera conic projections
cadat=cell(size(cdat));

% loop views
for i=1:length(cdat)
    
    % perturbation
    d=projcoord(tr,Pmat{i})-xdat{i};
    
    % conic vector
    cadat{i}=contrans(cdat{i},d,'displace');%*(0.5-rand(1)); % randly scale
    
end

%% solve reconstruction
qr=quadrec(cadat,Pmat);

disp([['Reconstructed Q = ' ;'                  ';'                  ' ;'                  ' ],num2str(qvec2qmat(qr))])

%% plot the reconstructed quadric
% make figure
figure
hold on;

% implicit plot
if isnan(sum(qr(:)))==0
    fimplicit3(@(x,y,z)(qr(1).*(x).^2+...
        qr(2).*(x).*(y)+...
        qr(3).*(x).*(z)+...
        qr(4).*(y).^2+...
        qr(5).*(y).*(z)+...
        qr(6).*(z).^2+...
        qr(7).*(x)+...
        qr(8).*(y)+...
        qr(9).*(z)+...
        qr(10)),[tr(1)-5 tr(1)+5 tr(2)-5 tr(2)+5 tr(3)-5 tr(3)+5],'FaceAlpha',0.5,'EdgeColor','none')
end
% explicit
X=qvec2pnts(qr);
plot3(X(1,:),X(2,:),X(3,:),'r.')

% midpoint
plot3(tr(1,:),tr(2,:),tr(3,:),'m*')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
xlim([-5 5])
ylim([-5 5])
zlim([0 20])
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Reconstructed Quadric','interpreter','latex')
drawnow

%% plot projection
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic data
    c=cdat{i};
    x=xdat{i};
    
    % get projection data
    cp=projquad(qr,Pmat{i});
    xp=projcoord(qvec2eshp(qr),Pmat{i});
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % implicit plot
    fimplicit(@(x,y)(c(1).*(x).^2+...
        c(2).*(x).*(y)+...
        c(3).*(y).^2+...
        c(4).*(x)+...
        c(5).*(y)+...
        c(6)),[-1 1 -1 1],'-b')%0 size(sub,1) 0 size(sub,2)
    
    % implicit plot
    fimplicit(@(x,y)(cp(1).*(x).^2+...
        cp(2).*(x).*(y)+...
        cp(3).*(y).^2+...
        cp(4).*(x)+...
        cp(5).*(y)+...
        cp(6)),[-1 1 -1 1],'-r')%0 size(sub,1) 0 size(sub,2)
    
    % explicit plot
    X=cvec2pnts(c);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % explicit plot
    X=cvec2pnts(cp);
    plot(X(1,:),X(2,:),'r.')
    
    % explicit plot
    plot(xp(1,:),xp(2,:),'mo')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title(['Ellipse Camera ',num2str(i)],'interpreter','latex')
    
end

% figure
drawnow

%% solve reconstruction
qr=quadrec(cadat,Pmat,tr,'free');

disp([['Reconstructed Q = ' ;'                  ';'                  ' ;'                  ' ],num2str(qvec2qmat(qr))])

%% plot reconstructed quadric with constrained triangulation vector
% make figure
figure
hold on;

% implicit plot
if isnan(sum(qr(:)))==0
    fimplicit3(@(x,y,z)(qr(1).*(x).^2+...
        qr(2).*(x).*(y)+...
        qr(3).*(x).*(z)+...
        qr(4).*(y).^2+...
        qr(5).*(y).*(z)+...
        qr(6).*(z).^2+...
        qr(7).*(x)+...
        qr(8).*(y)+...
        qr(9).*(z)+...
        qr(10)),[tr(1)-5 tr(1)+5 tr(2)-5 tr(2)+5 tr(3)-5 tr(3)+5],'FaceAlpha',0.5,'EdgeColor','none')
end

% explicit
X=qvec2pnts(qr);
plot3(X(1,:),X(2,:),X(3,:),'r.')

% midpoint
plot3(tr(1,:),tr(2,:),tr(3,:),'m*')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
xlim([-5 5])
ylim([-5 5])
zlim([0 20])
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Spherical Quadric with Prescribed Triangulation','interpreter','latex')
drawnow

%% plot projection
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic data
    c=cdat{i};
    x=xdat{i};
    
    % get projection data
    cp=projquad(qr,Pmat{i});
    xp=projcoord(qvec2eshp(qr),Pmat{i});
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % implicit plot
    fimplicit(@(x,y)(c(1).*(x).^2+...
        c(2).*(x).*(y)+...
        c(3).*(y).^2+...
        c(4).*(x)+...
        c(5).*(y)+...
        c(6)),[-1 1 -1 1],'-b')%0 size(sub,1) 0 size(sub,2)
    
    % implicit plot
    fimplicit(@(x,y)(cp(1).*(x).^2+...
        cp(2).*(x).*(y)+...
        cp(3).*(y).^2+...
        cp(4).*(x)+...
        cp(5).*(y)+...
        cp(6)),[-1 1 -1 1],'-r')%0 size(sub,1) 0 size(sub,2)
    
    % explicit plot
    X=cvec2pnts(c);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % explicit plot
    X=cvec2pnts(cp);
    plot(X(1,:),X(2,:),'r.')
    
    % explicit plot
    plot(xp(1,:),xp(2,:),'mo')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title(['Ellipse Camera ',num2str(i)],'interpreter','latex')
    
end

% figure
drawnow

%% solve reconstruction
qr=quadrec(cadat,Pmat,tr,'diag');

disp([['Reconstructed Q = ' ;'                  ';'                  ' ;'                  ' ],num2str(qvec2qmat(qr))])

%% plot reconstructed quadric with constrained triangulation vector
% make figure
figure
hold on;

% implicit plot
if isnan(sum(qr(:)))==0
    fimplicit3(@(x,y,z)(qr(1).*(x).^2+...
        qr(2).*(x).*(y)+...
        qr(3).*(x).*(z)+...
        qr(4).*(y).^2+...
        qr(5).*(y).*(z)+...
        qr(6).*(z).^2+...
        qr(7).*(x)+...
        qr(8).*(y)+...
        qr(9).*(z)+...
        qr(10)),[tr(1)-5 tr(1)+5 tr(2)-5 tr(2)+5 tr(3)-5 tr(3)+5],'FaceAlpha',0.5,'EdgeColor','none')
end

% explicit
X=qvec2pnts(qr);
plot3(X(1,:),X(2,:),X(3,:),'r.')

% midpoint
plot3(tr(1,:),tr(2,:),tr(3,:),'m*')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
xlim([-5 5])
ylim([-5 5])
zlim([0 20])
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Diagonal Quadric with Prescribed Triangulation','interpreter','latex')
drawnow

%% plot projection
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic data
    c=cdat{i};
    x=xdat{i};
    
    % get projection data
    cp=projquad(qr,Pmat{i});
    xp=projcoord(qvec2eshp(qr),Pmat{i});
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % implicit plot
    fimplicit(@(x,y)(c(1).*(x).^2+...
        c(2).*(x).*(y)+...
        c(3).*(y).^2+...
        c(4).*(x)+...
        c(5).*(y)+...
        c(6)),[-1 1 -1 1],'-b')%0 size(sub,1) 0 size(sub,2)
    
    % implicit plot
    fimplicit(@(x,y)(cp(1).*(x).^2+...
        cp(2).*(x).*(y)+...
        cp(3).*(y).^2+...
        cp(4).*(x)+...
        cp(5).*(y)+...
        cp(6)),[-1 1 -1 1],'-r')%0 size(sub,1) 0 size(sub,2)
    
    % explicit plot
    X=cvec2pnts(c);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % explicit plot
    X=cvec2pnts(cp);
    plot(X(1,:),X(2,:),'r.')
    
    % explicit plot
    plot(xp(1,:),xp(2,:),'mo')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title(['Ellipse Camera ',num2str(i)],'interpreter','latex')
    
end

% figure
drawnow

%% solve reconstruction
qr=quadrec(cadat,Pmat,tr,'spher');

disp([['Reconstructed Q = ' ;'                  ';'                  ' ;'                  ' ],num2str(qvec2qmat(qr))])

%% plot reconstructed quadric with constrained triangulation vector
% make figure
figure
hold on;

% implicit plot
if isnan(sum(qr(:)))==0
    fimplicit3(@(x,y,z)(qr(1).*(x).^2+...
        qr(2).*(x).*(y)+...
        qr(3).*(x).*(z)+...
        qr(4).*(y).^2+...
        qr(5).*(y).*(z)+...
        qr(6).*(z).^2+...
        qr(7).*(x)+...
        qr(8).*(y)+...
        qr(9).*(z)+...
        qr(10)),[tr(1)-5 tr(1)+5 tr(2)-5 tr(2)+5 tr(3)-5 tr(3)+5],'FaceAlpha',0.5,'EdgeColor','none')
end

% explicit
X=qvec2pnts(qr);
plot3(X(1,:),X(2,:),X(3,:),'r.')

% midpoint
plot3(tr(1,:),tr(2,:),tr(3,:),'m*')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
xlim([-5 5])
ylim([-5 5])
zlim([0 20])
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Spherical Quadric with Prescribed Triangulation','interpreter','latex')
drawnow

%% plot projection
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic data
    c=cdat{i};
    x=xdat{i};
    
    % get projection data
    cp=projquad(qr,Pmat{i});
    xp=projcoord(qvec2eshp(qr),Pmat{i});
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % implicit plot
    fimplicit(@(x,y)(c(1).*(x).^2+...
        c(2).*(x).*(y)+...
        c(3).*(y).^2+...
        c(4).*(x)+...
        c(5).*(y)+...
        c(6)),[-1 1 -1 1],'-b')%0 size(sub,1) 0 size(sub,2)
    
    % implicit plot
    fimplicit(@(x,y)(cp(1).*(x).^2+...
        cp(2).*(x).*(y)+...
        cp(3).*(y).^2+...
        cp(4).*(x)+...
        cp(5).*(y)+...
        cp(6)),[-1 1 -1 1],'-r')%0 size(sub,1) 0 size(sub,2)
    
    % explicit plot
    X=cvec2pnts(c);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % explicit plot
    X=cvec2pnts(cp);
    plot(X(1,:),X(2,:),'r.')
    
    % explicit plot
    plot(xp(1,:),xp(2,:),'mo')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title(['Ellipse Camera ',num2str(i)],'interpreter','latex')
    
end

% figure
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define multiple quadrics to reconstruct %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define quadrics

% initiate
q=zeros(10,N);
t=zeros(3,N);

% generate
for n=1:N
    
    % define shape
    Lam=inv(diag((1*rand(1,3)).^2)); % semi axis sort(rand(1,3),'descend')
    lam=-1; % ellipsoid
    Q0=[Lam zeros(3,1) ; zeros(1,3) lam]; % principle axis
    
    % apply some coordinate transformation
    R=eula2rotm(pi*2*rand(1,3),'ZYX');%pi*[0 0 0]
    t(:,n)=[0 0 15]'+10*(1/2-rand(3,1));
    H=[R t(:,n) ; zeros(1,3) 1];
    
    % transform quadric and compute its inverse
    Q=inv(H)'*Q0*inv(H);
    q(:,n)=qmat2qvec(Q);
    
end

%% compute camera conic projections
% initiate
cdat=cell(size(Pmat)); % conics
xdat=cell(size(Pmat)); % midpoints

% loop views
for i=1:length(cdat)
    
    % initiate
    cdat{i}=zeros(6,N);
    xdat{i}=zeros(2,N);
    
    % conic vector
    cdat{i}=projquad(q,Pmat{i});
    
    % image midpoint
    xdat{i}=cvec2eshp(cdat{i});
    
end

%% plot original quadric
% make figure
figure
hold on;

% explicit
X=reshape(qvec2pnts(q),3,[]);
plot3(X(1,:),X(2,:),X(3,:),'b.')

% midpoint
plot3(t(1,:),t(2,:),t(3,:),'c*')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
axis tight
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Original Quadrics','interpreter','latex')
drawnow

%% plot projection
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic
    c=cdat{i};
    x=xdat{i};
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % explicit plot
    X=reshape(cvec2pnts(c),2,[]);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title('Original Ellipse Projections','interpreter','latex')
    
end

%% object triangulation
tr=objtriang(xdat,Pmat);

%% Align reprojection camera conic projections
cadat=cell(size(cdat));

% loop views
for i=1:length(cdat)
    
    % perturbation
    d=projcoord(tr,Pmat{i})-xdat{i};
    
    % conic vector
    cadat{i}=contrans(cdat{i},d,'displace');%*(0.5-rand(1)); % randly scale
    
end

%% solve reconstruction
tic
qr=quadrec(cadat,Pmat);
toc

%% plot the reconstructed quadric
% make figure
figure
hold on;

% explicit
X=reshape(qvec2pnts(q),3,[]);
plot3(X(1,:),X(2,:),X(3,:),'b.')

% midpoint
plot3(t(1,:),t(2,:),t(3,:),'c*')

% explicit
X=reshape(qvec2pnts(qr),3,[]);
plot3(X(1,:),X(2,:),X(3,:),'r.')

% midpoint
plot3(tr(1,:),tr(2,:),tr(3,:),'mo')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
axis tight
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Reconstructed Quadrics','interpreter','latex')
drawnow

%% plot projection
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic data
    c=cdat{i};
    x=xdat{i};
    
    % get projection data
    cp=projquad(qr,Pmat{i});
    xp=projcoord(qvec2eshp(qr),Pmat{i});
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % implicit plot
    fimplicit(@(x,y)(c(1).*(x).^2+...
        c(2).*(x).*(y)+...
        c(3).*(y).^2+...
        c(4).*(x)+...
        c(5).*(y)+...
        c(6)),[-1 1 -1 1],'-b')%0 size(sub,1) 0 size(sub,2)
    
    % implicit plot
    fimplicit(@(x,y)(cp(1).*(x).^2+...
        cp(2).*(x).*(y)+...
        cp(3).*(y).^2+...
        cp(4).*(x)+...
        cp(5).*(y)+...
        cp(6)),[-1 1 -1 1],'-r')%0 size(sub,1) 0 size(sub,2)
    
    % explicit plot
    X=cvec2pnts(c);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % explicit plot
    X=cvec2pnts(cp);
    plot(X(1,:),X(2,:),'r.')
    
    % explicit plot
    plot(xp(1,:),xp(2,:),'mo')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title(['Ellipse Camera ',num2str(i)],'interpreter','latex')
    
end

% figure
drawnow

%% solve reconstruction
tic
qr=quadrec(cadat,Pmat,tr,'spher');
toc

%% plot the reconstructed quadric
% make figure
figure
hold on;

% explicit
X=reshape(qvec2pnts(q),3,[]);
plot3(X(1,:),X(2,:),X(3,:),'b.')

% midpoint
plot3(t(1,:),t(2,:),t(3,:),'c*')

% explicit
X=reshape(qvec2pnts(qr),3,[]);
plot3(X(1,:),X(2,:),X(3,:),'r.')

% midpoint
plot3(tr(1,:),tr(2,:),tr(3,:),'mo')

% layout
hold off
axis equal
grid minor
box on
camproj('perspective')
view(3)
axis tight
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
zlabel('z','interpreter','latex')
title('Contrained Reconstructed Quadrics','interpreter','latex')
drawnow

%% plot projection
% make figure
figure

% loop views
for i=1:length(cdat)
    
    % get conic data
    c=cdat{i};
    x=xdat{i};
    
    % get projection data
    cp=projquad(qr,Pmat{i});
    xp=projcoord(qvec2eshp(qr),Pmat{i});
    
    % figure
    subplot(floor(length(Pmat)/2)+1,ceil(length(Pmat)/2),i)
    hold on;
    
    % implicit plot
    fimplicit(@(x,y)(c(1).*(x).^2+...
        c(2).*(x).*(y)+...
        c(3).*(y).^2+...
        c(4).*(x)+...
        c(5).*(y)+...
        c(6)),[-1 1 -1 1],'-b')%0 size(sub,1) 0 size(sub,2)
    
    % implicit plot
    fimplicit(@(x,y)(cp(1).*(x).^2+...
        cp(2).*(x).*(y)+...
        cp(3).*(y).^2+...
        cp(4).*(x)+...
        cp(5).*(y)+...
        cp(6)),[-1 1 -1 1],'-r')%0 size(sub,1) 0 size(sub,2)
    
    % explicit plot
    X=cvec2pnts(c);
    plot(X(1,:),X(2,:),'b.')
    
    % explicit plot
    plot(x(1,:),x(2,:),'c*')
    
    % explicit plot
    X=cvec2pnts(cp);
    plot(X(1,:),X(2,:),'r.')
    
    % explicit plot
    plot(xp(1,:),xp(2,:),'mo')
    
    % layout
    hold off
    axis equal
    grid minor
    box on
    xlabel('x','interpreter','latex')
    ylabel('y','interpreter','latex')
    title(['Ellipse Camera ',num2str(i)],'interpreter','latex')
    
end

% figure
drawnow

%% Sandbox

