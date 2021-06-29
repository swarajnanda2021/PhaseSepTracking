function postproc_CavityCondAvg




%% get globals
global folder date rec prop mset post%plotting

%% Load data
% load general data
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
% Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Time.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
load([folder date rec vsl 'Postproc3DTracking' vsl 'TrackQualityMetric.mat'],'t_spurious');
% ( identification index | camera index | frame number | conic vector )
Cpln=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Cpln.dat']);
% ( track index | indentification index | camera index | frame number | camera conic | camera velocity | object position | object velocity | res )
Plink=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Plink.dat']);

% Load VH Data

% Data tabulation is as follows
% Visual_hull_reduced(1,i,n) = n;
% Visual_hull_reduced(2,i,n) = Centx(i);
% Visual_hull_reduced(3,i,n) = Centy(i);
% Visual_hull_reduced(4,i,n) = Orient;
% Visual_hull_reduced(5,i,n) = MajAx;
% Visual_hull_reduced(6,i,n) = MinAx;
% Visual_hull_reduced(7,i,n) = Area;
% Visual_hull_reduced(8,i,n) = Z_vox(1,1,i);
VHData = load([folder date rec vsl 'Postproc3DReco' vsl 'VHReduced.mat']);

% Find mean cavity radius
r_c_m = nanmean(0.5*(4*VHData.Visual_hull_reduced(7,:,:)./pi).^0.5,[2 3]);
r_c_std = nanstd(0.5*(4*VHData.Visual_hull_reduced(7,:,:)./pi).^0.5,0,[2 3]);


%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'CavityReference'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'CavityReference'])
end

%%

Rdata = zeros(1,0);
Timedata = zeros(1,0);
Zdata = zeros(1,0);
Urdata = zeros(1,0);
Uthdata = zeros(1,0);
Uzdata = zeros(1,0);    
Thetadata = zeros(1,0);    
Focalitydata = zeros(1,0);     
PosXdata = zeros(1,0);
PosYdata = zeros(1,0);


% loop frames
for n=post.tproc(1):post.tproc(2)
    disp(['Starting tstep ' num2str(n)])
    
    
    %% For visual hull processing
    % Load visual hull and estimate axis
    load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(n) '.mat'],'Centroid_cav')
    % Reshape 3D matrix to 2D with transverse dimensions squashed    

        
    
    
    %% For velocity processing
    % Time set
    N=(ismember(Index(2,:),n));  
    outl = ~(ismember(Index(1,:),unique(t_spurious)));
    % Gather position and velocities
    pos = Position(:,N&outl);
    vel = Velocity(:,N&outl);
    tindx = Index(1,N&outl);
    timeindx = Index(2,N&outl);
%     dummie = size(pos)
    
    % determine radial position
    centre_x = interp1(Centroid_cav(3,:),Centroid_cav(1,:),pos(3,:));
    centre_y = interp1(Centroid_cav(3,:),Centroid_cav(2,:),pos(3,:));
    
    r_pos = ((pos(1,:)-centre_x).^2 + (pos(2,:)-centre_y).^2).^0.5;
    theta_pos = atan2((pos(2,:)-centre_y),(pos(1,:)-centre_x));
%     u_theta_pos = -(vel(1,:).*sin(theta_pos)./r_pos) + (vel(2,:).*cos(theta_pos)./r_pos); % in radians/sec
    u_theta_pos = -(vel(1,:).*sin(theta_pos)) + (vel(2,:).*cos(theta_pos));
    u_r_pos = (vel(1,:).*cos(theta_pos)) + (vel(2,:).*sin(theta_pos));
    u_z_pos = vel(3,:);
    
    %% Extract camera focality
    % Extract relevant information from Plink
    Plink_sub = Plink(1:2,ismember(Plink(4,:),n));
    Focality = zeros(size(tindx));
    for kk=1:length(tindx) % loop over track indices
        % Find identification indices
        Focality(kk) = length(Plink_sub(2,ismember(Plink_sub(1,:),tindx(kk))   ));
    end
    %% Gather data
    Rdata = cat(2,Rdata,r_pos);
    Zdata = cat(2,Zdata,pos(3,:));
    Thetadata = cat(2,Thetadata,theta_pos);
    Urdata = cat(2,Urdata,u_r_pos);
    Uthdata = cat(2,Uthdata,u_theta_pos);
    Uzdata = cat(2,Uzdata,u_z_pos);
    Focalitydata = cat(2,Focalitydata,Focality);
    Timedata = cat(2,Timedata,timeindx);
    
    PosXdata = cat(2,PosXdata,(pos(1,:)-centre_y));
    PosYdata = cat(2,PosYdata,(pos(2,:)-centre_x));
    
    if 0 %plotcheck
        figure(20)
        hold all
        scatter(PosXdata,PosYdata)
        drawnow
%         pause
        
    end
    
    
    %% save data at end of each timestep
    if 0
        disp(['Storing tstep ' num2str(n) ' completed!'])
        save([folder date rec vsl 'Postproc3DTracking' vsl 'CavityReference' vsl 'Tstep_' num2str(n) '.mat'],'pos_new','vec_accum_mat')
    end
    
    
end

%% Begin binning in cartesian coordinates

q_x = -25:1:25;
q_y=q_x;
tempqx = 0.5*(q_x(2)-q_x(1));

[~,bin_posx_edges ,idx_posx]  = histcounts(PosYdata,q_x); % distances in units of object space (here, mm)
[~,bin_posy_edges ,idx_posy]  = histcounts(PosXdata,q_y);


u_theta_cart = zeros(length(unique(idx_posx)),length(unique(idx_posy)));
u_theta_cart_std = u_theta_cart;
u_theta_cart_cnt = u_theta_cart;

u_z_cart = u_theta_cart;
u_z_cart_std = u_theta_cart;
u_z_cart_cnt = u_theta_cart;

u_r_cart = u_theta_cart;
u_r_cart_std = u_theta_cart;
u_r_cart_cnt = u_theta_cart;


nfoc_3 = u_theta_cart;
nfoc_4 = u_theta_cart;
nfoc_5 = u_theta_cart;

dum1 = unique(idx_posx);
dum2 = unique(idx_posy);
for i=1:length(dum1)
    for j=1:length(dum2)
        
        % Bin velocity components
        
        u_theta_cart(i,j) = nanmean(Uthdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))));
        u_theta_cart_std(i,j) = nanstd(Uthdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))));
        u_theta_cart_cnt(i,j) = length(find(~isnan(Uthdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))))));
        
        
        u_r_cart(i,j) = nanmean(Urdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))));
        u_r_cart_std(i,j) = nanstd(Urdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))));
        u_r_cart_cnt(i,j) = length(find(~isnan(Urdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))))));
        
        u_z_cart(i,j) = nanmean(Uzdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))));
        u_z_cart_std(i,j) = nanstd(Uzdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))));
        u_z_cart_cnt(i,j) = length(find(~isnan(Uzdata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j))))));
        
        % Bin particle focality
        
        Focality_sub = Focalitydata(ismember(idx_posx,dum1(i)) & ismember(idx_posy,dum2(j)));
        
        nfoc_3(i,j) = length(Focality_sub(ismember(Focality_sub,3)));
        nfoc_4(i,j) = length(Focality_sub(ismember(Focality_sub,4)));
        nfoc_5(i,j) = length(Focality_sub(ismember(Focality_sub,5)));
        
        
        
        
    end
end

%% Plotting velocities in cartesian coordinates

% Image 1: binning focality

figure(1)
subplot(2,3,1)
hold all
surf(q_x(2:end-1)+tempqx,q_y(2:end-1)+tempqx,nfoc_3((2:end-1),(2:end-1)),'edgecolor','none')
xlim([-20,20])
ylim([-20,20])
view(2)

axis equal
axis tight
colorbar;
% colormap(brewermap(10,'Spectral'))
box on
set(gca,'linewidth',1,'fontsize',15)



subplot(2,3,2)
surf(q_x(2:end-1)+tempqx,q_y(2:end-1)+tempqx,nfoc_4((2:end-1),(2:end-1)),'edgecolor','none')
xlim([-20,20])
ylim([-20,20])
view(2)
axis equal
axis tight
colorbar;
% colormap(brewermap(10,'Spectral'))
box on
set(gca,'linewidth',1,'fontsize',15)

subplot(2,3,3)
surf(q_x(2:end-1)+tempqx,q_y(2:end-1)+tempqx,nfoc_5((2:end-1),(2:end-1)),'edgecolor','none')
xlim([-20,20])
ylim([-20,20])
view(2)
axis equal
axis tight
colorbar;
% colormap(brewermap(10,'Spectral'))
box on
set(gca,'linewidth',1,'fontsize',15)






% Image 2: binning velocities

% figure(2)
subplot(2,3,4)
surf(q_x(2:end-1)+tempqx,q_y(2:end-1)+tempqx,u_r_cart((2:end-1),(2:end-1)),'edgecolor','none')
xlim([-20,20])
ylim([-20,20])
view(2)
axis equal
axis tight
colorbar;
% colormap(brewermap(10,'Spectral'))
% caxis([min(nfoc_3(:)) max(nfoc_3(:))])
box on
set(gca,'linewidth',1,'fontsize',15)



subplot(2,3,5)
surf(q_x(2:end-1)+tempqx,q_y(2:end-1)+tempqx,u_theta_cart((2:end-1),(2:end-1)),'edgecolor','none')
xlim([-20,20])
ylim([-20,20])
view(2)
axis equal
axis tight
colorbar;
% colormap(brewermap(10,'Spectral'))
% caxis([min(nfoc_4(:)) max(nfoc_4(:))])
box on
set(gca,'linewidth',1,'fontsize',15)

subplot(2,3,6)
surf(q_x(2:end-1)+tempqx,q_y(2:end-1)+tempqx,u_z_cart((2:end-1),(2:end-1)),'edgecolor','none')
xlim([-20,20])
ylim([-20,20])
view(2)
axis equal
axis tight
colorbar;
% colormap(brewermap(10,'Spectral'))
% caxis([min(nfoc_5(:)) max(nfoc_5(:))])
box on
set(gca,'linewidth',1,'fontsize',15)





%% Begin binning in cylindrical coordinates

q_r = 1:25;
q_z = 0:50;
tempqr = 0.5*(q_r(2)-q_r(1));
tempqz = 0.5*(q_z(2)-q_z(1));


[~,bin_rad_edges ,idx_rad]  = histcounts(Rdata,q_r); 
[~,bin_z_edges ,idx_z]     = histcounts(Zdata,q_z);

u_theta = zeros(length(unique(idx_z)),length(unique(idx_rad)));
u_theta_std = u_theta;
u_theta_cnt = u_theta;

u_z = u_theta;
u_z_std = u_theta;
u_z_cnt = u_theta;

u_r = u_theta;
u_r_std = u_theta;
u_r_cnt = u_theta;

dum1 = unique(idx_z);
dum2 = unique(idx_rad);
for i=1:length(dum1)
    for j=1:length(dum2)
        
        u_theta(i,j) = mean(Uthdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        u_theta_std(i,j) = std(Uthdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        u_theta_cnt(i,j) = length(Uthdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        
        u_z(i,j) = mean(Uzdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        u_z_std(i,j) = std(Uzdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        u_z_cnt(i,j) = length(Uzdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        
        
        u_r(i,j) = mean(Urdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        u_r_std(i,j) = std(Urdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));
        u_r_cnt(i,j) = length(Urdata(ismember(idx_z,dum1(i)) & ismember(idx_rad,dum2(j))));

    end
end

figure(2)
subplot(1,3,1)
surf(q_z(2:end-1)+tempqz,q_r(2:end-1)+tempqr,u_r((2:end-1),(2:end-1))','edgecolor','none')
view(2)
axis equal
axis tight
colorbar

subplot(1,3,2)
surf(q_z(2:end-1)+tempqz,q_r(2:end-1)+tempqr,u_theta((2:end-1),(2:end-1))','edgecolor','none')
view(2)
axis equal
axis tight
colorbar

subplot(1,3,3)
surf(q_z(2:end-1)+tempqz,q_r(2:end-1)+tempqr,u_z((2:end-1),(2:end-1))','edgecolor','none')
view(2)
axis equal
axis tight
colorbar


%% Begin volumetric binning in time


[~,bin_rad_edges ,idx_rad]  = histcounts(Rdata,1:2:25); 
[~,bin_z_edges ,idx_z]     = histcounts(Zdata);

u_theta = zeros(length(unique(idx_z)),length(unique(idx_rad)));
u_theta_std = u_theta;
u_theta_cnt = u_theta;

% dum1 = unique(idx_z);
dum2 = unique(idx_rad);
% for i=1:length(dum1)
time_u_theta_vol = [];
time_u_theta_vol_std = [];
time_u_theta_vol_cnt = [];


time_u_z_vol = [];
time_u_z_vol_std = [];
time_u_z_vol_cnt = [];

time_u_theta_vol_quad1 = [];
time_u_theta_vol_std_quad1 = [];
time_u_theta_vol_cnt_quad1 = [];

time_u_theta_vol_quad2 = [];
time_u_theta_vol_std_quad2 = [];
time_u_theta_vol_cnt_quad2 = [];

time_u_theta_vol_quad3 = [];
time_u_theta_vol_std_quad3 = [];
time_u_theta_vol_cnt_quad3 = [];

time_u_theta_vol_quad4 = [];
time_u_theta_vol_std_quad4 = [];
time_u_theta_vol_cnt_quad4 = [];

twin=10;
tstart = min(unique(Timedata)): twin : max(unique(Timedata));
for i=1:length(tstart)
    i
    dumtime=ismember(Timedata,tstart(i):(tstart(i)+twin-1));
for j=1:length(dum2)
    
    
    % Type 1) Indiscriminate binning
    time_u_theta_vol(i,j) = mean(Uthdata(ismember(idx_rad,dum2(j)) & dumtime ));
    time_u_theta_vol_std(i,j) = std(Uthdata(ismember(idx_rad,dum2(j)) & dumtime ));
    time_u_theta_vol_cnt(i,j) = length(Uthdata(ismember(idx_rad,dum2(j)) & dumtime ));
    
    time_u_z_vol(i,j) = mean(Uzdata(ismember(idx_z,dum2(j)) & dumtime ));
    time_u_z_vol_std(i,j) = std(Uzdata(ismember(idx_z,dum2(j)) & dumtime ));
    time_u_z_vol_cnt(i,j) = length(Uzdata(ismember(idx_z,dum2(j)) & dumtime ));

    
    % Type 2) Binning utheta by Quadrants
    time_u_theta_vol_quad1(i,j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=-pi & Thetadata<(-pi + pi/2))   & dumtime ));
    time_u_theta_vol_std_quad1(i,j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=-pi & Thetadata<(-pi + pi/2))   & dumtime ));
    time_u_theta_vol_cnt_quad1(i,j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=-pi & Thetadata<(-pi + pi/2))    & dumtime ));
    
    time_u_theta_vol_quad2(i,j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=(-pi + pi/2) & Thetadata<(-pi + pi/2 + pi/2))   & dumtime ));
    time_u_theta_vol_std_quad2(i,j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + pi/2) & Thetadata<(-pi + pi/2 + pi/2))   & dumtime ));
    time_u_theta_vol_cnt_quad2(i,j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + pi/2) & Thetadata<(-pi + pi/2 + pi/2))    & dumtime ));
    
    time_u_theta_vol_quad3(i,j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=(-pi + 2*pi/2) & Thetadata<(-pi + 2*pi/2 + pi/2))   & dumtime ));
    time_u_theta_vol_std_quad3(i,j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 2*pi/2) & Thetadata<(-pi + 2*pi/2 + pi/2))   & dumtime ));
    time_u_theta_vol_cnt_quad3(i,j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 2*pi/2) & Thetadata<(-pi + 2*pi/2 + pi/2))    & dumtime ));
    
    time_u_theta_vol_quad4(i,j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=(-pi + 3*pi/2) & Thetadata<(-pi + 3*pi/2 + pi/2))   & dumtime ));
    time_u_theta_vol_std_quad4(i,j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 3*pi/2) & Thetadata<(-pi + 3*pi/2 + pi/2))   & dumtime ));
    time_u_theta_vol_cnt_quad4(i,j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 3*pi/2) & Thetadata<(-pi + 3*pi/2 + pi/2))    & dumtime ));
    
    
    
    
end

end
%%
v = VideoWriter('newfile.avi','Uncompressed AVI');
v.FrameRate = 10;
open(v)
for i=1:size(time_u_theta_vol,1)
    
    figure(77)
    subplot(1,2,1)
    title({['Time: ' num2str(tstart(i)./5000) ' sec,'],['$\Delta t=$ ' num2str(twin./5000) ' seconds' ]},'interpreter','latex','fontsize',12)
    hold on
%     plot(bin_rad_edges(2:end-1)./1000,time_u_theta_vol(i,2:end-1)./1000,'m-.','linewidth',1.5)
    plot(bin_rad_edges(2:end-1)./1000,time_u_theta_vol_quad1(i,2:end-1)./1000,'r-.','linewidth',1.5)
    plot(bin_rad_edges(2:end-1)./1000,time_u_theta_vol_quad2(i,2:end-1)./1000,'b-.','linewidth',1.5)
    plot(bin_rad_edges(2:end-1)./1000,time_u_theta_vol_quad3(i,2:end-1)./1000,'g-.','linewidth',1.5)
    plot(bin_rad_edges(2:end-1)./1000,time_u_theta_vol_quad4(i,2:end-1)./1000,'k-.','linewidth',1.5)
    hold off
    box on
    set(gca,'linewidth',1,'fontsize',15)
    % legend('$N_{foc}=0$','$N_{foc}=3$','$N_{foc}=4$','$N_{foc}=5$','interpreter','latex','fontsize',20)
    xlabel('$r$ (m)','interpreter','latex','fontsize',16)
    ylabel('$u_\theta$ (m/s)','interpreter','latex','fontsize',16)
%     legend('$-\pi \leq \theta < -\pi/2 $','$-\pi/2 \leq \theta < 0$','$0 \leq \theta < \pi/2 $','$\pi/2 \leq \theta < \pi $','interpreter','latex','fontsize',12)
    axis square
    ylim([0 1.5])

    % Plot particle count

    subplot(1,2,2)
    hold on
    plot(bin_rad_edges(2:end-1)./1000,(time_u_theta_vol_cnt_quad1(i,2:end-1)),'r-.','linewidth',1.5)
    plot(bin_rad_edges(2:end-1)./1000,(time_u_theta_vol_cnt_quad2(i,2:end-1)),'b-.','linewidth',1.5)
    plot(bin_rad_edges(2:end-1)./1000,(time_u_theta_vol_cnt_quad3(i,2:end-1)),'g-.','linewidth',1.5)
    plot(bin_rad_edges(2:end-1)./1000,(time_u_theta_vol_cnt_quad4(i,2:end-1)),'k-.','linewidth',1.5)
    hold off
    box on
    set(gca,'linewidth',1,'fontsize',15)
    % legend('$N_{foc}=0$','$N_{foc}=3$','$N_{foc}=4$','$N_{foc}=5$','interpreter','latex','fontsize',20)
    xlabel('$r$ (m)','interpreter','latex','fontsize',12)
    ylabel('$[\#]$ Particles','interpreter','latex','fontsize',12)
    legend('$-\pi \leq \theta < -\pi/2 $','$-\pi/2 \leq \theta < 0$','$0 \leq \theta < \pi/2 $','$\pi/2 \leq \theta < \pi $','interpreter','latex','fontsize',12,'location','northoutside')
    axis square
    ylim([0 2000])
    drawnow
%     pause

    A = getframe(gcf);
    writeVideo(v,A)
    clf
end
close(v)




%% Begin volumetric binning


[~,bin_rad_edges ,idx_rad]  = histcounts(Rdata,1:25); 
[~,bin_z_edges ,idx_z]     = histcounts(Zdata, 0:50);

u_theta = zeros(length(unique(idx_z)),length(unique(idx_rad)));
u_theta_std = u_theta;
u_theta_cnt = u_theta;

% dum1 = unique(idx_z);
dum2 = unique(idx_rad);
% for i=1:length(dum1)
u_theta_vol = [];
u_theta_vol_std = [];
u_theta_vol_cnt = [];
u_z_vol = [];
u_z_vol_std = [];
u_z_vol_cnt = [];

for j=1:length(dum2)
    
    
    % Type 1) Indiscriminate binning
    u_theta_vol(j) = mean(Uthdata(ismember(idx_rad,dum2(j))));
    u_theta_vol_std(j) = std(Uthdata(ismember(idx_rad,dum2(j))));
    u_theta_vol_cnt(j) = length(Uthdata(ismember(idx_rad,dum2(j))));
    
    u_z_vol(j) = mean(Uzdata(ismember(idx_z,dum2(j))));
    u_z_vol_std(j) = std(Uzdata(ismember(idx_z,dum2(j))));
    u_z_vol_cnt(j) = length(Uzdata(ismember(idx_z,dum2(j))));

    
    % Type 2) Binning utheta by Quadrants
    u_theta_vol_quad1(j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=-pi & Thetadata<(-pi + pi/2))   ));
    u_theta_vol_std_quad1(j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=-pi & Thetadata<(-pi + pi/2))   ));
    u_theta_vol_cnt_quad1(j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=-pi & Thetadata<(-pi + pi/2))    ));
    
    u_theta_vol_quad2(j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=(-pi + pi/2) & Thetadata<(-pi + pi/2 + pi/2))   ));
    u_theta_vol_std_quad2(j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + pi/2) & Thetadata<(-pi + pi/2 + pi/2))   ));
    u_theta_vol_cnt_quad2(j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + pi/2) & Thetadata<(-pi + pi/2 + pi/2))    ));
    
    u_theta_vol_quad3(j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=(-pi + 2*pi/2) & Thetadata<(-pi + 2*pi/2 + pi/2))   ));
    u_theta_vol_std_quad3(j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 2*pi/2) & Thetadata<(-pi + 2*pi/2 + pi/2))   ));
    u_theta_vol_cnt_quad3(j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 2*pi/2) & Thetadata<(-pi + 2*pi/2 + pi/2))    ));
    
    u_theta_vol_quad4(j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & (Thetadata>=(-pi + 3*pi/2) & Thetadata<(-pi + 3*pi/2 + pi/2))   ));
    u_theta_vol_std_quad4(j) = std(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 3*pi/2) & Thetadata<(-pi + 3*pi/2 + pi/2))   ));
    u_theta_vol_cnt_quad4(j) = length(Uthdata(ismember(idx_rad,dum2(j)) & (Thetadata>=(-pi + 3*pi/2) & Thetadata<(-pi + 3*pi/2 + pi/2))    ));
    
    
    
    % Type 3) Binning utheta by Focality
    u_theta_vol_foc0(j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & ismember(Focalitydata,0)   ));
    u_theta_vol_std_foc0(j) = std(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,0)    ));
    u_theta_vol_cnt_foc0(j) = length(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,0)    ));
    
    
    u_theta_vol_foc3(j) = mean(Uthdata(ismember(idx_rad,dum2(j))   & ismember(Focalitydata,3)   ));
    u_theta_vol_std_foc3(j) = std(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,3)    ));
    u_theta_vol_cnt_foc3(j) = length(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,3)    ));
    
    u_theta_vol_foc4(j) = mean(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,4)   ));
    u_theta_vol_std_foc4(j) = std(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,4)  ));
    u_theta_vol_cnt_foc4(j) = length(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,4) ));
    
    u_theta_vol_foc5(j) = mean(Uthdata(ismember(idx_rad,dum2(j)) & ismember(Focalitydata,5) ));
    u_theta_vol_std_foc5(j) = std(Uthdata(ismember(idx_rad,dum2(j))& ismember(Focalitydata,5) ));
    u_theta_vol_cnt_foc5(j) = length(Uthdata(ismember(idx_rad,dum2(j))& ismember(Focalitydata,5) ));
   
    
end



%% Further analyses of velocity profiles binned indiscriminately


% Fit Lamb-Oseen model
alpha = 1.2526; % a constant for centering the peak of a vortex profile at r_v
uth_lo = fittype(@(lambda_inf,r_v,x)(lambda_inf./(2.*pi.*x).*(1 - exp(-alpha.*x.^2./r_v.^2)) ));% - (u_theta_vol(2:end)./1000));
init_lo = [0.05 7e-3];
[lo_fitted,gof] = fit(bin_rad_edges(2:end)'./1000,(u_theta_vol(2:end)./1000)',uth_lo,'StartPoint',init_lo)

% Fit Bosschers Analytical Model
uth_bo = fittype(@(lambda_inf,r_v,x)((lambda_inf./(2.*pi.*x)).*(1  -  ((r_v.^2./(r_v.^2 + alpha.*(r_c_m./1000).^2)).*exp(-alpha.*(x.^2 - (r_c_m./1000).^2)./r_v.^2)) )));
[bo_fitted,gof] = fit(bin_rad_edges(2:end)'./1000,(u_theta_vol(2:end)./1000)',uth_bo,'StartPoint',init_lo)

% Fit Proctor's model
% Part 1) The irrotational region
uth_pr1 = fittype(@(lambda_inf,BB,beta,x)(lambda_inf./(2.*pi.*x).*(1 - exp(-beta.* (x./BB).^0.75   )) ));% - (u_theta_vol(2:end)./1000));
init_pr1 = [(0.5*0.1256*0.65*u_z_vol(end)./1000) 0.3 14];
[pr1_fitted,gof] = fit(bin_rad_edges(end-15:end)'./1000,(u_theta_vol(end-15:end)./1000)',uth_pr1,'StartPoint',init_pr1)
% Part 2) The viscous core region
uth_pr2 = fittype(@(r_v,x)(  (1.0939.*pr1_fitted.lambda_inf./ (2.*pi.*x) )   .*(1 - exp(-pr1_fitted.beta.* (1.4.*r_v./pr1_fitted.BB).^0.75   ))  .*(1 - exp(-alpha.* (x./r_v).^2   ))   ));
init_pr2 = 1e-2;
[pr2_fitted,gof] = fit(bin_rad_edges(2:7)'./1000,(u_theta_vol(2:7)./1000)',uth_pr2,'StartPoint',init_pr2)

% Fit Bosschers Semi-empirical model (9-parameter, holy shit)
uth_bo2 = fittype(@(lambda_inf,alphacap,beta,BB,zeta_1,zeta_2,p,q,r_v , x) ...
                   ( (lambda_inf./(2.*pi.*x))  .*  (1 - alphacap.*exp(-zeta_1.* (x./BB).^p   ))  ...
                      .*  (1 - beta.*exp(-zeta_2.* ((x.^q - (r_c_m./1000).^q)./r_v.^q)   ))    ));
                  
init_bo2 = [(0.5*0.1256*0.65*u_z_vol(end)./1000) 0.5 1 0.2 6 1.2526 0.75 2 0.01];
[bo2_fitted,gof] = fit(bin_rad_edges(2:end)'./1000,(u_theta_vol(2:end)./1000)',uth_bo2,'StartPoint',init_bo2)

% Compare fitted beta with intended value for zero shear stress (eqns 3.33 and 3.34 of his PhD thesis)
% alphabar = bo2_fitted.alphacap*exp(-bo2_fitted.zeta_1*((r_c_m/1000)/bo2_fitted.BB)^bo2_fitted.p); 
% beta_zeroshear = ;
                  


% Define Pepijn's data (Winkelman Model) coefficients
% lambda_pp = 0.5*0.1256*0.65*u_z_vol(end)./1000;
% rad_space = linspace(r_c_m./1000,25e-3,100);
% beta_i = 1.3e4;
% beta_o = 13.5;
% B_pp   = 0.3;
% p_pp   = 4;
% uth_pp = lambda_pp./(2.*pi.*rad_space).*(1 - exp( (-beta_i.* (rad_space./B_pp).^2) ./ ( 1 + ((beta_i./beta_o) .* (rad_space./B_pp).^(5/4)  ).^p_pp     ).^(1./p_pp)        )) ;


% % Plot results
% scatter(bin_rad_edges(2:end)'./1000,(u_theta_vol(2:end)./1000)', 'r+')
% hold on
% plot(linspace(1e-4,25e-3,100)',fitted_curve(linspace(1e-4,25e-3,100))')
% hold off
        


figure(3)

% subplot(1,3,1)

hold all
yyaxis left
% Plot Data
% plot(bin_rad_edges./1000,(u_theta_vol./1000))
errorbar(bin_rad_edges./1000,u_theta_vol./1000,2*u_theta_vol_std./sqrt(u_theta_vol_cnt)./1000,'kd','linewidth',1.5)
% Plot Models
rad_space = linspace(r_c_m./1000,25e-3,100); % define radial space
rad_space_pr1 = linspace(1.4*pr2_fitted.r_v,25e-3,100); % define radial space
rad_space_pr2 = linspace(r_c_m./1000,1.4*pr2_fitted.r_v,100); % define radial space

plot(rad_space',lo_fitted(rad_space)','b-','linewidth',1.5) % lamb-oseen
plot(rad_space',bo_fitted(rad_space)','r-','linewidth',1.5) % bosscher
plot(rad_space',bo2_fitted(rad_space)','r-.','linewidth',1.5) % bosscher
plot(rad_space_pr1',pr1_fitted(rad_space_pr1)','m--','linewidth',1.5) % proctor
plot(rad_space_pr2',pr2_fitted(rad_space_pr2)','m-.','linewidth',1.5) % proctor

ylim([0 1.200])
vline(r_c_m./1000,'k-')
vline(r_c_m./1000 + 2*r_c_std./1000,'k--')
vline(r_c_m./1000 - 2*r_c_std./1000,'k--')
xlabel('$r$ (m)','interpreter','latex','fontsize',20)
ylabel('$u_\theta$ (m/s)','interpreter','latex','fontsize',20)

yyaxis right
bar(bin_rad_edges(2:end)./1000,u_theta_vol_cnt(2:end),'r','edgecolor','none','facealpha',0.3)
ylabel('$[\#]$ particles binned','interpreter','latex','fontsize',20)



legend('Data','Lamb-Oseen','Bosschers-Analytical','Bosschers-Semi-Empirical','Proctors-1','Proctors-2','interpreter','latex','fontsize',15)
box on
set(gca,'linewidth',1,'fontsize',15)




figure(4)

% subplot(1,3,2)

hold all
yyaxis left
plot(bin_rad_edges./1000,u_z_vol./1000,'linewidth',1.5)
errorbar(bin_rad_edges./1000,u_z_vol./1000,2*u_z_vol_std./sqrt(u_z_vol_cnt)./1000,'linewidth',1.5)
ylim([4.0 6.50])
vline(r_c_m./1000,'k-')
vline(r_c_m./1000 + 2*r_c_std./1000,'k--')
vline(r_c_m./1000 - 2*r_c_std./1000,'k--')
xlabel('$r$ (m)','interpreter','latex','fontsize',20)
ylabel('$u_z$ (m/s)','interpreter','latex','fontsize',20)



yyaxis right
bar(bin_rad_edges(2:end)./1000,u_theta_vol_cnt(2:end),'r','edgecolor','none','facealpha',0.3)
ylabel('$[\#]$ particles binned','interpreter','latex','fontsize',20)


box on
set(gca,'linewidth',1,'fontsize',15)





figure(5)
% subplot(1,3,3)
hold all
yyaxis left
plot(bin_rad_edges./1000,u_z_vol_std./1000,'linewidth',1.5)
% ylim([4.0 6.50])
vline(r_c_m./1000,'k-')
vline(r_c_m./1000 + 2*r_c_std./1000,'k--')
vline(r_c_m./1000 - 2*r_c_std./1000,'k--')
xlabel('$r$ (m)','interpreter','latex','fontsize',20)
ylabel('$u_z\prime$ (m/s)','interpreter','latex','fontsize',20)



yyaxis right
bar(bin_rad_edges(2:end)./1000,u_theta_vol_cnt(2:end),'r','edgecolor','none','facealpha',0.3)
ylabel('$[\#]$ particles binned','interpreter','latex','fontsize',20)


box on
set(gca,'linewidth',1,'fontsize',15)





%% Examine velocity profile by focality and region


figure(6)

% Plot tangential profile by focality

subplot(2,2,1)
hold all
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_foc0(2:end-1)./1000,2*u_theta_vol_std_foc0(2:end-1)./sqrt(u_theta_vol_cnt_foc0(2:end-1))./1000,'r-','linewidth',1.5)
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_foc3(2:end-1)./1000,2*u_theta_vol_std_foc3(2:end-1)./sqrt(u_theta_vol_cnt_foc3(2:end-1))./1000,'b-','linewidth',1.5)
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_foc4(2:end-1)./1000,2*u_theta_vol_std_foc4(2:end-1)./sqrt(u_theta_vol_cnt_foc4(2:end-1))./1000,'g-','linewidth',1.5)
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_foc5(2:end-1)./1000,2*u_theta_vol_std_foc5(2:end-1)./sqrt(u_theta_vol_cnt_foc5(2:end-1))./1000,'k-','linewidth',1.5)
box on
set(gca,'linewidth',1,'fontsize',15)
legend('$N_{foc}=0$','$N_{foc}=3$','$N_{foc}=4$','$N_{foc}=5$','interpreter','latex','fontsize',20)
xlabel('$r$ (m)','interpreter','latex','fontsize',20)
ylabel('$u_\theta$ (m/s)','interpreter','latex','fontsize',20)

% Plot particle count

subplot(2,2,2)
hold all
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_foc0(2:end-1)),'r-.','linewidth',1.5)
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_foc3(2:end-1)),'b-.','linewidth',1.5)
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_foc4(2:end-1)),'g-.','linewidth',1.5)
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_foc5(2:end-1)),'k-.','linewidth',1.5)
box on
set(gca,'linewidth',1,'fontsize',15)
% legend('$N_{foc}=0$','$N_{foc}=3$','$N_{foc}=4$','$N_{foc}=5$','interpreter','latex','fontsize',20)
xlabel('$r$ (m)','interpreter','latex','fontsize',20)
ylabel('$[\#]$ Particles','interpreter','latex','fontsize',20)
% legend('$-\pi \leq \theta < -\pi/2 $','$-\pi/2 \leq \theta < 0$','$0 \leq \theta < \pi/2 $','$\pi/2 \leq \theta < \pi $','interpreter','latex','fontsize',20)

% Plot tangential profile by quadrant

subplot(2,2,3)
hold all
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_quad1(2:end-1)./1000,2*u_theta_vol_std_quad1(2:end-1)./sqrt(u_theta_vol_cnt_quad1(2:end-1))./1000,'r-.','linewidth',1.5)
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_quad2(2:end-1)./1000,2*u_theta_vol_std_quad2(2:end-1)./sqrt(u_theta_vol_cnt_quad2(2:end-1))./1000,'b-.','linewidth',1.5)
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_quad3(2:end-1)./1000,2*u_theta_vol_std_quad3(2:end-1)./sqrt(u_theta_vol_cnt_quad3(2:end-1))./1000,'g-.','linewidth',1.5)
errorbar(bin_rad_edges(2:end-1)./1000,u_theta_vol_quad4(2:end-1)./1000,2*u_theta_vol_std_quad4(2:end-1)./sqrt(u_theta_vol_cnt_quad4(2:end-1))./1000,'k-.','linewidth',1.5)
box on
set(gca,'linewidth',1,'fontsize',15)
% legend('$N_{foc}=0$','$N_{foc}=3$','$N_{foc}=4$','$N_{foc}=5$','interpreter','latex','fontsize',20)
xlabel('$r$ (m)','interpreter','latex','fontsize',20)
ylabel('$u_\theta$ (m/s)','interpreter','latex','fontsize',20)
legend('$-\pi \leq \theta < -\pi/2 $','$-\pi/2 \leq \theta < 0$','$0 \leq \theta < \pi/2 $','$\pi/2 \leq \theta < \pi $','interpreter','latex','fontsize',20)

% Plot particle count

subplot(2,2,4)
hold all
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_quad1(2:end-1)),'r-.','linewidth',1.5)
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_quad2(2:end-1)),'b-.','linewidth',1.5)
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_quad3(2:end-1)),'g-.','linewidth',1.5)
plot(bin_rad_edges(2:end-1)./1000,(u_theta_vol_cnt_quad4(2:end-1)),'k-.','linewidth',1.5)
box on
set(gca,'linewidth',1,'fontsize',15)
% legend('$N_{foc}=0$','$N_{foc}=3$','$N_{foc}=4$','$N_{foc}=5$','interpreter','latex','fontsize',20)
xlabel('$r$ (m)','interpreter','latex','fontsize',20)
ylabel('$[\#]$ Particles','interpreter','latex','fontsize',20)
% legend('$-\pi \leq \theta < -\pi/2 $','$-\pi/2 \leq \theta < 0$','$0 \leq \theta < \pi/2 $','$\pi/2 \leq \theta < \pi $','interpreter','latex','fontsize',20)





end

