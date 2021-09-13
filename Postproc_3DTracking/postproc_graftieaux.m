function postproc_graftieaux


%% get globals
global folder date rec prop mset post%plotting

%% Load data
% load grid
vecgrid=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'vecgrid.mat']);

% load general data
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Time.dat']);
GridPosition=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Position.dat']);
load([folder date rec vsl 'Postproc3DTracking' vsl 'TrackQualityMetric.mat'],'t_spurious');

% Load data
BinnedVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedVelocity.dat']);



% Load visual hull data
VHData = load([folder date rec vsl 'Postproc3DReco' vsl 'VHReduced.mat']);

% Find mean cavity radius
r_c_m = nanmean(0.5*(4*VHData.Visual_hull_reduced(7,:,:)./pi).^0.5,[2 3]);
r_c_std = nanstd(0.5*(4*VHData.Visual_hull_reduced(7,:,:)./pi).^0.5,0,[2 3]);

%% create path for processing function storage
if 0
    if exist([folder date rec vsl 'Vortex3DTracking' vsl 'VortexReference'],'dir')~=7
        mkdir([folder date rec vsl 'Vortex3DTracking' vsl 'VortexReference'])
    end

    if exist([folder date rec vsl 'Vortex3DTracking' vsl 'CavityReference'],'dir')~=7
        mkdir([folder date rec vsl 'Vortex3DTracking' vsl 'CavityReference'])
    end

end


% loop frames
for n=post.tproc(1):post.tproc(2)%mset.tspan(1):mset.tspan(2)
    disp(['Starting tstep ' num2str(n)])
    % Time set
    N=Index(2,:)==n ;
    dumInd = ~ismember(Index(1,:),t_spurious);
    % Grid indexing
    I=ismember(vecgrid.IND,Index(1,N));
    
    
    velgrid = nan(size(vecgrid.X));
    velgrid(:,I) =  BinnedVelocity(1:3,N);
    
    %%% Following lines commented cuz I found buggy
%     % Grid node
%     gridnode=vecgrid.X(:,I);
%     % get grid positions
%     pos=gridnode;
%     % select volume
%     seg=true(1,size(gridnode,2));
%     
%     % get velocity
%     vel = BinnedVelocity(1:3,N);
    
    gridnode = vecgrid.X;
    pos = gridnode;
    vel = velgrid;
    
    
    
    %% Estimate graftieaux function

    
    % get unique z-positions
    unq_zpos = unique(pos(3,:));
    % get total number of z-planes
    n_zslices = length(unique(pos(3,:)));

    % create a new, larger set of grid coordinate (for the vortex data only)
    % step 1: make meshgrid
    gres = unq_zpos(2)-unq_zpos(1); % calculate grid-resolution from data (grid must be isometric)
    [posx_new,posy_new] = meshgrid(-40:gres:40,-40:gres:40); % 40 chosen for vortex test-case only
    % step 2: unwrap meshgrid to column-vector format (use 3rd
    % dimension to bin data over time and slice
    pos_new(1,:) = reshape(posx_new,[1 (size(posx_new,1)*size(posx_new,2))]);
    pos_new(2,:) = reshape(posy_new,[1 (size(posx_new,1)*size(posx_new,2))]);
    
    if n==post.tproc(1)
        Ux_bin = zeros(n_zslices,size(pos_new,2),length(post.tproc(1):post.tproc(2)));
        Uy_bin = Ux_bin;
        Uz_bin = Ux_bin;
        num_bin = Ux_bin;
    end
    
    
    % make accumulation matrix
    vec_accum_mat = nan(size(pos_new,1),size(pos_new,2),n_zslices);
    pos_accum_mat = nan(size(pos_new,1),size(pos_new,2),n_zslices);

    % estimate the graftieaux function for each slice and
    % reposition
    for iii=1:n_zslices
        
        
        %% Part 1: Graftieaux function
%         disp([ num2str(100*(iii-1)/n_zslices) '% slice-data done'])
        % step 1: find logical index of all vectors that belong to
        % this z-slice
        slc_idx = pos(3,:)== unq_zpos(iii);
        if 0
            % step 2: calculate the direction-vector matrix
            % (0-diagonal)
            dirx_mat = -minus(pos(1,slc_idx),pos(1,slc_idx)');
            diry_mat = -minus(pos(2,slc_idx),pos(2,slc_idx)');

            % step 3: reshape the velocity vector data
            velx_slc = vel(1,slc_idx);
            vely_slc = vel(2,slc_idx);


            velx_mat = repmat(velx_slc,length(velx_slc),1);
            vely_mat = repmat(vely_slc,length(vely_slc),1);

            vel_iii(:,:,1) = velx_mat;
            vel_iii(:,:,2) = vely_mat;
            vel_iii(:,:,3) = zeros(size(vely_mat));

            vel_iii = vel_iii./vecnorm(vel_iii,2,3); % normalize

            dir_iii(:,:,1) = dirx_mat;
            dir_iii(:,:,2) = diry_mat;
            dir_iii(:,:,3) = zeros(size(dirx_mat));

            dir_iii = dir_iii./vecnorm(dir_iii,2,3); % normalize

            % step 5.a.: perform cross_-product
            cross_prod = cross(dir_iii,vel_iii,3);

            % step 5.b.: calculate vector norm
            sine_theta = vecnorm(cross_prod,2,3);

            % step 6: sum over all each row entries
            fcn_graftieaux = nanmean(sine_theta,2);

            [~,indMax(iii)] = max(fcn_graftieaux); % this finds the index in only the slice data

            % reposition slice with vortex as centre
            slc_pos(1,:) = pos(1,slc_idx);
            slc_pos(2,:) = pos(2,slc_idx);
            vel_pos(1,:) = vel(1,slc_idx);
            vel_pos(2,:) = vel(2,slc_idx);
            vel_pos(3,:) = vel(3,slc_idx);

            slc_pos(1,:) = slc_pos(1,:)-slc_pos(1,indMax(iii));
            slc_pos(2,:) = slc_pos(2,:)-slc_pos(2,indMax(iii));

            if 0
                figure(1);hold all;%scatter(slc_pos(1,:),slc_pos(2,:))
                quiver(slc_pos(1,:),slc_pos(2,:),vel(1,slc_idx),vel(2,slc_idx),'b');
                drawnow
            end
    %         pause
    %         
            % intersect binning positions with shifted vortex positions
            [~,ia,ib] = intersect(slc_pos',pos_new','rows'); % slc_pos(ia) = pos_new(ib)

            if 0
                    figure(1)
                    scatter(slc_pos(1,:),slc_pos(2,:),'r+')
                    hold on
                    scatter(pos_new(1,ib),pos_new(2,ib),'bo')
                    drawnow
                    pause

                    figure(1);hold all;plot(ia),plot(ib);drawnow;pause
                    figure(1);hold all;scatter(slc_pos(1,:),slc_pos(2,:),'r+')
                    scatter(pos_new(1,:),pos_new(2,:),'bo')
                    drawnow
                    pause
                    clf
                    size(slc_pos)
                    size(pos_new)
                    pause
                    ia_correct = ia<size(vel_pos,2) ; % remove some spurious stuff (check later?)
                    ind_rem = find(ia<size(vel_pos,2));
                    size(ia_correct)
                    size(ib)
                    size(ia)
                    pause

                    plot(ia>size(vel_pos,2))
                    pause
                    plot(ia_correct>size(vel_pos,2));drawnow;pause

            end



            % Accumulate velocities     
            vec_accum_mat(1,ib,iii) = vel_pos(1,ia);
            vec_accum_mat(2,ib,iii) = vel_pos(2,ia);
            vec_accum_mat(3,ib,iii) = vel_pos(3,ia);

            % Make grid of positions
            pos_accum_mat(1,ib,iii) = pos_new(1,ia);
            pos_accum_mat(2,ib,iii) = pos_new(2,ia);
            pos_accum_mat(3,ib,iii) = ones(size(pos_new(1,ia))).*unq_zpos(iii);
        end
        %% Part 2: Cavity centered coordinate system 
        if 1
            %% For visual hull processing
            % Load visual hull and estimate axis
            load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(n) '.mat'],'Centroid_cav')
            % determine radial position
%             centre_x = interp1(Centroid_cav(3,:),Centroid_cav(1,:),unq_zpos(iii));
%             centre_y = interp1(Centroid_cav(3,:),Centroid_cav(2,:),unq_zpos(iii));
            
            xcenfit = fit(fillmissing(Centroid_cav(3,:)','linear'),fillmissing(Centroid_cav(1,:)','linear'),'poly1');
            ycenfit = fit(fillmissing(Centroid_cav(3,:)','linear'),fillmissing(Centroid_cav(2,:)','linear'),'poly1');
            
            centre_x = xcenfit(unq_zpos(iii));
            centre_y = ycenfit(unq_zpos(iii));
            
            slc_pos(1,:) = pos(1,slc_idx);
            slc_pos(2,:) = pos(2,slc_idx);
            vel_pos(1,:) = vel(1,slc_idx);
            vel_pos(2,:) = vel(2,slc_idx);
            vel_pos(3,:) = vel(3,slc_idx);

            slc_pos(1,:) = slc_pos(1,:)-centre_y;
            slc_pos(2,:) = slc_pos(2,:)-centre_x;
%             
            
            
            vec_accum_mat_cc(1,:,iii) = vel_pos(1,:);
            vec_accum_mat_cc(2,:,iii) = vel_pos(2,:);
            vec_accum_mat_cc(3,:,iii) = vel_pos(3,:);

            pos_accum_mat_cc(1,:,iii) = slc_pos(1,:);
            pos_accum_mat_cc(2,:,iii) = slc_pos(2,:);
            pos_accum_mat_cc(3,:,iii) = ones(size(slc_pos(1,:))).*unq_zpos(iii);

            if 0
                figure(1);hold all;scatter(pos_new(1,:),pos_new(2,:),'b+')
                quiver(slc_pos(1,:),slc_pos(2,:),vel(1,slc_idx),vel(2,slc_idx),'r');
                drawnow
                axis equal
%                 pause
            end
            %% Binning the centered velocity field in 3D to one grid
            % find index of slc_pos that are close to pos_new
            [Idx,d] = knnsearch(pos_new',slc_pos');
            
            
            Ux_bin(iii,Idx,n) = vel_pos(1,:);
            Uy_bin(iii,Idx,n) = vel_pos(2,:);
            Uz_bin(iii,Idx,n) = vel_pos(3,:);
            
            
        end
        %%
%         pause
        % clear vars for next slice
        vel_iii=[];
        dir_iii=[];
        slc_pos=[];
        vel_pos=[];
        ia=[];
        ib=[];
        
        
    end % n_zslices
%     pause
%     quiver3(pos_accum_mat_cc(1,:,:),pos_accum_mat_cc(2,:,:),pos_accum_mat_cc(3,:,:),vec_accum_mat_cc(1,:,:),vec_accum_mat_cc(2,:,:),vec_accum_mat_cc(3,:,:),100)

    
    %% save data at end of each timestep
%     pause
    if 0
        disp(['Storing tstep ' num2str(n) ' completed!'])
        if 0
            save([folder date rec vsl 'Vortex3DTracking' vsl 'VortexReference' vsl 'Tstep_' num2str(n) '.mat'],'pos_accum_mat','vec_accum_mat')
        end
        if 1
            save([folder date rec vsl 'Vortex3DTracking' vsl 'CavityReference' vsl 'Tstep_' num2str(n) '.mat'],'pos_accum_mat_cc','vec_accum_mat_cc')
        end
    end
    vec_accum_mat = [];
    vec_accum_mat_cc = [];
    pos_accum_mat = [];
    pos_accum_mat_cc = [];
    
    
    
end




%% Reshape mean velocity

% Calculate mean
mean_Ux = nanmean(Ux_bin,3);
mean_Uy = nanmean(Uy_bin,3);
mean_Uz = nanmean(Uz_bin,3);
% Calculate number of samples
dumx = sum(double(~isnan(Ux_bin)),3);
dumy = sum(double(~isnan(Uy_bin)),3);
dumz = sum(double(~isnan(Uz_bin)),3);

Nsamples_x = reshape(dumx,[size(mean_Ux,1) size(posx_new)]);
Nsamples_y = reshape(dumy,[size(mean_Uy,1) size(posx_new)]);
Nsamples_z = reshape(dumz,[size(mean_Uz,1) size(posx_new)]);


Ux_m = reshape(mean_Ux,[size(mean_Ux,1) size(posx_new)]).*0.001; % reshape and convert to mm
Uy_m = reshape(mean_Uy,[size(mean_Uy,1) size(posx_new)]).*0.001; % reshape and convert to mm
Uz_m = reshape(mean_Uz,[size(mean_Uz,1) size(posx_new)]).*0.001; % reshape and convert to mm
 

%% Calculate divergence of mean velocity

[dudx,dudy,dudz] = gradient(permute(Ux_m,[2 3 1]),squeeze(posx_new(1,:))./1000,squeeze(posy_new(:,1))./1000,squeeze(unq_zpos)./1000);

[dvdx,dvdy,dvdz] = gradient(permute(Uy_m,[2 3 1]),squeeze(posx_new(1,:))./1000,squeeze(posy_new(:,1))./1000,squeeze(unq_zpos)./1000);

[dwdx,dwdy,dwdz] = gradient(permute(Uz_m,[2 3 1]),squeeze(posx_new(1,:))./1000,squeeze(posy_new(:,1))./1000,squeeze(unq_zpos)./1000);

% find union of all nan positions
common_nans = double(~isnan(dudx)) .* double(~isnan(dvdy)) .* double(~isnan(dwdz));
common_zeros = double(~ismember(dudx,0)) .* double(~ismember(dvdy,0)) .* double(~ismember(dwdz,0));


div_xaxis =  dudx(common_nans==1 & common_zeros==1);
div_yaxis = -(dvdy(common_nans==1 & common_zeros==1)+dwdz(common_nans==1 & common_zeros==1));


figure(19)
hold all
histogram2(div_xaxis,div_yaxis,linspace(-100,100,50),linspace(-100,100,50),'facecolor','flat','displaystyle','tile');view(2)
plot(linspace(-100,100,50),linspace(-100,100,50),'k-','linewidth',2)
plot(linspace(-100,100,50),linspace(-100,100,50)+50,'k--','linewidth',2)
plot(linspace(-100,100,50),linspace(-100,100,50)-50,'k--','linewidth',2)
xlabel('$\partial \overline{u}/\partial x$','interpreter','latex','fontsize',14)
ylabel('$-(\partial \overline{v}/\partial y + \partial \overline{w}/\partial z)$','interpreter','latex','fontsize',14)
ylim([-100 100])
xlim([-100 100])
box on
cb69 = colorbar;
cb69.Location = 'northoutside';
colormap(brewermap([],'Blues'))
set(gca,'linewidth',2)


% Calculate mean pressure gradient assuming steady and inviscid case
rho_l=1000;
dpdx = -(permute(Ux_m,[2 3 1]).*dudx + permute(Uy_m,[2 3 1]).*dudy + permute(Uz_m,[2 3 1]).*dudz);
dpdy = -(permute(Ux_m,[2 3 1]).*dvdx + permute(Uy_m,[2 3 1]).*dvdy + permute(Uz_m,[2 3 1]).*dvdz);
% quiver(posy_new,posx_new,imgaussfilt(nanmean(dpdx,3),1),imgaussfilt(nanmean(dpdy,3),1),5,'k')





% Calculate mean vorticity field
omega_x = dwdy - dvdz;
omega_y = dudz - dwdx;
omega_z = dvdx - dudy;


% Integrate pressure using the spectral decomposition algorithm by Cheng et al.
[pres,~]=SD_FPI(rho_l.*nanmean(dpdx,3),rho_l.*nanmean(dpdy,3),0.001,0.001);

% Add constant to the perimeter
pres = pres - (pres(1,1)-24500);

figure(69)
contourf(squeeze(posx_new(1,:))./1000,squeeze(posy_new(:,1))./1000,pres,10)
xlabel('$x_{cav}$ (m)','interpreter','latex','fontsize',14)
ylabel('$y_{cav}$ (m)','interpreter','latex','fontsize',14)



%% Radial binning
 
% Convert to polar
r = (posx_new.^2+posy_new.^2).^0.5;
% theta = atan2(posy_new,posx_new);
theta = cart2pol(posx_new,posy_new);
r = permute(repmat(r,1,1,size(Ux_m,1)),[3 1 2]);
theta = permute(repmat(theta,1,1,size(Ux_m,1)),[3 1 2]);
 
ur      = Uy_m.*sin(theta) + Ux_m.*cos(theta);
utheta  = -Ux_m.*sin(theta) + Uy_m.*cos(theta);
 
% 
% THETA = cart2pol(X,Y);
% Ur = U .* cos(THETA) + V.* sin(THETA);
% Uth = V .* cos(THETA) - U.* sin(THETA);
 
% bin radius 
bin_rad = 0.5:1:25;
[~,~,Idx] = histcounts(reshape(r,[size(r,1) size(r,2)*size(r,3)]),bin_rad);
u_theta_bin_mean = [];
u_theta_bin_std = [];
u_theta_bin_cnt = [];
 
uth_temp = reshape(utheta,[size(utheta,1) size(utheta,2)*size(utheta,3)]);
ur_temp = reshape(ur,[size(ur,1) size(ur,2)*size(ur,3)]);
 
% convert theta to degrees for.... convenience 
theta_2pi = wrapTo2Pi(theta);
 
% Define quadrants
quad1 = theta_2pi>0  & theta_2pi<pi/2;
quad2 = theta_2pi>pi/2  & theta_2pi<pi;
quad3 = theta_2pi>pi  & theta_2pi<3*pi/2;
quad4 = theta_2pi>3*pi/2  & theta_2pi<2*pi;
 
quad1 = reshape(quad1,[size(quad1,1) size(quad1,2)*size(quad1,3)]); % reshape to 2D 
quad2 = reshape(quad2,[size(quad2,1) size(quad2,2)*size(quad2,3)]);
quad3 = reshape(quad3,[size(quad3,1) size(quad3,2)*size(quad3,3)]);
quad4 = reshape(quad4,[size(quad4,1) size(quad4,2)*size(quad4,3)]);
 
quad_data = cat(3,quad1,quad2,quad3,quad4); % concatenate to make use in loop
 
% loop over z-slices
for i=1:size(Idx,1)
    
    for j=1:length(bin_rad) % loop over bins in radial position
        
        for quadrant = 1:4
            % utheta
            u_theta_bin_mean(i,j,quadrant) = mean(uth_temp(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;        
            u_theta_bin_std(i,j,quadrant)  = std(uth_temp(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;               
            u_theta_bin_cnt(i,j,quadrant)  = length(uth_temp(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;        
 
            % ur
            u_r_bin_mean(i,j,quadrant) = mean(ur_temp(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;                
            u_r_bin_std(i,j,quadrant)  = std(ur_temp(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;               
            u_r_bin_cnt(i,j,quadrant)  = length(ur_temp(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;        
 
            % uz
            u_z_bin_mean(i,j,quadrant) = mean(mean_Uz(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;                
            u_z_bin_std(i,j,quadrant)  = std(mean_Uz(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;                
            u_z_bin_cnt(i,j,quadrant)  = length(mean_Uz(i,ismember(Idx(i,:),j )  & ismember(squeeze(quad_data(i,:,quadrant)),1)   ) ) ;                
        
        end
        
        
    end
    
    
end


r_pres = (posx_new.^2+posy_new.^2).^0.5;
[~,~,Idx_pres] = histcounts(reshape(r_pres,[1 size(r_pres,1)*size(r_pres,2)]),bin_rad);
pres_array = reshape(pres,[1 size(pres,1).*size(pres,2)]);

for j=1:length(bin_rad)
    for quadrant = 1:4
        p_bin_mean(j,quadrant) = nanmean(pres_array(ismember(Idx_pres,j) & ismember(squeeze(quad_data(i,:,quadrant)),1)  ));
        
    end
    
end

figure(69420)
title('M34')
hold all
scatter(bin_rad + 0.5*(max(gradient(bin_rad))), p_bin_mean(:,1))
scatter(bin_rad + 0.5*(max(gradient(bin_rad))), p_bin_mean(:,2))
scatter(bin_rad + 0.5*(max(gradient(bin_rad))), p_bin_mean(:,3))
scatter(bin_rad + 0.5*(max(gradient(bin_rad))), p_bin_mean(:,4))
xlabel('r (mm)','interpreter','latex','fontsize',14)
ylabel('$\overline{p}$ (Pa)','interpreter','latex','fontsize',14)
legend('Quadrant 1','Quadrant 2','Quadrant 3','Quadrant 4','interpreter','latex','fontsize',14,'location','southeast')
ylim([20000 25000])
box on

%% Plotting

slice=1;
for i=1:size(Uy_m,1)
    slice=i;
    figure(11)
    hold all
    contourf(posx_new(1,:),posy_new(:,1),(squeeze(Uy_m(slice,:,:))),-1.2:0.2:1.2);view(2)
    quiver(pos_new(1,:),pos_new(2,:),mean_Ux(slice,:),mean_Uy(slice,:),1,'k');
    xlim([-30 30])
    ylim([-30 30])
    set(gca,'linewidth',2,'FontSize',16)
    colormap(brewermap([],'RdBu'))
%     pause
%     clf
end
figure(29)
imagesc(posx_new(1,:),posy_new(:,1),flipud(squeeze(mean(Nsamples_z,1))));view(2)
xlim([-25 25])
ylim([-25 25])
colormap(brewermap([],'Spectral'))
% caxis([10 100])
set(gca,'linewidth',2,'FontSize',16)%,'colorscale','log')


figure(111);hold all;%scatter(pos_new(1,:),pos_new(2,:),'b+')
quiver(pos_new(1,:),pos_new(2,:),mean_Ux(10,:),mean_Uy(10,:),2,'r');
drawnow
axis equal
%               
figure(333);
subplot(1,3,1)
hold all;%scatter(pos_new(1,:),pos_new(2,:),'b+')
contourf(posx_new(1,:),posy_new(:,1),flipud(squeeze(nanmean(utheta,1))),10)
viscircles([0 0],r_c_m)
drawnow
colormap(brewermap([],'RdYlGn'))
axis equal
axis square
xlim([-25 25])
ylim([-25 25])
box on

subplot(1,3,2);hold all;%scatter(pos_new(1,:),pos_new(2,:),'b+')
contourf(posx_new(1,:),posy_new(:,1),flipud(squeeze(nanmean(ur,1))),10)
viscircles([0 0],r_c_m)
drawnow
colormap(brewermap([],'RdYlGn'))
axis equal
axis square
xlim([-25 25])
ylim([-25 25])
box on

subplot(1,3,3);hold all;%scatter(pos_new(1,:),pos_new(2,:),'b+')
contourf(posx_new(1,:),posy_new(:,1),flipud(squeeze(nanmean(Ux_m,1))),10)
viscircles([0 0],r_c_m)
drawnow
colormap(brewermap([],'RdYlGn'))
axis equal
axis square
xlim([-25 25])
ylim([-25 25])
box on




%% Plotting quadrant wise averaging



for i=1:size(u_theta_bin_mean,1)
    if (unq_zpos(i)>10) && (unq_zpos(i)<55)
        
       
        
        
        
        % Plotting     
        figure(421)
        subplot(2,2,1);hold all
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,1),'rs-')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,1) + (2*u_theta_bin_std(i,:,1)./sqrt(u_theta_bin_cnt(i,:,1))),'k--')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,1) - (2*u_theta_bin_std(i,:,1)./sqrt(u_theta_bin_cnt(i,:,1))),'k--')
%         axis square
%         ylim([0 20])
%         zlim([0 1.2])
        box on
        view([-1 1 1])
        subplot(2,2,2);hold all
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,2),'bs-')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,2) + (2*u_theta_bin_std(i,:,2)./sqrt(u_theta_bin_cnt(i,:,2))),'k--')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,2) - (2*u_theta_bin_std(i,:,2)./sqrt(u_theta_bin_cnt(i,:,2))),'k--')
%         axis square
%         ylim([3 20])
        box on
        view([-1 1 1])
        subplot(2,2,3);hold all
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,3),'gs-')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,3) + (2*u_theta_bin_std(i,:,3)./sqrt(u_theta_bin_cnt(i,:,3))),'k--')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,3) - (2*u_theta_bin_std(i,:,3)./sqrt(u_theta_bin_cnt(i,:,3))),'k--')
%         axis square
%         ylim([0 20])
%         ylim([0 1.2])
        box on
        view([-1 1 1])
        subplot(2,2,4);hold all
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,4),'ms-')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,4) + (2*u_theta_bin_std(i,:,4)./sqrt(u_theta_bin_cnt(i,:,4))),'k--')
        plot3(bin_rad+(0.5*max(gradient(bin_rad))),ones(size(bin_rad)).*unq_zpos(i),u_theta_bin_mean(i,:,4) - (2*u_theta_bin_std(i,:,4)./sqrt(u_theta_bin_cnt(i,:,4))),'k--')
%         axis square
%         ylim([0 20])
%         ylim([0 1.2])
        box on
        view([-1 1 1])
        drawnow
    
    end
end




%% Vortex model fitting and book keeping


% List models
% Model 1: Choi and Ceccio
uth_cc  = fittype(@(lambda_inf,r_v,gama,x) ((lambda_inf./(2.*pi.*(x - (gama*r_c_m./1000)))).*(1  -  (exp(-1.2526.*((x - (gama*r_c_m./1000))./(r_v - (gama*r_c_m./1000))).^2) ))));
% Model 2: Bosschers analytical
uth_bo  = fittype(@(lambda_inf,r_v,x) ((lambda_inf./(2.*pi.*x)).*(1  -  ((r_v.^2./(r_v.^2 + 1.2526.*(r_c_m./1000).^2)).*exp(-1.2526.*(x.^2 - (r_c_m./1000).^2)./r_v.^2)) )));
% Model 3: Bosschers analytical 2
uth_bo2 = fittype(@(lambda_inf,r_v,beta,x) ((lambda_inf./(2.*pi.*x)).*(1  -  beta*exp(-1.2526.*(x.^2 - (r_c_m./1000).^2)./r_v.^2)) ));

% preallocate residuals
res_cc = zeros(1,0);
res_bo = zeros(1,0);
res_bosem = zeros(1,0);

for station = 1:size(u_theta_bin_mean,1)
    for quadrant=1:4

%         % locate data
%         station = 20;
%         quadrant = 3;
        quad_r=((bin_rad+(0.5*max(gradient(bin_rad))))');
        quad_vel = squeeze(u_theta_bin_mean(station,:,quadrant))';
        quad_std = squeeze(u_theta_bin_std(station,:,quadrant))';
        quad_cnt = squeeze(u_theta_bin_cnt(station,:,quadrant))';
        % initialize params
        init_cc = [0.1 7e-3 0.1];
        init_bo = [0.1 7e-3];
        init_bo2 = [0.1 7e-3 1];
        % fit models
        [cc_fitted,gof] = fit(quad_r(~isnan(quad_vel))./1000,quad_vel(~isnan(quad_vel)),uth_cc,'StartPoint',init_cc)
        [bo_fitted,gof] = fit(quad_r(~isnan(quad_vel))./1000,quad_vel(~isnan(quad_vel)),uth_bo,'StartPoint',init_bo)
        [bo2_fitted,gof] = fit(quad_r(~isnan(quad_vel))./1000,quad_vel(~isnan(quad_vel)),uth_bo2,'StartPoint',init_bo2)


        % fit Bosschers semi-empirical model within 95% confidence bounds
        rcm = r_c_m/1000;
        % generate the data
        data_r = quad_r(~isnan(quad_vel))./1000;
        data_uth = quad_vel(~isnan(quad_vel));
        data_std = quad_std(~isnan(quad_vel)); % standard deviation of the mean considered here only
        n_trials = 100;
        data_uth_mat = [];
        for i=1:n_trials
            data_uth_mat(:,i) = normrnd(data_uth,data_std);
        end
        lsqminfun = @(parameters) minfun_bo3(parameters,data_r,data_uth_mat,rcm);
        init_bo3 = [0.1 100 2 1 0.8 10 0.006 1]; % (lambda_inf,zeta_1,zeta_2,p,alphacap,r_v)
        options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
        out_bo3 = lsqnonlin(lsqminfun,init_bo3,[0.05 50 0.1 0.5 0.5 5 0.002 0.5],[0.2 150 10 1.5 1.5 25 0.007 5],options);




        % Plotting of vortex models against the experimental data
        r_space_temp = linspace(1e-3,25e-3,100);
        lambda_inf= out_bo3(1);
        zeta_1= out_bo3(2);
        zeta_2= out_bo3(3);
        p= out_bo3(4);
        alphacap= out_bo3(5);
        K = out_bo3(6);
        r_v= out_bo3(7);
        B = K*r_v;
        q=2;
        beta = out_bo3(8);
        uth_bo3 = (lambda_inf./(2.*pi.*r_space_temp)) .* (1 - alphacap*exp(-zeta_1 .* ((r_space_temp/B).^p))) .* (1 - beta*exp(-zeta_2 .* ( (r_space_temp.^q - rcm.^q)  /r_v.^2)));
        uth_bo3_temp = (lambda_inf./(2.*pi.*data_r)) .* (1 - alphacap*exp(-zeta_1 .* ((data_r./B).^p))) .* (1 - beta*exp(-zeta_2 .* ( (data_r.^q - rcm.^q)  /r_v.^2)));

        
        if 0
            figure(423)
            hold all
            errorbar(quad_r(~isnan(quad_vel)),quad_vel(~isnan(quad_vel)),(2*quad_std(~isnan(quad_vel))./sqrt(quad_cnt(~isnan(quad_vel)))),'ko','linewidth',1)
            plot(r_space_temp.*1e3,cc_fitted(r_space_temp),'r-','linewidth',1)
            plot(r_space_temp.*1e3,bo_fitted(r_space_temp),'b-','linewidth',1)
    %         plot(r_space_temp.*1e3,bo2_fitted(r_space_temp),'m-','linewidth',1)
            plot(r_space_temp.*1e3,uth_bo3,'g-','linewidth',1)
    %         pause
            drawnow
            clf
        end
        
        % book keeping the residuals of the fits
        % (column) : quadrants
        % (row) : stations
        res_cc      = cat(2,res_cc,(abs(cc_fitted(data_r) - quad_vel(~isnan(quad_vel))))'  );
        res_bo      = cat(2,res_bo,(abs(bo_fitted(data_r) - quad_vel(~isnan(quad_vel))))'  );
        res_bosem   = cat(2,res_bosem,(abs(uth_bo3_temp - quad_vel(~isnan(quad_vel))))'  );
        
        
        % Book-keeping all parameters of the vortex models
        % dimension 1 (row): parameters
        % dimension 2 (column) : quadrants
        % dimension 3 (z) : stations
    
        % choi and ceccio params
        ccparams(1,quadrant,station) =  cc_fitted.lambda_inf;   % circulation strength
        ccparams(2,quadrant,station) =  cc_fitted.r_v;          % viscous core radius
        ccparams(3,quadrant,station) =  cc_fitted.gama;         % stress condition param
    
        % bosschers analytical
        boparams(1,quadrant,station) = bo_fitted.lambda_inf;    % circulation strength
        boparams(2,quadrant,station) = bo_fitted.r_v;           % viscous core radius
        
        % bosschers semi-empirical
        bosemparams(1,quadrant,station) = out_bo3(1);   % circulation strength
        bosemparams(2,quadrant,station) = out_bo3(2);   % zeta_1
        bosemparams(3,quadrant,station) = out_bo3(3);   % zeta_2
        bosemparams(4,quadrant,station) = out_bo3(4);   % p
        bosemparams(5,quadrant,station) = out_bo3(5);   % alphacap
        bosemparams(7,quadrant,station) = out_bo3(7);   % r_v
        bosemparams(6,quadrant,station) = out_bo3(6).*out_bo3(7);   % (B = K*r_v)
        bosemparams(8,quadrant,station) = out_bo3(8);   % beta (fitted)
        % bosschers semi-empirical beta magnitude estimate (as we disregarded it in our fit)
        alphabar = alphacap*exp(-zeta_1 * (rcm/B)^p);
        bosemparams(9,quadrant,station) = ( (2*(1-alphabar) * B^p * r_v^q)    -   (alphabar*zeta_1*p*rcm^p*r_v^q) ) / (( (2*(1-alphabar) * B^p * r_v^q)    -   (alphabar*zeta_1*p*rcm^p*r_v^q) )  + ((1-alphabar) * B^p * q * zeta_2 * rcm^q)    );
        
        
                
    end 
end


%% Save data

clear Index Time GridPosition BinnedVelocity VHData


save([folder date rec vsl 'Postproc3DTracking' vsl 'Centered_eulerian.mat'])




end
