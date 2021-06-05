function plot_grafiteaux


%% get globals
global folder date rec prop mset %plotting


%% Load data
% get directory of post_grafiteaux
vec_files = dir([folder date rec vsl 'Vortex3DTracking' vsl 'CavityReference' vsl 'Tstep_*']);
% load positions
% load([vec_files(1).folder vsl vec_files(1).name],'pos_new')

% loop over files
% vecs_store = nan([3 size(pos_new,2) 2]); % initiate variable
for ii=1:length(vec_files)
    ii
    % load file
    data = load([vec_files(ii).folder vsl vec_files(ii).name],'vec_accum_mat_cc','pos_accum_mat_cc');
    
%     % store file
%     vecs_store = permute([permute(vecs_store,[2 3 1])  permute(data.vec_accum_mat_cc,[2 3 1])],[3 1 2]);
%     
    
    r(ii,:,:) = (data.pos_accum_mat_cc(1,:,:).^2+data.pos_accum_mat_cc(2,:,:).^2).^0.5;
    z(ii,:,:) = data.pos_accum_mat_cc(3,:,:);
    u_theta(ii,:,:)= (data.vec_accum_mat_cc(1,:,:).^2+data.vec_accum_mat_cc(2,:,:).^2).^0.5;
    u_z(ii,:,:) = data.vec_accum_mat_cc(3,:,:);
    
end

radius = r(:);
tangvel = u_theta(:);

bin_radius = [2:1.5:25];
[~, binidx] = histcounts(radius, bin_radius );

pos_data = z(:);
utheta_mean = [];
utheta_std = [];
dum = unique(binidx);
dumz = unique(pos_data);
for i=1:length(dum)
    i
    dumI = ismember(binidx,dum(i));
    for j=1:length(dumz)
        j
        
        dumJ = ismember(pos_data(:),dumz(j));
        
        utheta_mean(i,j) = nanmean(tangvel(dumI&dumJ));
    
        utheta_std(i,j) = nanstd(tangvel(dumI&dumJ));
        
    end
    
end

if 1
    for i=1:size(utheta_std,2)
        figure(1)
        hold all
%         errorbar(bin_radius./6.5,utheta_mean(:,i).*0.001./1.35,utheta_std(:,i).*0.001./1.35)
        plot3(dumz(i)*ones(size(bin_radius(2:end))),(bin_radius(1:end-1)+ 0.5*(bin_radius(2)-bin_radius(1))   )./6.5,utheta_mean(2:end,i).*0.001./1.35)
        ylabel('$r/r_v$','interpreter','latex','fontsize',20)
        xlabel('$z$','interpreter','latex','fontsize',20)
        zlabel('$u_\theta/u_\theta(r_v)$','interpreter','latex','fontsize',20)
        xlim([0 30])
%         ylim([0 2])
        set(gca,'linewidth',1)
        drawnow
    end
end

if 1
    
    figure(2)
    hold all
    plot((bin_radius(1:end-1)+ 0.5*(bin_radius(2)-bin_radius(1))   ),mean(utheta_mean(2:end,:)').*0.001)
%         plot3(dumz(i)*ones(size(bin_radius)),bin_radius./6.5,utheta_mean(:,i).*0.001./1.35)
    ylabel('$r/r_v$','interpreter','latex','fontsize',20)
    xlabel('$z$','interpreter','latex','fontsize',20)
    zlabel('$u_\theta/u_\theta(r_v)$','interpreter','latex','fontsize',20)
    xlim([0 30])
    ylim([0 2])
    set(gca,'linewidth',1)
    box on
    drawnow
   
end

% Ready the fitting models
% lamb_oseen = @(r) Lambda_inf/(2*pi*r)*(1 - e 


% 
% f = fit(radius(~isnan(tangvel)),tangvel(~isnan(tangvel)),'poly3');
% 
% % vecs_store = vecs_store.*1e-3;
% 
% 
% 
% % plotting
% if 1
%     test_u = nanmean(vecs_store(1,:,:),3).*1e-3;
%     test_v = nanmean(vecs_store(2,:,:),3).*1e-3;
%     test_w = nanmean(vecs_store(3,:,:),3).*1e-3;
%     
%     
%     utheta = (test_u.^2 + test_v.^2).^0.5;
%     
%     figure(69)
%     plot(r,utheta,'.')
% 
%     figure(1)
%     subplot(1,3,1)
%     imagesc(unique(pos_new(1,:)),unique(pos_new(2,:)),reshape(test_u,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))
%     colorbar
%     colormap('jet')
%     axis square
%     axis tight
%     caxis([-1 1])
% %     surf(reshape(pos_new(1,:),[33 33]),reshape(pos_new(2,:),[33 33]),reshape(test_u,[33 33]))
% %     caxis([-800 800])
% %     xlim([-25 25])
% %     ylim([-25 25])
% 
% %     figure(2)
%     subplot(1,3,2)
%     imagesc(unique(pos_new(1,:)),unique(pos_new(2,:)),reshape(test_v,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))
%     colorbar
%     colormap('jet')
%     axis square
%     axis tight
%     caxis([-1 1])
% %     surf(reshape(test_v,[33 33]))
% %     caxis([-800 800])
% %     xlim([-25 25])
% %     ylim([-25 25])
% 
% %     figure(3)
%     subplot(1,3,3)
%     imagesc(unique(pos_new(1,:)),unique(pos_new(2,:)),reshape(test_w,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))
%     colorbar
%     colormap('jet')
%     axis square
%     axis tight
%     caxis([4.5 5.5])
% %     surf(reshape(test_w,[33 33]))
% %     caxis([4500 5500])
% %     xlim([-25 25])
% %     ylim([-25 25])
% 
%     figure(4)
% %     quiver(pos_new(1,:),pos_new(2,:),nanmean(vecs_filtered(1,:,:),3),nanmean(vecs_filtered(2,:,:),3),0)
%     quiver(pos_new(1,:),pos_new(2,:),nanmean(vecs_store(1,:,:),3),nanmean(vecs_store(2,:,:),3),4)
%     axis square
% %     xlim([-25 25])
% %     ylim([-25 25])
% 
% 
%     figure (5)
%     surf(reshape(test_v,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))
% 
% end

end

