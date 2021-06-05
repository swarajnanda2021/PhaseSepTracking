function plot_CavCondAvg


%% get globals
global folder date rec prop mset %plotting


%% Load data
% get directory of post_grafiteaux
vec_files = dir([folder date rec vsl 'Vortex3DTracking' vsl 'CavityReference' vsl 'Tstep_*']);
% load positions
load([vec_files(1).folder vsl vec_files(1).name],'pos_new')

% loop over files
vecs_store = nan([3 size(pos_new,2) 2]); % initiate variable
for ii=1:length(vec_files)
    % load file
    load([vec_files(ii).folder vsl vec_files(ii).name],'vec_accum_mat')
    
    % store file
    vecs_store = permute([permute(vecs_store,[2 3 1])  permute(vec_accum_mat,[2 3 1])],[3 1 2]);
    
    
end

% vecs_store = vecs_store.*1e-3;

%% remove outliers
vecs_filtered = nan(size(vecs_store));
vecs_filtered = permute(vecs_filtered,[3 1 2]);

for ii=1:size(vecs_store,2)
    
    u_vel = squeeze(vecs_store(1,ii,:));
    v_vel = squeeze(vecs_store(2,ii,:));
    w_vel = squeeze(vecs_store(3,ii,:));
    
    if 1 % outlier removal
        outlieridx_u = isoutlier(squeeze(u_vel));
        outlieridx_v = isoutlier(squeeze(v_vel));
        outlieridx_w = isoutlier(squeeze(w_vel));

        outlieridx = (outlieridx_v+outlieridx_u+outlieridx_w);
        u_vel(outlieridx>=1) = nan;
        v_vel(outlieridx>=1) = nan;
        w_vel(outlieridx>=1) = nan;
    end
    
    vecs_filtered(:,:,ii) = [u_vel v_vel w_vel];
    
end

vecs_filtered = permute(vecs_filtered,[2 3 1]);


% plotting
if 1
    test_u = nanmean(vecs_store(1,:,:),3).*1e-3;
    test_v = nanmean(vecs_store(2,:,:),3).*1e-3;
    test_w = nanmean(vecs_store(3,:,:),3).*1e-3;

    figure(1)
    subplot(1,3,1)
    imagesc(unique(pos_new(1,:)),unique(pos_new(2,:)),reshape(test_u,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))
    colorbar
    colormap('jet')
    axis square
    axis tight
    caxis([-1 1])
%     surf(reshape(pos_new(1,:),[33 33]),reshape(pos_new(2,:),[33 33]),reshape(test_u,[33 33]))
%     caxis([-800 800])
%     xlim([-25 25])
%     ylim([-25 25])

%     figure(2)
    subplot(1,3,2)
    imagesc(unique(pos_new(1,:)),unique(pos_new(2,:)),reshape(test_v,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))
    colorbar
    colormap('jet')
    axis square
    axis tight
    caxis([-1 1])
%     surf(reshape(test_v,[33 33]))
%     caxis([-800 800])
%     xlim([-25 25])
%     ylim([-25 25])

%     figure(3)
    subplot(1,3,3)
    imagesc(unique(pos_new(1,:)),unique(pos_new(2,:)),reshape(test_w,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))
    colorbar
    colormap('jet')
    axis square
    axis tight
    caxis([4.5 5.5])
%     surf(reshape(test_w,[33 33]))
%     caxis([4500 5500])
%     xlim([-25 25])
%     ylim([-25 25])

    figure(4)
%     quiver(pos_new(1,:),pos_new(2,:),nanmean(vecs_filtered(1,:,:),3),nanmean(vecs_filtered(2,:,:),3),0)
    quiver(pos_new(1,:),pos_new(2,:),nanmean(vecs_store(1,:,:),3),nanmean(vecs_store(2,:,:),3),4)
    axis square
%     xlim([-25 25])
%     ylim([-25 25])


    figure (5)
    surf(reshape(test_v,[sqrt(size(test_u,2)) sqrt(size(test_u,2))]))

end

end

