function postproc_graftieaux_stats


%% get globals
global folder date rec prop mset post%plotting


filenames = dir([folder date rec vsl 'Vortex3DTracking' vsl 'VortexReference' vsl 'Tstep_*']);

if length(filenames) > 1 % checks if data exists in the folder, then averages over all within the folder

    for n=1:length(filenames)
        
        %% Load
        data = load([filenames(n).folder vsl filenames(n).name]);
        
        
        %% Convert to cylindrical
        
        % r and theta (should be same for all timesteps
        r = (data.pos_new(1,:).^2 + data.pos_new(2,:).^2).^0.5;
        
        
%         scatter(data.pos_new(1,:),data.pos_new(2,:),1,r,'LineWidth',5)
%         drawnow
%         pause
        
        
%         theta = atan2(data.pos_new(2,:),data.pos_new(1,:));
        theta = cart2pol(data.pos_new(1,:),data.pos_new(2,:));
        
        
%         scatter(r,nanmean(0.001*data.vec_accum_mat(1,:,:),3))
%         drawnow
%         pause
        
%         scatter(data.pos_new(1,:),data.pos_new(2,:),1,nanmean(0.001*data.vec_accum_mat(2,:,:),3),'LineWidth',15)
%         colorbar
%         caxis([-1 1])
%         drawnow
%         pause
%       
%         figure()
%         quiver(data.pos_new(1,:),data.pos_new(2,:),data.vec_accum_mat(1,:,10),data.vec_accum_mat(2,:,10))
%         drawnow
%         pause
        
        % ur and utheta (makes rows the slice entries of each vel. comp.)
        for i=1:size(data.vec_accum_mat,3)
            ur(i,:) = squeeze(0.001*data.vec_accum_mat(1,:,i)).*cos(theta) + squeeze(0.001*data.vec_accum_mat(2,:,i)).*sin(theta);
        
            utheta(i,:) = -squeeze(0.001*data.vec_accum_mat(1,:,i)).*sin(theta) + squeeze(0.001*data.vec_accum_mat(2,:,i)).*cos(theta);
        end
        
%         size(data.pos_new)
%         size(ur)
        
%         scatter(data.pos_new(1,:),data.pos_new(2,:),1,nanmean(utheta,1),'o','LineWidth',10)
%         colorbar
%         caxis([0 7])
%         drawnow
%         pause
%         
        
%         size(utheta)
%         size(r)
%         
        figure(1)
        hold all
        ylim([0 2])
        scatter(r,nanmean(utheta,1),'.')
% %         for kk=1:size(utheta,1)
% %             scatter(r,utheta(kk,:),'.')
% %         end
        drawnow
%         pause
% %         
        
        % collect histogram data for radial position
%         [N,edges,bin] = histcounts(r,30); % bin r and index each bin position
        
        
%         % testing: calculate mean of each bin for u_r and u_theta
%         for i=1:size(utheta,1)
%             for j=1:max(bin)
%             
%             
%             
%                 bin_u_theta(j) = nanmean(utheta(i,find(bin==j)));
%         
%             
%             
%             end
%         end
% %         bin_u_theta
        
%         size(bin_u_theta)
%         size(N)
%         figure(1)
%         hold all
%         ylim([0 1])
%         plot(edges(2:end),bin_u_theta)
%         drawnow
%         pause
%         
        
%         figure(1)
%         hold all
%         
%         scatter(reshp_r,reshp_utheta,'r.')
%         pause
        
        
        
        
        
        
        
        % clear data variable
        data = [];
        
    end




end