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

% Load data
BinnedVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedVelocity.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Vortex3DTracking' vsl 'VortexReference'],'dir')~=7
    mkdir([folder date rec vsl 'Vortex3DTracking' vsl 'VortexReference'])
end

if exist([folder date rec vsl 'Vortex3DTracking' vsl 'CavityReference'],'dir')~=7
    mkdir([folder date rec vsl 'Vortex3DTracking' vsl 'CavityReference'])
end




% loop frames
for n=post.tproc(1):post.tproc(2)%mset.tspan(1):mset.tspan(2)
    disp(['Starting tstep ' num2str(n)])
    % Time set
    N=Index(2,:)==n ;
    
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
            centre_x = interp1(Centroid_cav(3,:),Centroid_cav(1,:),unq_zpos(iii));
            centre_y = interp1(Centroid_cav(3,:),Centroid_cav(2,:),unq_zpos(iii));

            slc_pos(1,:) = pos(1,slc_idx);
            slc_pos(2,:) = pos(2,slc_idx);
            vel_pos(1,:) = vel(1,slc_idx);
            vel_pos(2,:) = vel(2,slc_idx);
            vel_pos(3,:) = vel(3,slc_idx);

            slc_pos(1,:) = slc_pos(1,:)-centre_y;
            slc_pos(2,:) = slc_pos(2,:)-centre_x;



            vec_accum_mat_cc(1,:,iii) = vel_pos(1,:);
            vec_accum_mat_cc(2,:,iii) = vel_pos(2,:);
            vec_accum_mat_cc(3,:,iii) = vel_pos(3,:);

            pos_accum_mat_cc(1,:,iii) = slc_pos(1,:);
            pos_accum_mat_cc(2,:,iii) = slc_pos(2,:);
            pos_accum_mat_cc(3,:,iii) = ones(size(slc_pos(1,:))).*unq_zpos(iii);

            if 0
                figure(1);hold all;%scatter(slc_pos(1,:),slc_pos(2,:))
                quiver(slc_pos(1,:),slc_pos(2,:),vel(1,slc_idx),vel(2,slc_idx),'r');
                drawnow
            end
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
    if 1
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





end
