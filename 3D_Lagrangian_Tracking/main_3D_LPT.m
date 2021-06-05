%% Info
%-- main file for "Lagrangian Particle [Object] Tracking"
%   
%   $Author: Koen Muller 30/03/2020$
%   
%   Lagrangian tracking at: [+ pros, - cons]
%       - Large lengthscale
%       - Low camera baseline
%       - [Ultra] Wide angle lenses
%       - Relative high reprojection error (i.e. >= O(1 px) & <= bounding ellipse )
%       - High camera occlusion number (close packing/crossing)
%       - High camera occlusion rate (close packing/crossing)
%       - Varying available views
%       - Large scale occlusions
%       - Anisotropy in particle (fish/ellipsoidal shape)
%       - Flickering particle images (algae data)
%       - Varying intensity projection across multiple views
%       + No object crossing in object space (inphysical)
%       + High temporal resolution
%       + High spatial resolution
%
%   Particular finite size and fish shape anisotropy (ellipse) compromises
%   instateneous triangulation by fundamental non-zero error in object
%   point projection from ellipse midpoint. This error increase with lens
%   angle and fish anisotropy, and varies by camera pose, fish orientation
%   and object position.
%
%   This is resolved by tracking from seeding the partial or compromised
%   reconstruction.
%
%-- settings
%       [folder date rec (cal)]:  locate carried out experiment (calibration)
%       ctrl:               controls to processing
%       prop:               properties carried out experiments
%
%-- data management
%       .mat files  binary struct files: such as prop and ctrl
%       .dat files  binary data files for processing
%       .bin files  binary data size
%       .h5 .csv files  final output for sharing
%   
%-- sequential function scripts in proc_data
%       proc_newframes    Process new frame data
%       proc_feasiblesol  Process best feasible tracks over multiple views
%       proc_physobject   Reconstruct physical object optimal solution
%
%-- Commenting throughout the code includes
%       Functionality scripts in headers
%       Possible implementations to make
%       Possible improvements to make
%       Relevant figures to plot and test
%       Structure of file formats
%   look for them in file headers and code commenting
%
%-- Notes
%       This code is optimized for robustness. It runs on single CPU usage,
%       least RAM, and allows parralazation by running multiple processing
%       simultaneously. Therefore it can seen slow for single usage.
%
%       This code has been developed in MatLab 2019b and is tested for
%       2018b, other versions might run into issues which require 
%       backward compatibility and debugging.
%
%       Apart from the main input this code is compatible with non windows
%       operating systems with respect to file reading. 
%
%       Shell scripts to batch run on clusters by Junaid Mehmood handle the
%       main input externally.
%   
%-- TODO
%       clean programming: use vector functionality in for loops for v(:,..)
%       reindexing st all indexing is without gaps
%       explore paralization by e.g. GPU array implementation, Nvidia..
%       clean plotting functionality, relevant inspection and stats.
%       check correctness all print messages while processing

%% Matlab start
close all
clear all
clc

global folder date rec cal ctrl prop plotting ...
    Kmat Rmat tvec Pmat Dpar Dmap Pmap camprop

%% Additional Paths
direc= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(direc) %  back to original folder
clear direc % remove variable

%% Locations
folder='F:\swaraj_nanda'; % Desktop Location
date='\ParticleTracking'; % Date Folder
rec='\Measurement37'; % Recording
cal='\Calibration_200514_141513'; % Calibration files
%% load calibration and define problem
[Dmap,Pmap,Dpar,Kmat,Rmat,tvec,camprop]=load_Cal;

% compute projection matrix
Pmat=cell(size(Rmat));
for c=1:size(Pmat,2)
    Pmat{c}=krtm2pmat(eye(3),Rmat{c},tvec{c});
end

%% create path for processing function storage
if exist([folder date rec vsl 'Preproc3DTracking'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DTracking'])
end

%% Controls
plotting='off'; % switch of all figures for cluster! ( 'off' | 'on')
if 1
    [prop,ctrl] = proc_set;
else
    load([folder date rec vsl 'Preproc3DTracking' vsl 'prop.mat'])
    load([folder date rec vsl 'Preproc3DTracking' vsl 'ctrl.mat'])
end

%% Get data from images, or generator
if 1
    tic
    
    proc_data
    
    toc
end

%% plotting
if strcmp(plotting,'on')
    
    % load last feasible solution
    Fsol=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Fsol.dat']);
    
    % count occurence
    [~,ia,~]=unique(Fsol([1 5],:)','rows');
    l=unique(Fsol(1,ia));
    L=histcounts(Fsol(1,ia),[l-1/2,max(l)+1/2]);
    figure; histogram(L)
    
    % segment
    l=l(L>=10);
    L=ismember(Fsol(1,:),l);
%     L=ismember(Fsol(1,:),Fsol(1,Fsol(5,:)==100));
    
    % plot
    figure
    hold on
    for n=unique(Fsol(5,:))
%         cla
        N=L& Fsol(5,:)==n;%Fsol(5,:)>=n-7 & Fsol(5,:)<=n+7 ; 
        
        scatter3(Fsol(14,N),Fsol(15,N),Fsol(16,N),[],Fsol(1,N),'.')
        view(2)
        axis equal
        xlim([-5 5])%([-2 2])
        ylim([5 25])%([-2 2])
        zlim([-5 1])%([-1 3])
        colormap lines
        colorbar
        drawnow
        
    end
    
    
    %%
    % load physical link data stored from Fsol
    Plink=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Plink.dat']);
    
    % count occurence
    [~,ia,~]=unique(Plink([1 4],:)','rows');
    l=unique(Plink(1,ia));
    L=histcounts(Plink(1,ia),[l-1/2,max(l)+1/2]);
    
%     figure(2); histogram(L,1:30)
    
    % segment
    l=l(L>=5);
    L=ismember(Plink(1,:),l);
    
%     
%     vidfile = VideoWriter('simult_20tsteps2.mp4','MPEG-4');
%     vidfile.FrameRate = 5; % always write in 10 fps
% %     vidfile.LosslessCompression = 'true';
%     open(vidfile);
    %%
    % plot
    figure(1)
    hold all
    for n=unique(Plink(4,:))
        
        N=L&Plink(4,:)==n; 
        
        scatter3(Plink(13,N),Plink(14,N),Plink(15,N),[],Plink(4,N),'.')
        quiver3(Plink(13,N),Plink(14,N),Plink(15,N),...
            Plink(16,N),Plink(17,N),Plink(18,N),0,'.')
        view(2)
        axis equal
        axis tight
%         xlim([-30 30])%([-2 2])
%         ylim([-30 30])%([-2 2])
%         zlim([-10 40])%([-1 3])
        colormap(flipud(jet))
        
        
        
        % Plot visual hull
        
%         load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(n) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox')
%         [fo,vo] = isosurface(Y_vox,X_vox,Z_vox,Vox_int);
%         [fe,ve,ce] = isocaps(Y_vox,X_vox,Z_vox,Vox_int);
%         hfig = figure(1)
%         p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
%         p1.FaceColor = 'blue';
%         p1.EdgeColor = 'none';
% 
%         p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
%            'FaceVertexCData', ce);
%         p2.FaceColor = 'interp';
%         p2.EdgeColor = 'none';
% 
%         camlight(40,40)                                % create two lights 
%         camlight(-20,-20)
%         lighting gouraud
% %         xlim([-5 9])
% %         xlim([-5 9])
%         zlim([-10 50])
%         axis equal
%         axis tight
% %         set(gca,'xtick',[])
% %         set(gca,'xticklabel',[])
% %         set(gca,'ytick',[])
% %         set(gca,'yticklabel',[])
%         hfig.WindowState = 'maximized';
%         view(2)

%         Frm = getframe(gcf);

%         writeVideo(vidfile, Frm);
%         clf
        
        
        
        
        
        
        
        
        drawnow
        
        
        
%         pause
        
        
        
%         Frm = getframe(gcf);
%         clf
%         writeVideo(vidfile, Frm);
        
        
    end
%     close(vidfile)
    
    %%
    % load physical object data
    Pobj=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat']);
    
    % count occurence
    [~,ia,~]=unique(Pobj([1 2],:)','rows');
    l=unique(Pobj(1,ia));
    L=histcounts(Pobj(1,ia),[l-1/2,max(l)+1/2]);
%     figure; histogram(L)
    
  %%  
    track_length_bins = 1:20;
    disparity_bins = linspace(0,1.5,101);
    bin_centers = disparity_bins(1:end-1) + (disparity_bins(2) - disparity_bins(1))/2;
    bin_counts = zeros(length(track_length_bins),length(bin_centers));    
    
    figure(1)
    hold all
    view([1 1 1])
    title('M37','interpreter','Latex','FontSize',20)    
    
    for i=1:track_length_bins(end)        
        % segment
        i
        l=unique(Pobj(1,ia));
        
        l=l(L==i );
        L2=ismember(Pobj(1,:),l);
        
        
        bin_counts(i,:) = histcounts(Pobj(17,L2),disparity_bins,'normalization','probability');
        
        plot3(ones(size(bin_centers)).*i,bin_centers,bin_counts(i,:),'b','LineWidth',1)
        xticks([1 5 10 15 20 25 30])
        yticks([0 0.1 0.2 0.3 0.4 0.5])
        box on
        xlabel('$N_f$','interpreter','latex','FontSize',18)
        ylabel('$\epsilon_{disp}$','interpreter','latex','FontSize',18)
        zlabel('Prob.','interpreter','latex','FontSize',18)
        set(gca,'linewidth',1)
        drawnow

    end
%%
    % segment
    l=l(L>=5 );
    L=ismember(Pobj(1,:),l);
    
    % plot
    figure(1)
    hold all
    for n=1:50%unique(Pobj(2,:))
        
        N=L&Pobj(2,:)==n; 
        
%         X=reshape(qvec2pnts(Pobj(3:12,N),10),3,[]); %ellipse centers
        [ X_c, X_p, X_ang ] = qvec2eshp( Pobj(3:12,N) );
        
        
%         col=reshape(repmat(Pobj(2,N),10,1),1,[]);
%         scatter3(X(1,:),X(2,:),X(3,:),[],col,'.') 
        scatter3(X_c(1,:),X_c(2,:),X_c(3,:),'.') 
        for i=1:length(X_c(1,:))
            [X_el,Y_el,Z_el] = ellipsoid(X_c(1,i),X_c(2,i),X_c(3,i),0.5*X_p(1,i),0.5*X_p(2,i),0.5*X_p(3,i),11);
            % If you knew the angle the axis was rotated by, you could multiply x, y, and z by the rotation matrix.
            rotm = eula2rotm(X_ang(:,i));
            points = [X_el(:)-X_c(1,i), Y_el(:)-X_c(2,i), Z_el(:)-X_c(3,i)];
            points_new = points * rotm + X_c(:,i)';
            
            
            s = surf(reshape(points_new(:,1),size(X_el)),reshape(points_new(:,2),size(X_el)),...
                reshape(points_new(:,3),size(X_el)),ones(size((X_el))).*n,'edgecolor','none');
            
            drawnow
%             pause
        end
        
        view(2)
        axis equal
%         xlim([-5 5])%([-2 2])
%         ylim([5 25])%([-2 2])
%         zlim([-5 1])%([-1 3])
        drawnow
        
    end
    
end
