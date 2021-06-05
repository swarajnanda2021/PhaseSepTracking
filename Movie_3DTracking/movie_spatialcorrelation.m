function movie_spatialcorrelation(dom, slc, vpnt)
%movie_eulerianreference Make movie eulerian reference
%
%   Input
%       dat Particular data to plot
%           - Velocity [to be implemented]
%       dom Type of grid position
%           - correlation
%           - Spectrum
%       slc Type of plot
%           - Volume
%           - Slice_[#]
%       vpnt Viewpoint
%           - CameraMotion
%           - TopView
%           - FrontView
%           - SideView
%
%   Output
%       Generate the movie in subfolder
%
%   TODO: Make slc a name value pair, e.g. why is obvious see inline code

%% get globals
global folder date rec prop mset %plotting

%% get data

% load general data
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Time.dat']);

% load case sensitive data
switch dom
    case 'Correlation'
        Distance=fload([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Distance.dat']);
        Correlation=fload([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Correlation.dat']);
        objdom=[min(Distance,[],2) max(Distance,[],2)];
    case 'Spectrum'
        WaveNumber=fload([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'WaveNumber.dat']);
        Spectrum=fload([folder date rec vsl 'Postproc3DTracking' vsl 'SpatialCorrelation' vsl 'Spectrum.dat']);
        Spectrum=Spectrum(1:6,:)+1i*Spectrum(7:12,:);
        objdom=[min(WaveNumber,[],2) max(WaveNumber,[],2)];
    otherwise
        error('Movie Eulerian Reference: Case not listed')
end

%% create path for processing function storage
if exist([folder date rec vsl 'Movie3DTracking' vsl 'SpatialCorrelation'],'dir')~=7
    mkdir([folder date rec vsl 'Movie3DTracking' vsl 'SpatialCorrelation'])
end

%% start movie
% name
name=[vpnt dom slc];
        
% figure or movie
switch mset.ext
    case {'avi' 'mp4'} % movie
        
        % video object
        if strcmp(mset.ext,'avi')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'SpatialCorrelation' vsl name '.' mset.ext],'Motion JPEG AVI');
        elseif strcmp(mset.ext,'mp4')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'SpatialCorrelation' vsl name '.' mset.ext],'MPEG-4');
        end
        vidObj.Quality = mset.qual;
        vidObj.FrameRate = prop.fps;
        open(vidObj);
        
    case {'-dbmp' '-depsc'} % figure
        
        % folder
        if exist([folder date rec vsl 'Movie3DTracking' vsl 'SpatialCorrelation' vsl name],'dir')~=7
            mkdir([folder date rec vsl 'Movie3DTracking' vsl 'SpatialCorrelation' vsl name])
        end
        
end

%% Initiate figure
h=figure(1);
h.Color=[1 1 1];
h.Position=[50 50 800 600];
colormap parula

%% render movie of 3d tracks

% loop frames
for n=1:50%mset.tspan(1):mset.tspan(2)
    
    % Clear figure
    cla
    
    % Time set
    N=Index(2,:)==n ;
    
    % switch cases
    switch dom
        case 'Correlation'
            
            % distance
            dis=Distance(:,N);
            
            % Slice
            switch slc
                case 'Volume'
                    seg=true(1,size(dis,2));
                case 'XYsliceZp1'
                    seg=dis(3,:)==1;
                case 'XYsliceZ0'
                    seg=dis(3,:)==0;
                case 'XYsliceZm1'
                    seg=dis(3,:)==-1;
            end
            
            % get Correlation function
            cor=Correlation(:,N);
            
            % check data
            if nnz(seg)>0
                
                % segment
                dis=dis(:,seg);
                cor=cor(:,seg);
                
                % color
                col=sum(cor([1 4 6],:),1);
                
                % plotting
                scatter3(dis(1,:),dis(2,:),dis(3,:),500*ones(size(col)),col,'.')
                
                % colorbar
                c=colorbar;
                ylabel(c,[dom ' $ \rm{ [',mset.units,'^2 / s^2] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+3*nanstd(col))])
                
            end
            
        case 'Spectrum'
            wav=WaveNumber(:,N);
            
            % Slice
            switch slc
                case 'Volume'
                    seg=true(1,size(wav,2));
                case 'XYsliceZp1'
                    seg=wav(3,:)==1;
                case 'XYsliceZ0'
                    seg=wav(3,:)==0;
                case 'XYsliceZm1'
                    seg=wav(3,:)==-1;
            end
            
            % get velocity
            spe=Spectrum(:,N);
            
            % check data
            if nnz(seg)>0
                
                % segment
                wav=wav(:,seg);
                spe=spe(:,seg);
                
                % color
                col=sum(abs(spe([1 4 6],:)),1);
                
                % plotting
                scatter3(wav(1,:),wav(2,:),wav(3,:),500*ones(size(col)),col,'.')
                
                % colorbar
                c=colorbar;
                ylabel(c,[dom ' $ \rm{ [',mset.units,'^2 / s^2] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
    end
    
    % axis
    view(2)
    axis tight
    axis equal
    
    % Get viewpoint
    switch vpnt
        case 'TopView'
            view([0 90])
        case 'FrontView'
            view([0 0])
        case 'SideView'
            view([90 0])
        case 'CameraMotion'
            view([0 90*(2-erfc( 5*(n-(range(mset.tspan)+1)/2)/((range(mset.tspan)+1)/2)))/2])
    end
    
    % label title and drawing
    xlim(objdom(1,:))
    ylim(objdom(2,:))
    zlim(objdom(3,:))
    grid minor
    box on
    camproj('perspective')
    
    % title
    title([ vpnt ' ' dom ' ' slc ' Frame ' num2str(n) ],'interpreter','latex')
    switch dom
        case 'Correlation'
            xlabel(['$\Delta X \;[\rm{' mset.units '}]$'],'interpreter','latex','FontSize',14)
            ylabel(['$\Delta Y \;[\rm{' mset.units '}]$'],'interpreter','latex','FontSize',14)
            zlabel(['$\Delta Z \;[\rm{' mset.units '}]$'],'interpreter','latex','FontSize',14)
            
        case 'Spectrum'
            xlabel(['$k_x \;[\rm{' mset.units '^{-1}}]$'],'interpreter','latex','FontSize',14)
            ylabel(['$k_y \;[\rm{' mset.units '^{-1}}]$'],'interpreter','latex','FontSize',14)
            zlabel(['$k_z \;[\rm{' mset.units '^{-1}}]$'],'interpreter','latex','FontSize',14)
    end

    % print
    drawnow
    
    % write
    switch mset.ext
        case {'avi' 'mp4'} % movie
            
            % getframe for movie
            frm=getframe(h);
            
            % write video
            writeVideo(vidObj, frm);%
            
        case {'-dbmp' '-depsc'} % figure
            
            % print
            print(h,[folder date rec vsl 'Movie3DTracking' vsl 'SpatialCorrelation' vsl name vsl 'frame_' num2str(n)],mset.rend,mset.qual,mset.ext)
            
    end
    
    % display message
    disp(['Processed ' vpnt ' ' dom ' ' slc ' Frame ' num2str(n)])
    
end

%% close video
close(vidObj)

end

