function proc_background
%proc_background Background Image processing, image averaging patching etc.

%% Get global variables
global folder date rec ctrl prop plotting

%% check for available data create tspan and indexing
tspan=ctrl.tspan(1):ctrl.tspan(2); % span to proces

%% Get data from image
disp('process images')

% start index labeling
for c=1:prop.res(4) %loop different cameras
    
    if strcmp(plotting,'on')
        figure
    end
    
    % get existing background image if any
    if exist([folder date rec vsl 'bgr_' num2str(c) '.mat'],'file')~=0
        % get existing estimates
        load([folder date rec vsl 'bgr_' num2str(c) '.mat'],'bgr')
        
        % write variables
        avg=bgr.avg;
        dev=bgr.std;
        
        % write variables
        am1=avg*bgr.count;
        am2=(dev.^2+avg.^2).*bgr.count;
        
        % keep counting
        count=bgr.count;
    else
        % variables
        avg=ones(prop.res(1),prop.res(2)); % mean value
        dev=zeros(prop.res(1),prop.res(2)); % standard deviation
        
        % initiate accumalation
        am1=zeros(prop.res(1),prop.res(2)); % accumulate moment 1
        am2=zeros(prop.res(1),prop.res(2)); % accumulate moment 2
        
        % intiate estimates
        bgr.avg=avg; % fit pln
        bgr.std=dev; % no dev
        
        % intitiate counting
        count=0;
    end
    
    for n=randperm(length(tspan)) % loop time in camera frame
        % time windowing including boundaries
        imd=double(import_frames({folder date rec},prop.ext,tspan(n),c));
        
        % substract bgr
        [ ims , scl ] = imrembgr( imd , avg , 2 , 1 );
        
        % segmentation and compress intensity by automated treshold
        [ ~ , ~ , imc ] = imsegscl( ims , 2 , 1 );
        
        % fill image holes
        int=avg.*scl+imc;
        
        % count
        count=count+1;
        
        % update accumulants
        am1=am1 + int;
        am2=am2 + int.^2;
        
        % compute statistics
        avg=am1/count;
        dev=sqrt(abs(am2/count-avg.^2)); % real() .. abs() num in acc at eps
        
        % plotting
        if strcmp(plotting,'on')
            surf(avg')
            shading flat
            colormap gray
            axis tight
            view(2)
            axis equal
            title('Patched background image','interpreter','latex')
            xlabel('x $[\rm{px}]$','interpreter','latex')
            ylabel('y $[\rm{px}]$','interpreter','latex')
            drawnow
        end
        
        % display message
        disp(['processed camera ',num2str(c),' frame ',num2str(n)])
        
    end % n
    
    % save
    bgr.avg=avg;
    bgr.std=dev;
    bgr.count=count;
    save([folder date rec vsl 'bgr_' num2str(c) '.mat'],'bgr')
    
    % write an tiff 
    imwrite(rot90(bgr.avg/max(bgr.avg(:))),[folder date rec vsl 'bgr_' num2str(c) '.tiff'],'tiff')
    
end

end



        % mask
        % msk=mskp | mskm;
        %regionfill(imd,msk);
        
        
        % perception
%         imf=sqrt(imd);
        
        % noise filter
%         coef=SGfilt_coef(imd,ctrl.fobj,ctrl.xtord,1);
        
        % remove noise eval
%         imd=SGfilt_eval(coef,[0 0 0]);
%         tic 

% 
%         % remove noise
%         seg = imopen(seg,strel('rectangle',[1 2])) & imopen(seg,strel('rectangle',[2 1]));%imopen(bw,strel('diamond',1)); % min nonzeros for 2nd ord polyfitbw & bwmorph(...
%         seg = imclose(seg,strel('rectangle',[1 2])) & imclose(seg,strel('rectangle',[2 1]));%imopen(bw,strel('diamond',1)); % min nonzeros for 2nd ord polyfitbw & bwmorph(...
%         seg = bwareaopen(seg,6); % min nonzeros for 2nd ord polyfit
%         seg = ~bwareaopen(~seg,6); % min nonzeros for 2nd ord polyfit
%         
%         % dilate shape
%         seg=imdilate(seg,ctrl.fobj);
        
        % dont filter noise by plane fit keep image res sharp
%         coef=SGfilt_coef(imd,strel('disk',2,0),[2 2 0],1);
%         imf=SGfilt_eval(coef,[0 0 0]);
        