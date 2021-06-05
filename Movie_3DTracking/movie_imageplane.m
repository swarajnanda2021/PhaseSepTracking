function movie_imageplane(name)
%process_images Processing images in query
%   here for image exploration

%% Get global variables
global folder date rec cal ctrl prop mset % plotting

%% create path for processing function storage
if exist([folder date rec vsl 'Movie3DTracking' vsl 'ImagePlane'],'dir')~=7
    mkdir([folder date rec vsl 'Movie3DTracking' vsl 'ImagePlane'])
end

%% load calibration and define problem
Rmat=importdata([folder date cal vsl 'Rmat.mat']);
tvec=importdata([folder date cal vsl 'tvec.mat']);
Pmat=cell(size(Rmat));
for c=1:size(Pmat,2)
    Pmat{c}=krtm2pmat(eye(3),Rmat{c},tvec{c});
end

%% load processed data
switch name
    case {'RawImage' 'DewarpedImage'}
        
        % no data to select
    case {'ObjectAndImageEllipsoids'}
        % 2D image data
        Cpln=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Cpln.dat']);
        % 3D object data
        Accuracy=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accuracy.dat']);
        Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
        Ttime=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat']);
        Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat']);
        Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
        Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
        
        
        Plink=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Plink.dat']);
        
        
    case {'ImageEllipses' 'ImageVelocimetry' 'ImageSingleEllipse'}
        
        % 2D image data
        Cpln=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Cpln.dat']);
        
    case {'ObjectEllipsoids' 'ObjectEllipsoidsHalfROI' 'ObjectVelocimetry' 'ObjectSingleEllipsoid'}
        
        % 3D object data
        Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
        Ttime=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat']);
        Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat']);
        Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
        Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
        
        Plink=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Plink.dat']);
        
    otherwise
        
        % message
        error('Processing not found')
        
end
    
%% Get data from image
disp('process images')

% start index labeling
h=figure(1);

% h.Color=[0 0 0];%
prop.fps = 10;
%[1 1 1];%
% inverted colors...

switch name
    case {'ImageSingleEllipse' 'ObjectSingleEllipsoid'}
        % bigger window for 2 figures
        h.Position=[50 50 1150 600];
        
    otherwise
        % size window
%         h.Position=[50 50 800 600];
        
end
for c=1:prop.res(4) %loop different cameras
    
    % initiate object to track
    switch name
        case 'ImageSingleEllipse'
            
            % initiate object to track
            ltrack=sort( Cpln(1,Cpln(1, Cpln(3,:)==mset.tspan(1) & Cpln(2,:)==c ) ));
            ltrack=ltrack(ceil(length(ltrack)/2));
            L=Cpln(1,:)==ltrack;
            
            % ellipse conics
            con=Cpln(4:9, L);%C&
            
            [~,ax,~]=cvec2eshp(con); % avoid imagionary
            wtrack=5*max(sqrt(4*prod(ax,1)/pi));
            
        case 'ObjectSingleEllipsoid'
            
            % initiate object to track
            ltrack=sort( Index( 1 , Index(2,:)==mset.tspan(1) ) );
            ltrack=ltrack(randi(length(ltrack)));
            L=Index(1,:)==ltrack;
            
            % ellipse conics
            quad=Quadric(:,L);
            
            % project
            con=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
                num2cell(quad,1),'UniformOutput',false));
            
            [~,ax,~]=cvec2eshp(con); % avoid imagionary
            wtrack=5*max(sqrt(4*prod(ax,1)/pi));
            
    end
    
    %%%%%%%%%%%%%%%%%%%
    %%% start movie %%%
    %%%%%%%%%%%%%%%%%%%
        
    % figure or movie
    switch mset.ext
        case {'avi' 'mp4'} % movie
            
            % video object
            if strcmp(mset.ext,'avi')
                vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'ImagePlane' vsl name '_Cam' num2str(c) '.' mset.ext],'Uncompressed AVI');
            elseif strcmp(mset.ext,'mp4')
                vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'ImagePlane' vsl name '_Cam' num2str(c) '.' mset.ext],'MPEG-4');
                vidObj.Quality = mset.qual;
            end
            
            vidObj.FrameRate = prop.fps;
            open(vidObj);
            
        case {'-dbmp' '-depsc'} % figure
            
            % folder
            if exist([folder date rec vsl 'Movie3DTracking' vsl 'ImagePlane' vsl name vsl '_Cam' num2str(c) ],'dir')~=7
                mkdir([folder date rec vsl 'Movie3DTracking' vsl 'ImagePlane' vsl name vsl '_Cam' num2str(c) ])
            end
            
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% load camera warpings %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load image warping
    if exist([folder date cal vsl 'wrp_' num2str(c) '.mat'],'file')~=0
        load([folder date cal vsl 'wrp_' num2str(c) '.mat'],'wrp')
    else
        [wrp.y,wrp.x]=meshgrid(1:prop.res(2),1:prop.res(1));
        wrp.X=wrp.x;
        wrp.Y=wrp.y;
        wrp.XX=wrp.X;
        wrp.YY=wrp.Y;
        wrp.xx=wrp.XX;
        wrp.yy=wrp.YY;
        wrp.J=ones(prop.res(1:2));
    end
    
    %%%%%%%%%%%%%%%%
    %%% crop roi %%%
    %%%%%%%%%%%%%%%%
    
    % get roi
    roi=prop.roi;
    
    % background avg std and image grids
    [wrp.J]=imcroproi(wrp.J,roi);
    [wrp.x,wrp.y,wrp.X,wrp.Y]=croproigrid(wrp.x,wrp.y,wrp.X,wrp.Y,roi);
    
    % dewarped grid for interpolation
    droi=[find(wrp.XX(:,1)>=min(wrp.X(:)),1,'first') ...
        find(wrp.YY(1,:)>=min(wrp.Y(:)),1,'first') ...
        find(wrp.XX(:,1)<=max(wrp.X(:)),1,'last') ...
        find(wrp.YY(1,:)<=max(wrp.Y(:)),1,'last')];
    [wrp.XX,wrp.YY,wrp.xx,wrp.yy]=croproigrid(wrp.XX,wrp.YY,wrp.xx,wrp.yy,droi);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% load background image %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load image warping
    if exist([folder date rec vsl 'bgr_' num2str(c) '.mat'],'file')~=0
        load([folder date rec vsl 'bgr_' num2str(c) '.mat'],'bgr')
    else
        bgr.avg=ones(size(wrp.x));
        bgr.std=zeros(size(wrp.x));
        bgr.count=0;
    end
    
    % intensity correction
    bgr.avg=bgr.avg./wrp.J;
    bgr.std=bgr.std./wrp.J;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% crop camera roi %%% % only crop dewarped image
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % get roi
    roi=ctrl.croi(c,:);
    
    % background avg std and image grids
    [~,~,dumX,dumY]=croproigrid(wrp.x,wrp.y,wrp.X,wrp.Y,roi);
    
    % dewarped grid for interpolation
    droi=[min(dumX(:)) ...
        min(dumY(:)) ...
        max(dumX(:))-min(dumX(:)) ...
        max(dumY(:))-min(dumY(:))];
    
    %%%%%%%%%%%%%%%%%%
    %%% downsample %%%
    %%%%%%%%%%%%%%%%%%
    
    % downsample warping
    [wrp.x,wrp.y,wrp.X,wrp.Y] = boxgrid( wrp.x,wrp.y,wrp.X,wrp.Y ,ctrl.rbox);
    [wrp.XX,wrp.YY,wrp.xx,wrp.yy] = boxgrid( wrp.XX,wrp.YY,wrp.xx,wrp.yy ,ctrl.rbox);
    
    % downsample background
    bgr.avg = imbox( bgr.avg , ctrl.rbox );
    bgr.std = imbox( bgr.std , ctrl.rbox );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% loop over all time span %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % brightness perception
    switch mset.cmod
        case 'sqrt'
            bgr.avg=sqrt(bgr.avg);
            bgr.std=sqrt(bgr.std);
        case 'log'
            bgr.avg=log(bgr.avg+1);
            bgr.std=log(bgr.std+1);
    end
    
    for n=mset.tspan(1):mset.tspan(2) % loop time in camera frame
        % clean figure
        cla
        
        % start image data frame query
        imd=double(import_frames({folder date rec},prop.ext,n,c));
        
        % intensity correction
        imf=imd./wrp.J; % image frm
        
        % downsample image
        imf = imbox( imf , ctrl.rbox );
        
        % brightness perspection
        switch mset.cmod
            case 'sqrt'
                imf=sqrt(imf);
            case 'log'
                imf=log(imf+1);
            case 'minmax'
                imf=imminmax(imf,strel('disk',ctrl.fwin(end),0),0.1);
        end
        
        % background
        if 0
            % substract background
            ims = imrembgr( imf , bgr.avg , [2 2 0] );
            
            % segment image
%             [segp,segm]=imsegscl(ims,0,3,5);
            
            % filtering etc..
        end
        
        % image to plot
        plt=imf;
        
        % case by name
        switch name % differs per case
            
            case 'RawImage'
                
                % plot dewarped grid
                surf(wrp.x',wrp.y',0*plt',plt')
                
                % layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
                colormap gray
                xlim([min(wrp.x(:)) max(wrp.x(:))])
                ylim([min(wrp.x(:)) max(wrp.x(:))])
                
            case 'DewarpedImage'
                
                % plot dewarped grid
                surf(wrp.X',wrp.Y',0*plt',plt')
                
                % make up
                view(2)
                shading interp
                %%axis off
                axis equal
                caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+5*std(plt(:))])
                
                % layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
                colormap gray
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                
            case 'ImageEllipses'
                
                % plot dewarped grid
                surf(wrp.X',wrp.Y',0*plt',plt')
                
                % initiate layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-2*std(plt(:)) mean(plt(:))+3*std(plt(:))])
                
                % data to overlay on image
                hold on % append
                
                % plot region of interest
                rectangle('Position',droi,'EdgeColor',[1 0 0],'Linewidth',1)
                
                % select
                C=Cpln(2,:)==c;
                N=Cpln(3,:)==n;
                
                % ellipse conics
                con=Cpln(4:9,C&N);
                
                % plot ellipses
                Xe=cvec2pnts(-con); % avoid imagionary
                plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                    [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'Color',[0 1 0],'LineWidth',1)
                
                % finish layout
                colormap gray
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                
            case 'ImageVelocimetry'
                
                % plot dewarped grid
                surf(wrp.X',wrp.Y',0*plt',plt')
                
                % initiate layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
                
                % data to overlay on image
                hold on % append
                
                % plot region of interest
                rectangle('Position',droi,'EdgeColor',[1 0 0],'Linewidth',1)
                                
                % select
                C=Cpln(2,:)==c;
                N=Cpln(3,:)>=n-floor(mset.tleng/2) & Cpln(3,:)<=n+floor(mset.tleng/2);
                
                % position
                pos=cvec2eshp(Cpln(4:9,C&N));
                
                % velocity
                vel=Cpln(10:11,C&N);
                
                % colormap
                cmap=colormap('parula');
                
                % color coding
                col=sqrt(sum(vel.^2,1));
                col(col>(mean(col)+3*std(col)))=(mean(col)+std(col));
                col=floor((size(cmap,1)-1)*col/max(col))+1;
                
                % color code plotting
                scatter(pos(1,:)',pos(2,:)',[],cmap(col',:),'.')
                
                % finish layout
                colormap gray
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                 
            case 'ImageSingleEllipse'
                
                % select
                C=Cpln(2,:)==c;
                N=Cpln(3,:)==n;
                
                % is tracked object still there?
                if nnz(L&C&N)==0
                    
                    % initiate object to track
                    ltrack=sort( Cpln(1, N & C ) );
                    ltrack=ltrack(ceil(length(ltrack)/2));
                    
                    L=Cpln(1,:)==ltrack;
                    
                    % ellipse conics
                    con=Cpln(4:9,C& L);
                    
                    [~,ax,~]=cvec2eshp(con); % avoid imagionary
                    wtrack=5*max(sqrt(4*prod(ax,1)/pi));
                end
                
                % generate images
                for i=1:2
                    % subplot and clear axis
                    subplot(1,2,i)
                    cla
                    
                    % plot dewarped grid
                    surf(wrp.X',wrp.Y',0*plt',plt')
                    
                    % initiate layout
                    view(2)
                    shading interp
                    %axis off
                    axis equal
                    caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
                    
                    % data to overlay on image
                    hold on % append
                    
                    % ellipse conics
                    con=Cpln(4:9,C&N & ~L);
                    
                    % plot ellipses
                    Xe=cvec2pnts(-con); % avoid imagionary
                    plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                        [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'Color',[0 1 0],'LineWidth',1)
                    
                    % ellipse conics
                    con=Cpln(4:9,C&N & L);
                    
                    % plot ellipses
                    Xe=cvec2pnts(-con); % avoid imagionary
                    plot([Xe(1,:) Xe(1,1)],[Xe(2,:) Xe(2,1)],'Color',[1 0 0],'LineWidth',1)
                    
                    % bounding box
                    pos=cvec2eshp(con); % avoid imagionary
                    bbox=[pos'-wtrack 2*wtrack 2*wtrack];
                    
                    % plot region of interest
                    subplot(1,2,i)
                    rectangle('Position',bbox,'EdgeColor',[1 0 0],'Linewidth',1)
                    
                end
                
                % finish axis 
                subplot(1,2,1)
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                hold off
                
                % finish axis 
                subplot(1,2,2)
                xlim([bbox(1) bbox(1)+bbox(3)])
                ylim([bbox(2) bbox(2)+bbox(4)])
                hold off
                
                % finish layout
                colormap gray
                
            case 'ObjectEllipsoids' 
                
                % plot dewarped grid
                surf(wrp.X',wrp.Y',0*plt',plt')
                
                % initiate layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-2*std(plt(:)) mean(plt(:))+3*std(plt(:))])
                
                % data to overlay on image
                hold on % append
                
                % plot region of interest
                rectangle('Position',droi,'EdgeColor',[1 0 0],'Linewidth',1)
                
                dumI=ismember(Index(1,:),Plink(1,Plink(3,:)==c & Plink(4,:)==n));
                
                
                
                % select
                N=Index(2,:)==n ;%& dumI;
                
                % Long tracks
                [~,test,~] = unique(Index(1,:),'rows');
                l = unique(Index(1,:));
                L=histcounts(Index(1,:),[l-1/2,max(l)+1/2]);
                l=l(L>=mset.lmin);
                L=ismember(Index(1,:),l);
                
                
                % quadric
                quad=Quadric(:,N&L);
                
                % project
                con=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
                    num2cell(quad,1),'UniformOutput',false));
                
                % plot projected ellipses
                Xe=cvec2pnts(con); % avoid imagionary
                if numel(size(Xe))==3
                    plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                        [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'Color',[0 1 0],'LineWidth',1)
                elseif numel(size(Xe))==2
                    plot([Xe(1,:) Xe(1,1)],...
                        [Xe(2,:) Xe(2,1)],'Color',[0 1 0],'LineWidth',1)
                end
                % add index
%                 pos=cvec2eshp(con);
%                 ind=Index(1,N);
%                 text(pos(1,:)',pos(2,:)',num2str(ind'),'Color','green')
                
                % finish layout
                colormap gray
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
            
                
            case 'ObjectAndImageEllipsoids' 
                
                % plot dewarped grid
                surf(wrp.X',wrp.Y',0*plt',plt')
                
                % initiate layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-2*std(plt(:)) mean(plt(:))+3*std(plt(:))])
                
                % data to overlay on image
                hold on % append
                
                % plot region of interest
                rectangle('Position',droi,'EdgeColor',[1 0 0],'Linewidth',1)
                
                dumI=ismember(Index(1,:),Plink(1,Plink(3,:)==c & Plink(4,:)==n));
                
                % select
                N=Index(2,:)==n & dumI;
               
                % Long tracks
                [~,test,~] = unique(Index(1,:),'rows');
                l = unique(Index(1,:));
                L=histcounts(Index(1,:),[l-1/2,max(l)+1/2]);
                l=l(L>=mset.lmin);
                L=ismember(Index(1,:),l);
                
                % quadric
                quad1=Quadric(:,N&L);
            
                
                % project
                con1=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
                    num2cell(quad1,1),'UniformOutput',false));
%                 con2=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
%                     num2cell(quad2,1),'UniformOutput',false));
%                 con3=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
%                     num2cell(quad3,1),'UniformOutput',false));
%                 con4=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
%                     num2cell(quad4,1),'UniformOutput',false));
                
                
                % plot projected ellipses
                Xe1=cvec2pnts(con1); % avoid imagionary
                if numel(size(Xe1))==3
                    plot([squeeze(Xe1(1,:,:));squeeze(Xe1(1,1,:))'],...
                        [squeeze(Xe1(2,:,:));squeeze(Xe1(2,1,:))'],'Color','g','LineWidth',1.5)
                elseif numel(size(Xe1))==2
                    plot([Xe1(1,:) Xe1(1,1)],...
                        [Xe1(2,:) Xe1(2,1)],'Color','g','LineWidth',1.5)
                end
                
%                 
%                  % plot projected ellipses
%                 Xe2=cvec2pnts(con2); % avoid imagionary
%                 if numel(size(Xe2))==3
%                     plot([squeeze(Xe2(1,:,:));squeeze(Xe2(1,1,:))'],...
%                         [squeeze(Xe2(2,:,:));squeeze(Xe2(2,1,:))'],'Color','g','LineWidth',1.5)
%                 elseif numel(size(Xe1))==2
%                     plot([Xe2(1,:) Xe2(1,1)],...
%                         [Xe2(2,:) Xe2(2,1)],'Color','g','LineWidth',1.5)
%                 end
%                 
%                  % plot projected ellipses
%                 Xe3=cvec2pnts(con3); % avoid imagionary
%                 if numel(size(Xe3))==3
%                     plot([squeeze(Xe3(1,:,:));squeeze(Xe3(1,1,:))'],...
%                         [squeeze(Xe3(2,:,:));squeeze(Xe3(2,1,:))'],'Color','g','LineWidth',1.5)
%                 elseif numel(size(Xe3))==2
%                     plot([Xe3(1,:) Xe3(1,1)],...
%                         [Xe3(2,:) Xe3(2,1)],'Color','g','LineWidth',1.5)
%                 end
%                 
               
                
                
                
                
                % add index
%                 pos=cvec2eshp(con);
%                 ind=Index(1,N);
%                 text(pos(1,:)',pos(2,:)',num2str(ind'),'Color','green')
                

                % select
                C=Cpln(2,:)==c;
                N=Cpln(3,:)==n;
                
                % ellipse conics
                con2=Cpln(4:9,C&N);
                
                % plot ellipses
                Xe2=cvec2pnts(-con2); % avoid imagionary
                plot([squeeze(Xe2(1,:,:));squeeze(Xe2(1,1,:))'],...
                    [squeeze(Xe2(2,:,:));squeeze(Xe2(2,1,:))'],'Color','r','LineWidth',1.5)
                


                % finish layout
                colormap gray
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                
                
                
            case 'ObjectEllipsoidsHalfROI' 
                
                % plot dewarped grid
                surf(wrp.X',wrp.Y',0*plt',plt')
                
                % initiate layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
                
                % data to overlay on image
                hold on % append
                
                % plot region of interest
                rectangle('Position',droi,'EdgeColor',[1 0 0],'Linewidth',1)
                
                % plot half line
                plot(droi(1)+droi(3)/2*[1 1],[min(wrp.Y(:)) max(wrp.Y(:))],'g-','Linewidth',1)
                
                % select
                N=Index(2,:)==n;
                
                % quadric
                quad=Quadric(:,N);
                
                % project
                con=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
                    num2cell(quad,1),'UniformOutput',false));
                
                % segment half-plane
                x=cvec2eshp(con);
                S=x(1,:)>(droi(1)+droi(3)/2);
                con=con(:,S);
                
                % plot projected ellipses
                Xe=cvec2pnts(con); % avoid imagionary
                plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                    [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'Color',[0 1 0],'LineWidth',1)
                
                % finish layout
                colormap gray
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                
            case 'ObjectVelocimetry'
%                 N=Index(2,:)>0 & Index(2,:)<=n;
                % plot dewarped grid
                surf(wrp.X',wrp.Y',0*plt',plt')
                
                % initiate layout
                view(2)
                shading interp
                %axis off
                axis equal
                caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
                
                % data to overlay on image
                hold on % append
                
                % plot region of interest
                rectangle('Position',droi,'EdgeColor',[1 0 0],'Linewidth',1)
                
                % select
                N=Index(2,:)>=n-floor(mset.tleng/2) & Index(2,:)<=n+floor(mset.tleng/2);
                
                % quadric
                pos=Position(:,N);
                vel=Velocity(:,N);
                
                % project
%                 vel=homc2inhc(Pmat{c}*inhc2homc(pos+vel)) - ...
%                     homc2inhc(Pmat{c}*inhc2homc(pos)) ;
                pos=homc2inhc(Pmat{c}*inhc2homc(pos)) ;
                
                % colormap
                cmap=colormap('parula');
                
                % color coding
                col=sqrt(sum(vel.^2,1));
                col(col>(mean(col)+3*std(col)))=(mean(col)+3*std(col));
                col=floor((size(cmap,1)-1)*col/max(col))+1;
                
                % color code plotting
                scatter(pos(1,:)',pos(2,:)',[],cmap(col',:),'.')
                
                % finish layout
                colormap gray
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                
            case 'ObjectSingleEllipsoid'
                
                % select
                N=Index(2,:)==n;
                
                % is tracked object still there?
                if nnz(L&N)==0
                    
                    % initiate object to track
                    ltrack=sort( Index( 1 , N ) );
                    ltrack=ltrack(ceil(length(ltrack)/2));
                    L=Index(1,:)==ltrack;
                    
                    % ellipse conics
                    quad=Quadric(:,L);
                    
                    % project
                    con=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
                        num2cell(quad,1),'UniformOutput',false));
                    
                    [~,ax,~]=cvec2eshp(con); % avoid imagionary
                    wtrack=5*max(sqrt(4*prod(ax,1)/pi));
                    
                end
                
                % generate images
                for i=1:2
                    % subplot and clear axis
                    subplot(1,2,i)
                    cla
                    
                    % plot dewarped grid
                    surf(wrp.X',wrp.Y',0*plt',plt')
                    
                    % initiate layout
                    view(2)
                    shading interp
                    %axis off
                    axis equal
                    caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
                    
                    % data to overlay on image
                    hold on % append
                    
                    % ellipse conics
                    quad=Quadric(:,N & ~L);
                    
                    % project
                    con=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
                        num2cell(quad,1),'UniformOutput',false));
                    
                    % plot ellipses
                    Xe=cvec2pnts(con); % avoid imagionary
                    plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                        [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'Color',[0 1 0],'LineWidth',1)
                    
                    % ellipse conics
                    quad=Quadric(:,N & L);
                    
                    % project
                    con=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
                        num2cell(quad,1),'UniformOutput',false));
                    
                    % plot ellipses
                    Xe=cvec2pnts(con); % avoid imagionary
                    plot([Xe(1,:) Xe(1,1)],[Xe(2,:) Xe(2,1)],'Color',[1 0 0],'LineWidth',1)
                    
                    % bounding box
                    pos=cvec2eshp(con); % avoid imagionary
                    bbox=[pos'-wtrack 2*wtrack 2*wtrack];
                    
                    % plot region of interest
                    rectangle('Position',bbox,'EdgeColor',[1 0 0],'Linewidth',1)
                    
                end
                
                % finish axis 
                subplot(1,2,1)
                xlim([min(wrp.X(:)) max(wrp.X(:))])
                ylim([min(wrp.Y(:)) max(wrp.Y(:))])
                hold off
                
                % finish axis 
                subplot(1,2,2)
                xlim([bbox(1) bbox(1)+bbox(3)])
                ylim([bbox(2) bbox(2)+bbox(4)])
                hold off
                
                % finish layout
                colormap gray
                
        end
%         pause
        % finish
        hold off
        
        % adjust axes depending on frame
        if 1
            if c==1
                xlim([-0.08 -0.04])
                ylim([-0.03 0.01])
            elseif c==2
                xlim([0.02 0.06])
                ylim([-0.025 0.015])
            elseif c==3
                xlim([-0.22 -0.0])
                ylim([-0.04 -0.00])
            elseif c==4
                xlim([-0.05 -0.01])
                ylim([-0.09 -0.05])
            elseif c==5
                xlim([0.015 0.055])
                ylim([-0.025 0.005])
            end
        end
        
        h.WindowState = 'maximized';
        
        % draw frame 
        drawnow
        
        % zoom axis
        % make one reference figure and one following figure
        
        % write
        switch mset.ext
            case {'avi' 'mp4'} % movie
                
                % getframe for movie
                frm=getframe(h);
                
                % write video
                writeVideo(vidObj, frm);%
                
            case {'-dbmp' '-depsc'} % figure
                
                % print
                print(h,[folder date rec vsl 'Movie3DTracking' vsl 'ImagePlane' vsl name vsl '_Cam' num2str(c) vsl 'frame_' num2str(n)],mset.rend,mset.qual,mset.ext)
                
        end
        
        % display message
        disp(['processed ' name ' camera ',num2str(c),' frame ',num2str(n)])
        
        
%         pause
        
        set(gca, 'CameraUpVector', [1, 0, 0]);
        
    end % n
    
    % close object and write to file
    close(vidObj)
    
end

end

            % insert text
%             pos=reshape(max(Xe,[],2),2,[]);
%             text(pos(1,:),pos(2,:),...
%                 cellfun(@(x)num2str(x),num2cell(ind,1),'UniformOutput',false),...
%                 'Color',[0 1 0],'FontSize',4,'HorizontalAlignment','left')
            
%         % 1) image filtering: noise and gradient estimate
%         if 0 % SG filter
%             % filter image
%             coef=SGfilt_coef(plt,ctrl.fobj,ctrl.xtord,1); %
%             
%             % im frm
%             plt=SGfilt_eval(coef,[0 0 0]);
%             
% %             dx=SGfilt_eval(coef,[1 0 0]);
% %             dy=SGfilt_eval(coef,[0 1 0]);
%             
% %             imf=sqrt(dx.^2+dy.^2);
%             
%             % plot this
% %             plt=imf;
%         end
        

%         % 2) background modeling and segmentation
%         if 0
%             % remove background
%             [subp,subm]=imrembgr(imf,bgr.avg,0,1);
%             
%             % combine
%             sub=subp+subm;
%             
%             % binarize
%             bw = sub~=0;
%             
%             % remove noise
%             bw = imopen(bw,strel('rectangle',[1 2])) & imopen(bw,strel('rectangle',[2 1]));%imopen(bw,strel('diamond',1)); % min nonzeros for 2nd ord polyfitbw & bwmorph(...
%             bw = imclose(bw,strel('rectangle',[1 2])) & imclose(bw,strel('rectangle',[2 1]));%imopen(bw,strel('diamond',1)); % min nonzeros for 2nd ord polyfitbw & bwmorph(...
%             bw = bwareaopen(bw,6); % min nonzeros for 2nd ord polyfit
%             bw = ~bwareaopen(~bw,6); % min nonzeros for 2nd ord polyfit
%             
%             % plot this
%             plt=sub.*bw; % | bw etc..
%         end

% if 0
%             plt=zeros(size(imf,2),size(imf,1),3);
%             plt(:,:,1)=bgr.avg'/(2*max(bgr.avg(:)))+B';
%             plt(:,:,2)=(bgr.avg.*(1-B))'/(2*max(bgr.avg(:)))+(sub.*bw.*(1-B))'/(2*max(sub(sub>0)));
%             plt(:,:,3)=(bgr.avg.*(1-B))'/(2*max(bgr.avg(:)));
%             plt(plt>1)=1;
%             plt(plt<0)=0;
%             plt=flipud(plt);
%             plt=im2frame(plt);
%         end
%         
%         % colormap conv.
%         if size(plt,3)==1
%             plt=flipud(1+plt'*63); % plot in extension colormap
%             plt=im2frame(plt,cmap);
%         end
%         