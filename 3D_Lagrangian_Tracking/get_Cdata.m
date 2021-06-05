function [Cpln,frm_new] = get_Cdata(frm_rem,frm_new)
%get_Cdata Get camera data and proces images to identify particles and 
%   objects represented as bounding ellipsoids over multiple scales.
%   
%   Input,
%       frm_rem Frames to remove from window data
%       frm_new Frames to append to window data
%   
%   Output,
%       Cpln Camera plane data
%   
%   Notes at image processing:
%   
%      - there is no definite setting to distinguish between correct
%        and maxima, therefore now only maxima.
%      - There is no definite setting for imsegscl in featpoints, which
%        can greatly affect correct and missing identifycation, for now
%        optimized for speed.
%      - This script can be benchmarked with regard to a deeplearning
%        approach (e.g. CNN RNN RCNN LSTM), and if better make
%        improvements or replace by an external program. The ellipsoid 
%        implementation can then represent a statistical likelyhood for 
%        the position incase the position accuracy goes up.
%    
%    solution to the above can perhaps solve both together
%    
%      + Perhaps looking for maxima and minima that not contain saddles in
%        their body span, vice versa for saddles - unicity
%      + preset ellipse shape in complement to find the ellipse from a
%        circular filters, maximizing the contrast as subloop in size
%      + use a histogram mean and std filter instead of imsegscl, or
%        to make improvements inside insegscl.
%    

%% Get global variables
global folder date rec cal ctrl prop plotting

%% Initiate and load data
if exist([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Cpln.dat'],'file')~=0 % load
    
    % load previous image data
    Cpln=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Cpln.dat']);
    
    % define maximum index
    cmax=max(Cpln(1,:)); % if empty, camera loop empty, thus no problem
    
    % remove processed frame if already passed this step
    frm_new=frm_new(:,~ismember(frm_new,Cpln(3,:)));
    
else % initiate
    
    % start image data
    Cpln=zeros(14,0); % [zero-index camera frame ellipse velocity intensity]
    
    % max index
    cmax=0;
    
end

%% needed image filters
rbox=ctrl.dres;
fobj=cell(1,length(ctrl.fwin));
for k=1:length(ctrl.fwin)
    nhood=repmat(getnhood(strel('disk',floor(ctrl.fwin(k)/2),0)),1,1,1); % window [px px frm]
    fobj{k}=strel('arbitrary',nhood);
end

%% Get data from image
disp('Get data cameras')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extract tracking data from camera views %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop different cameras
for c=1:prop.res(4)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% load camera warpings %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load image warping
    if exist([folder date cal vsl 'wrp_',num2str(c),'.mat'],'file')~=0
        load([folder date cal vsl 'wrp_',num2str(c),'.mat'],'wrp')
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
    if exist([folder date rec vsl 'bgr_',num2str(c),'.mat'],'file')~=0
        load([folder date rec vsl 'bgr_',num2str(c),'.mat'],'bgr')
    else
        bgr.avg=ones(size(wrp.x));
        bgr.std=zeros(size(wrp.x));
        bgr.count=0;
    end
    
    % intensity correction
    bgr.avg=bgr.avg./wrp.J;
    bgr.std=bgr.std./wrp.J;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% crop camera dewarped roi %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get roi
    roi=ctrl.croi(c,:);
    
    % background avg std and image grids
    [~,~,dumX,dumY]=croproigrid(wrp.x,wrp.y,wrp.X,wrp.Y,roi);
    
    % dewarped grid for interpolation
    droi=[find(wrp.XX(:,1)>=min(dumX(:)),1,'first') ...
        find(wrp.YY(1,:)>=min(dumY(:)),1,'first') ...
        find(wrp.XX(:,1)<=max(dumX(:)),1,'last') ...
        find(wrp.YY(1,:)<=max(dumY(:)),1,'last')];
    [wrp.XX,wrp.YY,wrp.xx,wrp.yy]=croproigrid(wrp.XX,wrp.YY,wrp.xx,wrp.yy,droi);
    
    %%%%%%%%%%%%%%%%%%
    %%% downsample %%%
    %%%%%%%%%%%%%%%%%%
    
    % downsample warping
    [wrp.x,wrp.y,wrp.X,wrp.Y] = boxgrid( wrp.x,wrp.y,wrp.X,wrp.Y ,ctrl.rbox);
    [wrp.XX,wrp.YY,wrp.xx,wrp.yy] = boxgrid( wrp.XX,wrp.YY,wrp.xx,wrp.yy ,ctrl.rbox);
    
    % downsample background
    bgr.avg = imbox( bgr.avg , ctrl.rbox );
    bgr.std = imbox( bgr.std , ctrl.rbox );
            
%     %figure
%     figure
%     surf(wrp.x',wrp.y',bgr.avg')
%     view(2)
%     shading flat
%     axis tight
%     axis equal
%     colormap gray
%     hold on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Loop over all time span %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % loop over image query
    for n=frm_new % loop time for image query
        
        %%%%%%%%%%%%%%%%%%%%%
        %%% begin message %%%
        %%%%%%%%%%%%%%%%%%%%%
        
        % display message
        disp(['Process new frames camera ',num2str(c),' frame ',num2str(n)])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Load image data in frame queue %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % Initiate and update image queue
        if n==min(frm_new) % start image queue at new frames
            
            % get image data
            imd=double( import_frames( {folder date rec} , prop.ext , ...
                max( ctrl.tproc(1) , n-ctrl.tres ) : min( n+ctrl.tres , ctrl.tproc(2) ) , c ) ); % within feasible frames span
            
            % intensity correct
            imd=imd./repmat(wrp.J,1,1,size(imd,3));
            
            % downsample and start image frame query
            imf = imbox( imd , ctrl.rbox );
            
            % initial frame of interest in query
            N=size(imf,3)-ctrl.tres;
            
%             %figure
%             figure
%             surf(wrp.x',wrp.y',log(imf(:,:,1)'+1))
%             view(2)
%             shading flat
%             axis tight
%             axis equal
%             colormap gray
%             hold on
            
        elseif n+ctrl.tres<=ctrl.tproc(2) % update query at interior dataset
            
            % get new image data
            imd=double(import_frames({folder date rec},prop.ext,n+ctrl.tres,c));
            
            % intensity correct
            imd=imd./wrp.J;
            
            % downsample and append to image frame query
            imf= cat(3, imf , imbox( imd , ctrl.rbox ) );
            
            % remove frist frame or update frame of interest
            if size(imf,3)==2*ctrl.tres+2
                
                % remove first frame
                imf=imf(:,:,2:end);
                
            else
                
                % update the frame indicator
                N=N+1;
                
            end
            
        else % end of dataset
            
             % get image data
            imd=double( import_frames( {folder date rec} , prop.ext , ...
                n : n+1 , c ) ); % within feasible frames span
            
            % intensity correct
            imd=imd./repmat(wrp.J,1,1,size(imd,3));
            
            % downsample and start image frame query
            imf = imbox( imd , ctrl.rbox );
            
            % initial frame of interest in query
            N=size(imf,3)-ctrl.tres;
            
            % remove first frame queue
%             imf(:,:,1)=[];
%             N=1;
            
        end
        
%         figure; surf(imf(:,:,1)'); shading flat; view(2); axis tight; axis equal
%         figure; surf(imf(:,:,2)'); shading flat; view(2); axis tight; axis equal
%         figure; surf(imf(:,:,3)'); shading flat; view(2); axis tight; axis equal
%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% fit background and dewarp %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % substract background - most stable
        ims = imrembgr( log(imf+1) , log(bgr.avg+1) , [2 2 ctrl.tres] );
        
%         figure; surf(ims(:,:,1)'); shading flat; view(2); axis tight; axis equal
%         figure; surf(ims(:,:,2)'); shading flat; view(2); axis tight; axis equal
%         figure; surf(ims(:,:,3)'); shading flat; view(2); axis tight; axis equal
        
        % dewarp image for processing using the cropped grid of interest
        dew = imdewarp( wrp.x , wrp.y , ims , wrp.xx , wrp.yy , 'cubic' , 'bval'); 
        
%         figure; surf(dew(:,:,1)'); shading flat; view(2); axis tight; axis equal
%         figure; surf(dew(:,:,2)'); shading flat; view(2); axis tight; axis equal
%         figure; surf(dew(:,:,3)'); shading flat; view(2); axis tight; axis equal
        
        Pos_fin=[];
        Con_fin=[];
        Vel_fin=[];
        Int_fin=[];
        
         %% Section written by Swaraj, preshifts image based on correlation peak value (corrects advection component ALONE)

        if N==1 % If frame is first frame in the window
            
            
            % Store first two frames of the window in dummy variable dum
            dum = imminmax(dew,strel('disk',ctrl.fwin(end),0),0.1);

            % Perform correlation to estimate inter frame-displacement
            cor=xcorimp(dum(:,:,1),dum(:,:,2));
            
            mpk=cor==max(cor(:));
            xpk_1=wrp.XX(mpk)-mean(wrp.XX(:));
            ypk_1=wrp.YY(mpk)-mean(wrp.YY(:));

            if 1
                % erase first half to get second peak
                if c==1 %peaks up n down, up being the prominent one
                    cor(1:ceil(size(cor,1)/2),:) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==2 %peaks left and right
                    cor(:,ceil(size(cor,2)/2):end) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==3 %peaks up and down, up is dimmer peak
                    cor(ceil(size(cor,1)/2):end,:) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==4 % left right, left is dimmer peak
                    cor(:,ceil(size(cor,2)/2):end) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==5 % left and right, left is dimmer
                    cor(:,ceil(size(cor,2)/2):end) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                end    
            else
                xpk_2=xpk_1;
                ypk_2=ypk_1;
            end

%             imagesc(cor); colorbar;drawnow

%             mpk=cor==max(cor(:));
            
            xpk = 0.5*(xpk_1+xpk_2);
            ypk = 0.5*(ypk_1+ypk_2);



%             pause

            % Examine correlation to estimate shift contribution from
            % highest peak
%             mpk=cor==max(cor(:));
%             xpk=wrp.XX(mpk)-mean(wrp.XX(:));
%             ypk=wrp.YY(mpk)-mean(wrp.YY(:));

            

            % Make dummy variable shft and perform shift in both x and y
            % (therefore biased towards peak displacement or brightest
            % particles)
            shft=dew;
            shft(:,:,2) = imdewarp( wrp.XX , wrp.YY , shft(:,:,2) , wrp.XX+xpk , wrp.YY+ypk , 'cubic' , 'bval'); 

            dew=shft;
%             
%             figure(1); surf(shft(:,:,1)); shading flat; view(2); axis tight; axis equal;colormap(flipud(gray))
%             figure(2); surf(shft(:,:,2)); shading flat; view(2); axis tight; axis equal;colormap(flipud(gray))
%             drawnow
%             pause

        elseif N==2 && n+ctrl.tres<ctrl.tproc(2) % if n is some intermediate frame number in the frame window

            % Store first two frames of the window in dummy variable dum
            dum = imminmax(dew,strel('disk',ctrl.fwin(end),0),0.1);

            % Perform 2 correlations to estimate inter frame-displacement
            % in the two sets {dum(:,:,1:2)} and {dum(:,:,2:3)}
            cor1=xcorimp(dum(:,:,1),dum(:,:,2));
            cor2=xcorimp(dum(:,:,2),dum(:,:,3));
            
            
%             figure; imagesc(cor2)
%             drawnow
%             pause
            
            
            coravg = 0.5*(cor1+cor2);
            
%             figure;imagesc(coravg);drawnow;pause
            
            mpk=coravg==max(coravg(:));
            xpk_1=wrp.XX(mpk)-mean(wrp.XX(:));
            ypk_1=wrp.YY(mpk)-mean(wrp.YY(:));
            
            if 1 % modify cor only if the shifter_loop variable is 2
                if c==1 %peaks up n down, up being the prominent one
                    coravg(1:ceil(size(coravg,1)/2),:) = 0;
                    mpk=coravg==max(coravg(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==2 %peaks left and right
                    coravg(:,ceil(size(coravg,2)/2):end) = 0;
                    mpk=coravg==max(coravg(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==3 %peaks up and down, up is dimmer peak
                    coravg(ceil(size(coravg,1)/2):end,:) = 0;
                    mpk=coravg==max(coravg(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==4 % left right, left is dimmer peak
                    coravg(:,ceil(size(coravg,2)/2):end) = 0;
                    mpk=coravg==max(coravg(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==5 % left and right, left is dimmer
                    coravg(:,ceil(size(coravg,2)/2):end) = 0;
                    mpk=coravg==max(coravg(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                end  
            else
                xpk_2=xpk_1;
                ypk_2=ypk_1;
            end


%             mpk1=coravg==max(coravg(:));
%             xpkmean=wrp.XX(mpk1)-mean(wrp.XX(:));
%             ypkmean=wrp.YY(mpk1)-mean(wrp.YY(:));


%             xpk=xpkmean;
%             ypk=ypkmean;
%             xpk_1
%             xpk_2

            xpk = 0.5*(xpk_1+xpk_2);
            ypk = 0.5*(ypk_1+ypk_2);
            
            
            
            % Make dummy variable shft and perform shift in both x and y
            % (therefore biased towards peak displacement or brightest
            % particles)
            shft=dew;
            shft(:,:,1) = imdewarp( wrp.XX , wrp.YY , shft(:,:,1) , wrp.XX-xpk , wrp.YY-ypk , 'cubic' , 'bval'); 
            shft(:,:,3) = imdewarp( wrp.XX , wrp.YY , shft(:,:,3) , wrp.XX+xpk , wrp.YY+ypk , 'cubic' , 'bval'); 

            
            
            dew=shft;
            
        elseif N==2 && n+ctrl.tres==ctrl.tproc(2) % end of dataset
            
            % Store first two frames of the window in dummy variable dum
            dum = imminmax(dew,strel('disk',ctrl.fwin(end),0),0.1);

            % Perform correlation to estimate inter frame-displacement
            cor=xcorimp(dum(:,:,1),dum(:,:,2));
            
            mpk=cor==max(cor(:));
            xpk_1=wrp.XX(mpk)-mean(wrp.XX(:));
            ypk_1=wrp.YY(mpk)-mean(wrp.YY(:));

            if 1
                % erase first half to get second peak
                if c==1 %peaks up n down, up being the prominent one
                    cor(1:ceil(size(cor,1)/2),:) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==2 %peaks left and right
                    cor(:,ceil(size(cor,2)/2):end) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==3 %peaks up and down, up is dimmer peak
                    cor(ceil(size(cor,1)/2):end,:) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==4 % left right, left is dimmer peak
                    cor(:,ceil(size(cor,2)/2):end) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                elseif c==5 % left and right, left is dimmer
                    cor(:,ceil(size(cor,2)/2):end) = 0;
                    mpk=cor==max(cor(:));
                    xpk_2=wrp.XX(mpk)-mean(wrp.XX(:));
                    ypk_2=wrp.YY(mpk)-mean(wrp.YY(:));
                end    
            else
                xpk_2=xpk_1;
                ypk_2=ypk_1;
            end

%             imagesc(cor); colorbar;drawnow

%             mpk=cor==max(cor(:));
            
            xpk = 0.5*(xpk_1+xpk_2);
            ypk = 0.5*(ypk_1+ypk_2);



%             pause

            % Examine correlation to estimate shift contribution from
            % highest peak
%             mpk=cor==max(cor(:));
%             xpk=wrp.XX(mpk)-mean(wrp.XX(:));
%             ypk=wrp.YY(mpk)-mean(wrp.YY(:));

            

            % Make dummy variable shft and perform shift in both x and y
            % (therefore biased towards peak displacement or brightest
            % particles)
%             shft=dew;
%             shft(:,:,2) = imdewarp( wrp.XX , wrp.YY , shft(:,:,2) , wrp.XX+xpk , wrp.YY+ypk , 'cubic' , 'bval'); 
% 
%             dew=shft;
            shft=dew;
            shft(:,:,1) = imdewarp( wrp.XX , wrp.YY , shft(:,:,1) , wrp.XX-xpk , wrp.YY-ypk , 'cubic' , 'bval'); 
            shft(:,:,3) = imdewarp( wrp.XX , wrp.YY , shft(:,:,3) , wrp.XX+xpk , wrp.YY+ypk , 'cubic' , 'bval'); 

            
            
            dew=shft;
            
            
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% loop over multiple scales in image to extract data %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Initiate mask and image grid
        msk=false(size(dew)); % initiate mask for removal 
        XX=wrp.XX; 
        YY=wrp.YY;
        Tp=[diff(XX(1:2,1)) 0 XX(1,1)-diff(XX(1:2,1))
            0 diff(YY(1,1:2)) YY(1,1)-diff(YY(1,1:2))
            0 0 1]; % pixel to scale transformation matrix -1 pix shift (T=Tcal*Tpix, did the math)
        
        % Initaite idenitification variables
        Pos=zeros(2,0); % position
        Vel=zeros(2,0); % velocity
        Con=zeros(6,0); % conic
        Int=zeros(3,0); % peak intensity
        Levl=zeros(1,0); % scale level
        Wght=zeros(1,0); % identification weight
        Feas=false(1,0); % final/feasible ellipse       
        
        % maximize local image constrast identifications
        sclr=1;
        while numel(dew)>0 && numel(dew)>numel(fobj{end}.Neighborhood) && sclr<=ctrl.nscl
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% define objective scale and spatial mapping %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % pixel to scale transformation matrix
            Tc=[diff(XX(1:2,1)) 0 XX(1,1)-diff(XX(1:2,1))
                0 diff(YY(1,1:2)) YY(1,1)-diff(YY(1,1:2))
                0 0 1]; % -1 pix shift (T=Tcal*Tpix, did the math)
            
            % loop filters
            for f=1:length(fobj)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% identify ellipses by local max (or min / both) %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % define objects to interest
                [scls,scll]=imrelscl(dew,fobj{f});
                
                % Define object
                obj=(scls - scll); % use abs to include minima max(cat(4,scls,scll),[],4);%
                
%                 figure; surf(obj(:,:,1)'); shading flat; view(2); axis tight; axis equal
%                 figure; surf(obj(:,:,2)'); shading flat; view(2); axis tight; axis equal
%                 figure; surf(obj(:,:,3)'); shading flat; view(2); axis tight; axis equal
%                 
                % detect new features abs=min and max    (:,:,1:2)
                [pos,vel,con,int]=featpoints(obj(:,:,N),fobj{f},'maxima',1); % ,wghtuse ellipse/minima to include minima
                
                % no definite strategy to include minima, which in turn can
                % enchance maxima
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% transform to spatial mapping %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % transform  position
                pos=homc2inhc(Tc*inhc2homc(pos'));
                
                % transform velocity
                vel=Tc(1:2,1:2)*vel';
                
                % transform conics
                con= contrans( con',Tc );
                
                % transform peaks
                int=int';
                
                % weight
                [~,ax,~]=cvec2eshp(con);
                wght=prod(ax,1); %sqrt(sum(ax.^2,1))sqrt(prod(ax,1)); % sqrt(sum(ax.^2,1));  wght'; %
                
                % define feasible set
                feas=true(1,size(int,2));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% inter scale adjacency %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % selection from previous
                sel= Feas & Levl>=( (sclr-1)*length(fobj) + f - 1 ) ;%length(fobj)
                seg= feas ;
                
                % peak contrast
                T=( (1:size(imd,3))' - N ).^(0:2);
                pnp1=mean(T*Int(:,sel),1); % [],abs value for ellipses
                pint=mean(T*int(:,seg),1); % [],
                
                % weight
                wnp1=Wght(1,sel);
                wint=wght(1,seg);
                
                % coherent adjacency in displacement
                AdjI=sparse(nnz(sel),nnz(seg));
                for t=0:ctrl.tres % consistency obj rec ((1:size(dew,3))-N) % 
                    
                    % transform the conic
                    Cdis=contrans(Con(:,sel),t*Vel(:,sel),'displace');
                    
                    % distance on unit ellipse
                    Es=ellipsedist(Cdis,pos(:,seg)+t*vel(:,seg),'matrix');
                    
                    % trasnform the conic
                    cdis=contrans(con(:,seg),t*vel(:,seg),'displace');
                    
                    % distance on unit ellipse
                    El=ellipsedist(cdis,Pos(:,sel)+t*Vel(:,sel),'matrix');
                    
                    % unique back forth within ellipse adjacency
                    AdjI=AdjI + sparse( double( Es<=ctrl.ovlp & El'<=ctrl.ovlp ) );
                    
                end
                AdjI=AdjI==(1+2*ctrl.tres);% basicly & statement w.r.t above
                
                % 0) old
                wrtF=full(sum(AdjI,2)'==0);
                
                % 1) new
                wrtf=full(sum(AdjI,1)==0);
                
                % find adjacency pairs over scale decompose A=I*J'
                [I,J]=find(AdjI);
                selI=sparse(I,1:length(I),ones(size(I)),size(AdjI,1),length(I));
                selJ=sparse(J,1:length(J),ones(size(J)),size(AdjI,2),length(J));
                
                % average/max/min peak value for robust vote
%                 dvot=pnp1(I)>pint(J);
                dvot=pnp1(I)./wnp1(I) > pint(J)./wint(J);
%                 dvot=pnp1(I).*wnp1(I) > pint(J).*wint(J);
                
                % get consistent best matches
                sI=(selI*dvot')>0.5;
                sJ=(selJ*~dvot')>0.5;
                
                % write
                wrtF(sI)=1; % uI(sI)
                wrtf(sJ)=1; % uJ(sJ)
                
                % segmentation
                Feas(sel)=Feas(sel) & wrtF; % demand unique to write
                feas(seg)=feas(seg) & wrtf; % allow new pos to be found
                
                %%%%%%%%%%%%%
                %%% write %%%
                %%%%%%%%%%%%%
                
                % new points at spatial grid
                Pos=cat(2,Pos,pos);
                
                % velocity vectors at spatial grid
                Vel=cat(2,Vel,vel);
                
                % peaks
                Int=cat(2,Int,int);
                
                % adjust conics to scale grid
                Con=cat(2,Con,con);
                
                % level
                Levl=cat(2,Levl,ones(size(seg))*( (sclr-1)*length(fobj) + f )); %
                
                % weight
                Wght=cat(2,Wght,wght); %
                
                % selected feasible
                Feas=cat(2,Feas,feas);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% remove selected data from image %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % remove written from previousSnp1 &
                rem= Feas & Levl==( (sclr-1)*length(fobj) + f - 1 );%- - 
                for k=1:size(obj,3)
                    % displace ellipses
                    crem=contrans(Con(:,rem),...
                        Vel(:,rem).*(k-N),'displace');
                    
%                     crem=contrans(crem,1/sqrt(2),'resize'); % 0.75
                    
                    % transform back to pixels scale
                    crem=contrans(crem,inv(Tp));
                    
                    % mask valid ellipses
                    if ~isempty(crem)
                        prem=cvec2pnts(-crem,10); % high sampling
                        
                        erem=immappol( [size(obj,1) size(obj,2)] , prem )>0;
                        
                        %figure; surf(double(erem)'); shading flat; view(2); axis tight; axis equal; drawnow
                        
                        msk(:,:,k)=msk(:,:,k) | bwmorph(erem,'dilate',1);%erem;%
                        
                        if numel(msk(:,:,k))~=nnz(msk(:,:,k))
                            dew(:,:,k)=regionfill(dew(:,:,k),msk(:,:,k));
                        else
                            dew(:,:,k)=0;
                        end
                        
                    end
                    
                end
                Tp=Tc;
                
%                 figure
%                 surf(wrp.X',wrp.Y',0*wrp.X',log(imf(:,:,N)'+1)); shading flat; view(2);axis tight; axis square; drawnow
%                 Xe=reshape(cvec2pnts(-Con(:,Feas)),2,[]);
%                 hold on; plot(Xe(1,:),Xe(2,:),'g.','MarkerSize',2)
%                 plot(Pos(1,Feas),Pos(2,Feas),'r.','MarkerSize',5); 
%                 quiver(Pos(1,Feas),Pos(2,Feas),Vel(1,Feas),Vel(2,Feas),0,'r','LineWidth',1);hold off
%                 title(['n id = ',num2str(nnz(Feas))])
%                 xlim([min(wrp.X(:)) max(wrp.X(:))])
%                 ylim([min(wrp.Y(:)) max(wrp.Y(:))])% xlim([-0.1 -0.0]);ylim([-0.05 0.05])
%                 drawnow
%                 colormap gray
%                 caxis([4 6])
%                 pause
                
            end % f
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% downgrade image to next spatial scale %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % box the average image to desired scale
            dew=imbox(dew,rbox,'int'); % integrate to continue maximize intensity to lower resolution
            msk=imbox(double(msk),rbox,'avg')>0; % average to lower resolution
            [XX,YY] = boxgrid( XX,YY,[], [] ,rbox);
            
            % raise image scale
            sclr=sclr+1;
            
        end % s
        
        % keep valid
        Pos=Pos(:,Feas);
        Con=Con(:,Feas);
        Vel=Vel(:,Feas)+[xpk ; ypk];
        Int=Int(:,Feas);
        
        if strcmp(ctrl.brightsel,'true')==1 % If you have decided to segment out bubbles only
        % Segment N brightest
            
            [~,index_sort] = sort(Int(1,:));
            
%             size(Int)
%             
%             size(Int,2)-ctrl.nparts+1
            if size(Int,2)>ctrl.nparts(c) % if detections are greater than intended bright particles
                indices_2 = index_sort(max([1 size(Int,2)-ctrl.nparts(c)+1]):end); % track ctrl.nparts brightest only 
            else
                indices_2 = index_sort(max([1 size(Int,2)]):end); % use all particles if less than threshold
            end
            % Segment
            Pos_fin = Pos(:,indices_2);
            Con_fin = Con(:,indices_2);
            Vel_fin = Vel(:,indices_2);
            Int_fin = Int(:,indices_2);
        
        else
        
            % append data
            Pos_fin=[Pos_fin Pos];
            Con_fin=[Con_fin Con];
            Vel_fin=[Vel_fin Vel];
            Int_fin=[Int_fin Int];

        end
        
        %% Local correlation
        
        
        
        nKernel = ctrl.nKernel; % Set interrogation for local correlation estimation, ODD NUMBERS ONLY!! 7 IS TOO SMALL APPARENTLY
        [x,y] = meshgrid(-nKernel:nKernel,-nKernel:nKernel); % Make meshgrid 5X5
        Grd_coord = [x(:),y(:)];

        % Convert physical location of detections in
        % [Pos_fin(1,:)::x,Pos_fin(2,:)::y] into pixel indices in dewarped
        % image coordinate system
        
        
        [~,idxLoc_y] = min(abs(Pos_fin(1,:)-wrp.XX(:,1)));
        [~,idxLoc_x] = min(abs(Pos_fin(2,:)-wrp.YY(1,:)'));
        
     
        % Remove detected particles that are too close to boundary
        RemPart_x = [find((idxLoc_x-(nKernel+1)<=0)) find((idxLoc_x+(nKernel+1)>=size(shft,2)))];
        RemPart_y = [find((idxLoc_y-(nKernel+1)<=0)) find((idxLoc_y+(nKernel+1)>=size(shft,1)))];
        
        Pos_fin_2 = Pos_fin;
        Vel_fin_2 = Vel_fin;
        Con_fin_2 = Con_fin;
        Int_fin_2 = Int_fin;
        idxLoc_x_2 = idxLoc_x;
        idxLoc_y_2 = idxLoc_y;
        
        idxLoc_x_2(:,[RemPart_x RemPart_y]) = [];
        idxLoc_y_2(:,[RemPart_x RemPart_y]) = [];
        Pos_fin_2(:,[RemPart_x  RemPart_y]) = [];
        Con_fin_2(:,[RemPart_x  RemPart_y]) = [];
        Vel_fin_2(:,[RemPart_x  RemPart_y]) = [];
        Int_fin_2(:,[RemPart_x  RemPart_y]) = [];

        % Convert subscript indices to linear indices on dewarped
        % pixel coordinate system

        Locs = [idxLoc_x_2;idxLoc_y_2]';
        
        
        gridPts = bsxfun(@plus,Grd_coord,permute(repmat(Locs,[1,1,size(Grd_coord,1)]),[3 2 1])); % adds the generic box to all our points of interest using element-by-element operations (bsxfun)
        idx = squeeze(sub2ind(size(squeeze(shft(:,:,1))),gridPts(:,2,:),gridPts(:,1,:))); % converts the pixel coordinates from subscript to linear indices
        

        im_new_N = squeeze(shft(:,:,1));
        im_new_Np1 = squeeze(shft(:,:,2));

      
        
        
        % Get pixel positions on indices
        vals_X = wrp.XX(idx);
        vals_Y = wrp.YY(idx);
        vals_J_N = im_new_N(idx);
        vals_J_Np1 = im_new_Np1(idx);
        
        test2_J_N = reshape(vals_J_N,[size(x,1) size(x,2) size(Locs,1)]);
        test2_J_Np1 = reshape(vals_J_Np1,[size(x,1) size(x,2) size(Locs,1)]);
        test2_X = reshape(vals_X,[size(x,1) size(x,2) size(Locs,1)]);
        test2_Y = reshape(vals_Y,[size(x,1) size(x,2) size(Locs,1)]);
        
        % Calculate cross-correlation between both image frames
        R = [];
        disp(['Localized Correlation over ' (2*num2str(nKernel) +1) ' pixels^2 (dewarped)'])
        for i=1:size(test2_J_N,3)
            
%             R(:,:,i) = xcorimp(imminmax(test2_J_N(:,:,i),strel('disk',ctrl.fwin(end),0),0.1) , imminmax(test2_J_Np1(:,:,i),strel('disk',ctrl.fwin(end),0),0.1) );
            R(:,:,i) = xcorimp(test2_J_N(:,:,i) , test2_J_Np1(:,:,i) );
        end
     
        mpk_loc = R == max(R,[],[1,2]); % logical indexing to pick out the maximum peak only. 
        
        
        ypk_loc=test2_X(mpk_loc)-squeeze(mean(test2_X,[1 2]));
        xpk_loc=test2_Y(mpk_loc)-squeeze(mean(test2_Y,[1 2]));
       
        
        % Append vectors to Vel function
        Vel_fin_2 = Vel_fin_2+[ypk_loc' ; xpk_loc'];

        
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%
        %%% plot results %%%
        %%%%%%%%%%%%%%%%%%%%
         %%
        
        
        %%%%%%%%%%%%%%%%%%%%
        %%% plot results %%%
        %%%%%%%%%%%%%%%%%%%%
%         plotting = 'on';
        if strcmp(plotting,'on')
            
            %figure
            figure(1)
            
            % clear axis
            cla
            
            % quickfix beneath
            
            
            
            plt=log(imf(:,:,N)+1);
            
%             imshow(plt);drawnow;pause
            
            surf(wrp.X',wrp.Y',0*plt',plt')
            view(2)
            shading flat
            axis tight
            axis equal
            colormap gray
            hold on
%             caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
%             caxis([4 6])
            caxis([4.5 5.5])
            % ellipses
            Xe=cvec2pnts(-Con_fin_2); % avoid imagionary
            plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-g','LineWidth',1)
            
            % position
            plot(Pos_fin_2(1,:),Pos_fin_2(2,:),'.r','MarkerSize',5)
            
            % velocity
            quiver(Pos_fin_2(1,:),Pos_fin_2(2,:),Vel_fin_2(1,:),Vel_fin_2(2,:),0,'r','LineWidth',1)
            
            % axis
            hold off
            xlabel('x-ref [m]','interpreter','latex')
            ylabel('y-ref [m]','interpreter','latex')
            title(['tracking cam ',num2str(c),' frm ',num2str(n)],'interpreter','latex')
%             xlim([min(wrp.X(:)) max(wrp.X(:))])
%             ylim([min(wrp.Y(:)) max(wrp.Y(:))])
%             xlim([-0.125 -0.075]);ylim([-0.025 0.025])
%             xlim([-0.1 -0.0]);ylim([-0.05 0.05])
%             xlim([-0.1 0.1]); ylim([-0.1 0.1])
%             xlim([-0.08 -0.07]);ylim([0.0 0.01])
            drawnow
%             pause
            
            
            figure(2)
            
            % clear axis
            cla
            
            % quickfix beneath
            
            
            
            plt=log(imf(:,:,N+1)+1);
            
%             imshow(plt);drawnow;pause
            
            surf(wrp.X',wrp.Y',0*plt',plt')
            view(2)
            shading flat
            axis tight
            axis equal
            colormap gray
            hold on
%             caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
%             caxis([4 6])
            caxis([4.5 5.5])
            % ellipses
            Xe=cvec2pnts(-Con_fin_2); % avoid imagionary
            plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-g','LineWidth',1)
            
            % position
            plot(Pos_fin_2(1,:),Pos_fin_2(2,:),'.r','MarkerSize',5)
            
            % velocity
            quiver(Pos_fin_2(1,:),Pos_fin_2(2,:),Vel_fin_2(1,:),Vel_fin_2(2,:),0,'r','LineWidth',1)
            
            % axis
            hold off
            xlabel('x-ref [m]','interpreter','latex')
            ylabel('y-ref [m]','interpreter','latex')
            title(['tracking cam ',num2str(c),' frm ',num2str(n)],'interpreter','latex')
%             xlim([min(wrp.X(:)) max(wrp.X(:))])
%             ylim([min(wrp.Y(:)) max(wrp.Y(:))])
%             xlim([-0.125 -0.075]);ylim([-0.025 0.025])
%             xlim([-0.1 -0.0]);ylim([-0.05 0.05])
%             xlim([-0.1 0.1]); ylim([-0.1 0.1])
%             xlim([-0.08 -0.07]);ylim([0.0 0.01])
            drawnow
            pause
            
            
        end
%         plotting = 'off';
        
        %%%%%%%%%%%%%%%%%%%%%
        %%% write results %%%
        %%%%%%%%%%%%%%%%%%%%%
                
        % write new camera data
        cpln=[cmax+(1:size(Con_fin_2,2)) % zeros(1,size(Con,2)) 
                c*ones(1,size(Con_fin_2,2)) 
                n*ones(1,size(Con_fin_2,2)) 
                Con_fin_2
                Vel_fin_2
                Int_fin_2];

        
        % max index
        cmax=cmax+size(Con,2);
        
        % write
        Cpln=cat(2,Cpln,cpln);
        
        %%%%%%%%%%%%%%%%%%%
        %%% end message %%%
        %%%%%%%%%%%%%%%%%%%
        
        % display message
        disp(['Identified ',num2str(size(cpln,2)),' bounding ellipse(s)'])
        
    end % n
    
end % c

%% Remove data

% trash processed frames previous data
Cpln=Cpln(:,~ismember(Cpln(3,:),frm_rem));
    
%% camera stats
if strcmp(plotting,'on')
    
    for c=1:prop.res(4)
        % make figure
        figure(c)
        
        % select camera
        C=Cpln(2,:)==c;
        
        % object size distribution
        subplot(1,2,1)
        [~,ax,~]=cvec2eshp(Cpln(4:9,C));
        L=sqrt(abs(ax(1,:).*ax(2,:)));
        histogram(log(L),linspace(log(min(L)),log(max(L)),50))
        xlabel('Length scale (log) $[\rm{m}]$','interpreter','latex')
        ylabel('No object occ. $[\rm{\#}]$','interpreter','latex')
        title('Object size','interpreter','latex')
        
        % object velocity distribution
        subplot(1,2,2)
        V=sqrt(sum(Cpln(10:11,C).^2,1));
        histogram(log(V),linspace(log(min(V)),log(max(V)),50))
        xlabel('Velocity scale (log) $[\rm{\#}]$','interpreter','latex')
        ylabel('No object occ. $[\rm{\#}]$','interpreter','latex')
        title('Object velocity','interpreter','latex')
        
        drawnow
        
    end
    
end

end

% %% play small movie of ellipse
% if strcmp(plotting,'on')
%     
%     % loop camera
%     for c=1:prop.res(4)
%         % select camera
%         C=Cpln(2,:)==c;
%         
%         % make figure
%         figure
%         
%         % loop frames
%         for n=frm_new
%             % get frame(s)
%             N=Cpln(3,:)>=n-ctrl.tres & Cpln(3,:)<=n+ctrl.tres ;
%             
%             % ellipses
%             E=reshape(cvec2pnts(-Cpln(4:9,N&C)),2,[]);
%             
%             % midpoint
%             X=cvec2eshp(Cpln(4:9,N&C));
%             
%             % velocity vector
%             U=Cpln(10:11,N&C);
%             
%             % plotting
%             hold on
%             plot(E(1,:),E(2,:),'.g','Markersize',0.5)
%             
%             % position
%             plot(X(1,:),X(2,:),'.b','MarkerSize',1)
%             
%             % velocity
%             quiver(X(1,:),X(2,:),U(1,:),U(2,:),0,'r','LineWidth',1)
%             
%             hold off
%             axis equal
%             xlabel('x-ref [m]','interpreter','latex')
%             ylabel('y-ref [m]','interpreter','latex')
%             title(['tracking cam ',num2str(c),' frm ',num2str(n)],'interpreter','latex')
%             
%             drawnow
%         end
%         
%     end
%     
% end