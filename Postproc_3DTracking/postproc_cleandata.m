function postproc_cleandata
%postproc_cleandata Process data cleaning on different segmentation criteria
%   
%   Here point are removed from the tracking, invalid velocity can be
%   reinterpolated in processing the trajectories. improving the traj qual.
%   

global folder date rec cal prop ctrl post traj plotting

%% load calibration 
Rmat=importdata([folder date cal vsl 'Rmat.mat']);
tvec=importdata([folder date cal vsl 'tvec.mat']);
Pmat=cell(size(Rmat));
for c=1:size(Pmat,2)
    Pmat{c}=krtm2pmat(eye(3),Rmat{c},tvec{c});
end
%% Create memory maps
% ( track index | frame number | object quadric | object velocity | peak intensity | res )
Pobj=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat']);

% remove small tracks entirely, source of most of the issue
% count occurence
[~,ia,~]=unique(Pobj([1 2],:)','rows');
l=unique(Pobj(1,ia));
L=histcounts(Pobj(1,ia),[l-1/2,max(l)+1/2]);
% segment
l=l(L>=5 );
L=ismember(Pobj(1,:),l);





% ( track index | indentification index | camera index | frame number | camera conic | camera velocity | object position | object velocity | res )
Plink = fload([folder date rec vsl 'Preproc3DTracking' vsl 'Plink.dat']);

% ( identification index | camera index | frame number | conic vector )
Cpln = fload([folder date rec vsl 'Preproc3DTracking' vsl 'Cpln.dat']);


for i=post.tproc(1):post.tproc(2)

    trackIndRem = zeros(1,0);
    for c=1:5
        disp(['Gathering duplicates from frame ' num2str(i) ' and view ' num2str(c)])
        % Find all Pobj corresponding to i
        Cplni = Cpln(:,ismember(Cpln(3,:),i) & ismember(Cpln(2,:),c));

        dumI=ismember(Pobj(2,:),i);
        dumJ=ismember(Pobj(1,:),Plink(1,Plink(3,:)==c));
        [~,ia,~]=unique(Pobj([1 2],:)','rows');
        l=unique(Pobj(1,ia));
        L=histcounts(Pobj(1,ia),[l-1/2,max(l)+1/2]);
        % segment
        l=l(L>=5 );
        L=ismember(Pobj(1,:),l);
        Pobji = Pobj(:,dumI&dumJ&L);

%             size(Pobji)

        % Project 
        quad=Pobji(3:12,:);

        % project
        con_pobj=cell2mat(cellfun(@(x)cmat2cvec(inv(Pmat{c}*(qvec2qmat(x)\Pmat{c}'))),...
            num2cell(quad,1),'UniformOutput',false));

        con_cpln = Cplni(4:9,:);
        [con_cpln_xcen,~,~] = cvec2eshp(con_cpln);
        [con_pobj_xcen,~,~] = cvec2eshp(con_pobj);

%             if 0
% 
%                figure()
%                hold all
%                scatter(con_cpln_xcen(1,:),con_cpln_xcen(2,:),'o')
%                scatter(con_pobj_xcen(1,:),con_pobj_xcen(2,:),'s')
%                axis equal
%                axis tight
% 
%             end

        % Find k nearest neighbors of Pobj

        [Idx,D] = knnsearch(con_cpln_xcen',con_pobj_xcen','K',1);
        % Find indices that are common
        [dumK,~,ic] = unique(Idx');

%             tic
        l=unique(Pobj(1,:));
        cnt=histcounts(Pobj(1,:),[l-1/2,max(l)+1/2]);

%             % track lengths of detected duplicates
%             tracklength = cnt(ismember(l,trackInd));
%             toc


        % loop over unique Indices
        for k=1:length(dumK) % for each unique detection in the image
%                k

%                tic
           trackInd = Pobji(1,ismember(Idx,dumK(k))); % find track indices that project to a unique detection in the image
           tracklength = cnt(ismember(l,trackInd));
%                toc



%                tic
           if length(trackInd)>1
               trackIndRem = cat(2,trackIndRem,trackInd(~ismember(tracklength,max(tracklength))));                   
           end
%                toc

% YADA YADA approaches for deletion, but ain't got the time/compute to push through this one by one

           % if length of track index is greater than 1, then
           % Pobji(:,ismember(Idx,dumK(k))) are Pobj columns that
           % correspond to the same ellipse detection in the image
%                tic
%                if length(trackInd)>1
% %                    trackInd
%                    % find track lengths                   
%                    l=unique(Pobj(1,:));
%                    cnt=histcounts(Pobj(1,:),[l-1/2,max(l)+1/2]);
%                    
%                    % track lengths of detected duplicates
%                    tracklength = cnt(ismember(l,trackInd));
%                    
%                    % indiscriminate track removal, keep only long tracks,
%                    % even if there are 2 of them. Faster, but
%                    % indiscriminate.
%                    trackIndRem = trackInd(~ismember(tracklength,max(tracklength)));
%                    Pobj(:,ismember(Pobj(1,:),trackIndRem)) = [];
%                    
%                    % find tracks that are longest
% %                    if length(trackInd(tracklength==max(tracklength))) == 2 % if there are two tracks equal to the maximum track length of ambiguous matches
% %                        trackIndreproj = zeros(size(trackInd));
% %                        for iii=1:length(trackInd)
% %                            trackIndreproj(iii) = mean(Pobj(17,ismember(Pobj(1,:),trackInd(iii))));
% %                        end
% %                        
% %                        trackIndRem = trackInd(~(ismember(tracklength,max(tracklength)) & trackIndreproj==min(trackIndreproj(ismember(tracklength,max(tracklength))))));
% %                        Pobj(:,ismember(Pobj(1,:),trackIndRem)) = [];
% %                        
% %                    else % delete remaining track indices from Pobj
% %                        
% %                        trackIndRem = trackInd(~ismember(tracklength,max(tracklength)));
% %                        Pobj(:,ismember(Pobj(1,:),trackIndRem)) = [];
% %                        
% %                    end
%                        
% %                    if  0%plot the selected track and the duplicates
% %                        [x_test,~,~]=qvec2eshp(Pobj(3:12,ismember(Pobj(1,:),trackInd)));
% %                        [x_sel,~,~]=qvec2eshp(Pobj(3:12,ismember(Pobj(1,:),trackInd(~ismember(trackInd,trackIndRem)))));
% %                        figure(1)
% %                        hold all
% %                        scatter3(x_test(1,:),x_test(2,:),x_test(3,:),'ro')
% %                        scatter3(x_sel(1,:),x_sel(2,:),x_sel(3,:),'b+')
% %                        axis equal
% %                        axis tight
% %                        drawnow
% %     %                    pause
% %     %                    clf
% %                    end
%                end
%                toc

        end
%             pause





    end
    % remove tracks
    disp(['Removing ' num2str(length(unique(trackIndRem))) ' short-tracked duplicates from frame ' num2str(i) ' and view ' num2str(c)])
    Pobj(:,ismember(Pobj(1,:),unique(trackIndRem))) = [];
            

end






%% Imposed overall segmentation

%-- #) remove nans from preprocessing
N=~isnan(sum(Pobj,1)) & ~isinf(sum(Pobj,1));

%-- #) remove short tracks
l=unique(Pobj(1,:));
cnt=histcounts(Pobj(1,:),[l-1/2,max(l)+1/2]);
%figure; histogram(cnt)
L=ismember(Pobj(1,:),l(cnt>=post.tlen)); % at least a second/ characteristic timescale
% 
% %-- #) segment domain
% X=nan*zeros(3,size(Pobj,2));
% ax=nan*zeros(3,size(Pobj,2));
% [X(:, N&L ),ax(:,N&L)]=qvec2eshp(Pobj(3:12, N&L ));
% X=~( sum( X<post.dom(:,1) | X>post.dom(:,2) ,1) >=1 ) ;
% 
% %-- #) segment time
% T=Pobj(2,:)>=post.tproc(1) & Pobj(2,:)<=post.tproc(2);
% 
% %-- #) Segment minimal camera cover triangulation/tracks

%to be added, not big fan of that, because thats a principle solution of the tracking
% 
% %-- #) Segment Reprojection error
% E=Pobj(end,:);
% %figure; histogram(E,50);
% E=E <= post.rpe;
% 
% %-- #) Segment Velocity
% V=sqrt(sum(Pobj(13:15,:).^2,1))*prop.fps;
% %figure; histogram(V,50);
% V=V <= post.vel;
% 
% %-- #) Segement Size
% S=(prod(ax,1)).^1/3;
% %figure; histogram(S,50);
% S=S <= post.siz;

%% Initiate overall segmentation

% valid data
val=N & L ;%& X & T & E & V & S ;

% write tracking and noise data
Tdata=Pobj( : ,  val ); % tracking
Ndata=Pobj( : , ~val ); % noise

%% Iterative segmentation trajectory kinematics

% initiate removal
G_val=true(1,size(Tdata,2));

% outlier filter
if ~strcmp('none',traj.outl)
    
    % loop frames
    for n=post.tproc(1):post.tproc(2)
        
        %-- #) message
        
        disp(['Filter trajectory dynamics frame ',num2str(n)])
        
        %-- #) get fit
        
        % get traj frame data
        N= ismember(Tdata(1,:),Tdata(1,Tdata(2,:)==n)) ...
            & Tdata(2,:)>=n-floor(traj.tfit/2) & Tdata(2,:)<=n+floor(traj.tfit/2) ;
        
        % get indexing
        ind=Tdata(1:2,N); % index and frame
        
        % get position 
        X = qvec2eshp( Tdata(3:12,N) );
        
        %figure; plot3(X(1,:),X(2,:),X(3,:),'.')
        
        % get velocity
        V=Tdata(13:15,N);
        
        % fit data
        Xdat=zeros(5,0);
        for t=-ctrl.tres:ctrl.tres
            Xdat=cat(2,Xdat,[ind+[0 t/(2*ctrl.tres+1)]'
                X+V*t/(2*ctrl.tres+1)]);
        end
        
        %figure; plot3(Xdat(3,:),Xdat(4,:),Xdat(5,:),'.')
        
        % compute trajectory
        coef=polytraj(Xdat, traj.tord); % front back stencil velocity
        
        %if there is velocity vector present
        
        %-- #) segment ghost dynamics
        
        % intiate removal
        g_val=true(1,nnz(N)); % indepencency bootstrapping G_val(N);
        
        % first and second derivative
        for k=1:2
            
            % derivative dynamics
            D=polyeval(coef,ind,k);
            
            % velocity object space
            d=sqrt(sum( D.^2 ,1));
            
            %figure; histogram(d,50); drawnow
            
            % iterate outlier removal
            [~,val]=stattest(d(g_val),traj.outl,traj.conf,traj.iter);
            
            % update removal
            g_val(g_val)=val;
            
            %figure; plot3(X(1,g_val),X(2,g_val),X(3,g_val),'.')
            %figure; histogram(d(g_val),50); drawnow
            
            % loop over dof views
            for c=1:prop.res(4)
                
                % derivative dynamics in camera reference
                d=[0 0 1]*Rmat{c}*D;
                
                %figure; histogram(d,50); drawnow
                
                % iterate outlier removal
                [~,val]=stattest(d(g_val),traj.outl,traj.conf,traj.iter);
                
                % update valid set
                g_val(g_val)=val;
                
                %figure; plot3(X(1,g_val),X(2,g_val),X(3,g_val),'.')
                %figure; histogram(d(g_val),50); drawnow
                
            end % c
            
        end % k
        
        %-- #) enforce g_val only on current frame
        
        g_val=g_val(ind(2,:)==n);
        
        % get traj frame data
        N=  Tdata(2,:)==n;
        
        %-- #) Update removal
        
        % then remove those points on the trajectory
        G_val(N)=G_val(N) & g_val;
        
        %-- #) message
        
        disp(['Identified ',num2str((1-nnz(g_val)/length(g_val))*100),' [%] outlier dynamics'])
        
    end % n
    
    %-- #) remove short tracks
    
    % indices
    l=unique(Tdata(1,G_val));
    
    % count index occurence, length frames
    cnt=histcounts(Tdata(1, G_val ),[l-1/2,max(l)+1/2]);
    
    %figure; histogram(cnt)
    
    % then remove whole trajectory
    G_val=G_val & ismember(Tdata(1,:),l(cnt>=post.tlen)); % remove also
    
end % if

%% Refine iterative segmentation

% valid data
val= G_val ; % 

% write tracking and noise data
Ndata=cat(2,Ndata,Tdata( : , ~val )); % noise
Tdata=Tdata( : , val ); % tracking

%% plotting
if strcmp(plotting,'on')
    
    figure
    
    for n=unique(Tdata(2,:))
        
        N=Tdata(2,:)>=n-post.tlen & Tdata(2,:)<=n+post.tlen;
        
        x = qvec2eshp(Tdata(3:12,N));
        col=sqrt(sum(Tdata(13:15,N).^2,1))*prop.fps;
        
        scatter3(x(1,:),x(2,:),x(3,:),[],col,'.')
        
        axis tight
        axis equal
        view(2)
        camproj('Perspective')
        xlim(post.dom(1,:))
        ylim(post.dom(2,:))
        zlim(post.dom(3,:))
        colorbar
        
        drawnow
        
    end
    
end

%% save
% write valid tracking data
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Tdata.dat'],Tdata,'w');

% write invalid noise data
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Ndata.dat'],Ndata,'w');

end
