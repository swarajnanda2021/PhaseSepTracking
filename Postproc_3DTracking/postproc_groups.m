function postproc_groups
%postproc_groups Process Groups in tracking data
%
%   Todo, 
%       Implement advanced segmentation and outlier filtering
%       Implement rotation axis slicing
%       Implement generating multiple groups

global folder date rec prop post grps plotting

%% Get data

Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);

% Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
% Acc=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);
% Quad=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat']);
% Tangent=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Tangent.dat']);

% also load eulerian reference / or rebin in group frame

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'Groups'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'Groups'])
end

%% Initiate variables
GroupIndex=zeros(2,0); % Group Index / Frame index
SegmentationIndex=zeros(3,0); % group Index / Traj Index /Frame Index
CenterGravity=zeros(3,0); % Center of Gravity
PrincipleAxisOrientation=zeros(3,0); % Principle Axis
PrincipleAxisSize=zeros(3,0); % Principle Size
GroupQuadric=zeros(10,0); % Group Quadric
GroupDisplacementField=zeros(7*5,0); % Displacement Field Group
FitResiduals=zeros(2,0); % Track Initial Conditions and Residuals
GroupVelocity=zeros(3,0); % Group Velocity
GroupRotation=zeros(3,0); % Group Rotation
GroupStrainRate=zeros(6,0); % Group Strain
GroupDivergence=zeros(1,0); % Group Divergence
CenterRotation=zeros(3,0); % Center Rotation
RotationAxis=zeros(3,0); % Rotation Axis as unit vector

%% Find the instateneous center and center of rotation
disp('post proc')
tic

% loop over time frames in the data
for n=post.tproc(1):post.tproc(2)
    
    % select time frame
    N=Index(2,:)>=n-(grps.tres-1)/2 & Index(2,:)<=n+(grps.tres-1)/2 ...
        & ismember(Index(1,:),Index(1,Index(2,:)==n));
    
    % indexing
    ind=Index(:,N);
    
    % select position data
    pos=Position( :,N );
    
    % figure; plot3(pos(1,:),pos(2,:),pos(3,:),'.'); axis equal; view(2)
    
    %%% outlier filtering will affect this
    
    % segment depth of field (kmeans is unstable, giving completely different result each time)
    [~,seg]=stattest(pos,grps.outl,grps.conf,grps.iter);
    seg=~ismember(ind(1,:),ind(1,~seg));
    
    % group indexing
    l=1; % only one group now
    
    %-- #) outlier filtering and segmentation
    
    % plot segmentation
    if strcmp(plotting,'on')
        
        % figure
        figure(1)
        cla
        
        % plot
        scatter3(pos(1,:),pos(2,:),pos(3,:),[],seg,'.')
        
        % layout
        axis equal
        view(2)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('segmentation')
        
    end
    
    % Group index
    GroupIndex=cat(2,GroupIndex,[l;n]); % prepared for indexing multiple groups
    
    % segment school
    pos=pos(:,seg);
    ind=ind(:,seg);
    
    % index mapping
    [ui,~,ii]=unique(ind(1,:));
    G_dat=sparse(ii,1:length(ii),ones(size(ii)),length(ui),length(ii));
    
    % Write segmentation
    SegmentationIndex=cat(2,SegmentationIndex,[l*ones(size(ui)) ; ui ; n*ones(size(ui))] ); % Segmentation to Tracking, what belong to this group?
    
    %-- #) Center of Gravity
    
    %%% convexhull would perhaps make bounding shape better, or shape alpha
    
%     % find midpoint and principle axis to shape
%     [xmid,ax,ang,c]=fitellipsepnts(pos(1:2,:),'free'); % fit ellipse
%     
%     %xe=cvec2pnts(c,100);
%     %hold on; plot3(xe(1,:),xe(2,:),0*xe(1,:),'.')
    
    % find midpoint and principle axis to shape
    [xmid,ax,ang,q]=fitellipsoidpnts(pos,'free'); % fit ellipsoid, sensitive to density points
    
    % write center of gravity
    CenterGravity=cat(2,CenterGravity,xmid);
    
    % write principle axis orientation
    PrincipleAxisOrientation=cat(2,PrincipleAxisOrientation,ang);
    
    % write principle axis
    PrincipleAxisSize=cat(2,PrincipleAxisSize,ax);
    
    % write group quadric
    GroupQuadric=cat(2,GroupQuadric,q);
    
    % plot center of gravity
    if strcmp(plotting,'on')
        
        % figure
        figure(2)
        
        % ellipsoid points
        xe=qvec2pnts(q,1000);
        
        % plot
        scatter3(pos(1,:),pos(2,:),pos(3,:),[],ind(2,:),'.')
        hold on
        plot3(xe(1,:),xe(2,:),xe(3,:),'r.')
        hold off
        
        % layout
        axis equal
        view(2)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Bounding Ellipsoid')
        
    end
    
    %-- #) Fit displacement field
    
    % map data to center of gravity pointcloud
    tdat=[ind-[0
                n]
        (pos-xmid)]; % centered 
    
    % fit linear instateneous displacement field
    [coef,fdat]=polydisp(tdat,[1 1 1 2],4); % all crossterms, few will unnessary drop by dt
    
    % write displacement field
    GroupDisplacementField=cat(2,GroupDisplacementField,reshape(coef,[],1));
    
    % compute mean and std residual per track
    resavg=fdat(5,:)*G_dat'./sum(G_dat,2)'; % avg
    resstd=sqrt((fdat(5,:)-resavg*G_dat).^2*G_dat'./sum(G_dat,2)'); % std
    
    % write residuals
    FitResiduals=cat(2,FitResiduals,[resavg ; resstd]);
    
    % plot displacement field fit
    if strcmp(plotting,'on')
        
        % figure
        figure(3)
        subplot(1,2,1)
        cla
        
        % evaluate displacement to check fit
        dX=dispeval(coef,fdat(1:4,:),[0 0 0 0]); % on normalized coordinates
        
        % fit positions
        X=fdat(2:4,:)+dX;
        
        %plot
        plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.','Color',[0.75 0.75 0.75])
        hold on
        scatter3(X(1,:),X(2,:),X(3,:),[],fdat(5,:),'.')
        quiver3(X(1,:),X(2,:),X(3,:),tdat(3,:)-X(1,:),tdat(4,:)-X(2,:),tdat(5,:)-X(3,:),1,'.r')
        hold off
        
        %layout
        view(2)
        axis equal
        colorbar
        drawnow
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Displacement Field Fit')
        
        % grid
        X=linspace(min(fdat(2,:)),max(fdat(2,:)),10);
        Y=linspace(min(fdat(3,:)),max(fdat(3,:)),10);
        Z=linspace(min(fdat(4,:)),max(fdat(4,:)),10);
        [X,Y,Z]=ndgrid(X,Y,Z);
        
        % indexed grid data
        xdat=[X(:)'
            Y(:)'
            Z(:)'];
        
        % displacement field
        dX=dispeval(coef,[zeros(1,size(xdat,2));xdat],[0 0 0 0]); % on normalized coordinates
        
        % displaced positions
        X=xdat+dX;
        
        % displacement field
        U=prop.fps*dispeval(coef,[zeros(1,size(xdat,2));xdat],[0 0 0 1]); % on normalized coordinates
        
        % plot
        subplot(1,2,2)
        scatter3(tdat(3,:),tdat(4,:),tdat(5,:),[],tdat(2,:),'.')
        hold on
        quiver3(X(1,:),X(2,:),X(3,:),U(1,:),U(2,:),0*U(3,:),1,'Color',[0.9290 0.6940 0.1250])
        
        % layout
        axis tight
        axis equal
        view(2)
        drawnow
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Points in Displacement Field')
        
    end
    
    % velocity at center of gravity of the group
    velcen=dispeval( coef, zeros(4,1) ,[0 0 0 1] )*prop.fps;
    
    % derive velocity gradient from displacement field
    velgrad=[dispeval( coef, zeros(4,1) ,[1 0 0 1] ),...
        dispeval( coef, zeros(4,1) ,[0 1 0 1] ),...
        dispeval( coef, zeros(4,1) ,[0 0 1 1] )]*prop.fps;
    
    % decompose symetric and skewsymetric part resp.
    [strain,rotation] = symskewdec(velgrad) ;
    
    % strain rate tensor
    strainvec = symmat2voightvec(strain); % voight vector [to be checked]
    
    % vorticity
    vorticity = skewmat2rotvec(rotation); % voricty vector from cross product matrix notation
    
    % compute divergence
    div=trace(velgrad); % sum dudx dudy..
    
    % Write Group Velocity
    GroupVelocity=cat(2,GroupVelocity,velcen);
    
    % Write Group Rotation
    GroupRotation=cat(2,GroupRotation,vorticity);
    
    % Write Group Strain
    GroupStrainRate=cat(2,GroupStrainRate,strainvec);
    
    % Write Group Divergence
    GroupDivergence=cat(2,GroupDivergence,div);
    
    % plot decomposed velocity field in uniform, rotation, strain
    if strcmp(plotting,'on')
        
        % figure
        figure(4)
        cla
        
        % velocity field rotation
        U=prop.fps*(velcen+0*xdat); % on normalized coordinates
        
        %plot
        subplot(1,3,1)
        plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.','Color',[0.75 0.75 0.75])
        hold on
        quiver3(X(1,:),X(2,:),X(3,:),U(1,:),U(2,:),U(3,:),1,'Color',[0.9290 0.6940 0.1250])
        
        % layout
        axis tight
        axis equal
        view(2)
        drawnow
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Uniform Component')
        
        % velocity field rotation
        U=prop.fps*(rotation*xdat); % on normalized coordinates
        
        %plot
        subplot(1,3,2)
        plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.','Color',[0.75 0.75 0.75])
        hold on
        quiver3(X(1,:),X(2,:),X(3,:),U(1,:),U(2,:),U(3,:),1,'Color',[0.9290 0.6940 0.1250])
        
        % layout
        axis tight
        axis equal
        view(2)
        drawnow
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Rotation Component')
        
        % velocity field rotation
        U=prop.fps*(strain*xdat); % on normalized coordinates
        
        %plot
        subplot(1,3,3)
        plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.','Color',[0.75 0.75 0.75])
        hold on
        quiver3(X(1,:),X(2,:),X(3,:),U(1,:),U(2,:),U(3,:),1,'Color',[0.9290 0.6940 0.1250])
        
        % layout
        axis tight
        axis equal
        view(2)
        drawnow
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Strain Component')
        
    end
    
    % vorticity frame, topview rotation
    b3=vorticity/norm(vorticity,2); % new z-dir in vorticity direction
    b1=cross(b3,[0 0 1]')./norm(cross(b3,[0 0 1]'),2);
    b2=cross(b3,b1);
    R=[b1 b2 b3]'; % project on that basisvector
    
    % find center rotation
    A=R*rotation*R'; % transform rotation field
    vel0=R*dispeval(coef,zeros(4,1),[0 0 0 1]);
    xrot=-A(1:2,1:2)\vel0(1:2);
    
    % data rotated according vorticity vector
    if strcmp(plotting,'on')
        
        % velocity field rotation
        U=prop.fps*(velcen+rotation*xdat); % on normalized coordinates
        
        % figure
        figure(5)
        cla
        
        %plot
        subplot(1,2,1)
        plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.','Color',[0.75 0.75 0.75])
        hold on
        quiver3(X(1,:),X(2,:),X(3,:),U(1,:),U(2,:),U(3,:),1,'Color',[0.9290 0.6940 0.1250])
        quiver3(0,0,0,R(1,1),R(1,2),R(1,3),0,'r')
        quiver3(0,0,0,R(2,1),R(2,2),R(2,3),0,'g')
        quiver3(0,0,0,R(3,1),R(3,2),R(3,3),0,'b')
        hold off
        
        % layout
        axis tight
        axis equal
        view(2)
        drawnow
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Rotation Component')
        
        % figure
        subplot(1,2,2)
        
        tacc=R*tdat(3:5,:);
        
        % displacement field
        Uacc=R*U; % on normalized coordinates
        Xacc=R*X;
        plot3(tacc(1,:),tacc(2,:),tacc(3,:),'.','Color',[0.75 0.75 0.75])
        hold on
        quiver3(Xacc(1,:),Xacc(2,:),Xacc(3,:),Uacc(1,:),Uacc(2,:),Uacc(3,:),1,'.','Color',[0.9290 0.6940 0.1250])
        quiver3(0,0,0,1,0,0,0,'r')
        quiver3(0,0,0,0,1,0,0,'g')
        quiver3(0,0,0,0,0,1,0,'b')
        plot3(xrot(1),xrot(2),0,'.r','MarkerSize',100)
        hold off
        
        % layout
        axis tight
        axis equal
        view(2)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Trajectories in Vortex frame')
        drawnow
        
    end
    
    % center of rotation in 3D space
    xrot=xmid+R\[xrot;0];
    axrot=b3;
    
    % plot center of gravity
    if strcmp(plotting,'on')
        
        % figure
        figure(6)
        
        % plot
        plot3(pos(1,:),pos(2,:),pos(3,:),'.')
        hold on
        quiver3(xrot(1,:)-2*axrot(1),xrot(2,:)-2*axrot(2),xrot(3,:)-2*axrot(3)...
            ,40*axrot(1),40*axrot(2),40*axrot(3),1,'r.','LineWidth',5)
        plot3(xrot(1,:),xrot(2,:),xrot(3,:),'g.','Markersize',25)
        hold off
        
        % layout
        axis equal
        view(2)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Bounding Ellipsoid')
        
    end
    
    % Center Rotation
    CenterRotation=cat(2,CenterRotation,xrot);
    
    % Unit Vector Axis
    RotationAxis=cat(2,RotationAxis,axrot);
    
    %-- #) message
    disp(['Processed group(s) at frame ',num2str(n)])
    
end


%% save

% Group Segmentation
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupIndex.dat'],GroupIndex,'w');

% Group Segmentation
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'SegmentationIndex.dat'],SegmentationIndex,'w');

% Center of Gravity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterGravity.dat'],CenterGravity,'w');

% Principle Axis
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'PrincipleAxisOrientation.dat'],PrincipleAxisOrientation,'w');

% Principle Size
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'PrincipleAxisSize.dat'],PrincipleAxisSize,'w');

% Group Quadric
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupQuadric.dat'],GroupQuadric,'w');

% Displacement Field Group
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupDisplacementField.dat'],GroupDisplacementField,'w');

% Fit Residuals
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'FitResiduals.dat'],FitResiduals,'w');

% Group Velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupVelocity.dat'],GroupVelocity,'w');

% Group Velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupRotation.dat'],GroupRotation,'w');

% Group Velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupStrainRate.dat'],GroupStrainRate,'w');

% Center Rotation
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterRotation.dat'],CenterRotation,'w');

% Group Divergence
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupDivergence.dat'],GroupDivergence,'w');

% Rotation Axis
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'RotationAxis.dat'],RotationAxis,'w');

end


%     % define indexing groups by branching adjacency in rank
%     
%     % loop
%     for i=1:Ngroup
%         
%         % create path for processing trajectory reference
%         if exist([folder date rec !'\postproc\groupframe\group_',num2str(i)],'dir')~=7
%                 mkdir([folder date rec !'\postproc\groupframe\group_',num2str(l)])
%                 end
%                 
%                 % save there
%                 Tangent=cell2mat(Tvec);
%                 fsave([folder,date,rec,...
%                 !'\postproc\\groupframe\group_',num2str(i),'\Tangent.dat']...
%                 ,Tangent,'w');
%             
%         end
%         
%     end
%     
%     toc
    