function postproc_trajectories
%postproc_trajectories Process Trajectories by smoothing and differentiation
%   filter. Root to other trajectory processings

global folder date rec prop ctrl post traj

%% Create memory maps
Tdata=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Tdata.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories'])
end

%% Initiate files to write
Index=zeros(2,0);
Time=zeros(1,0);
Quad=zeros(10,0);
Curve=zeros(3*(traj.tord+1),0);
Pos=zeros(3,0);
Vel=zeros(3,0);
Acc=zeros(3,0);
Axi=zeros(3,0);
Ang=zeros(3,0);
Accuracy = zeros(1,0);

FitRes = zeros(1,0);

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for n=post.tproc(1):post.tproc(2)%unique(Tdata(2,:))
    
	% define traj
    N=  ( ismember(Tdata(1,:),Tdata(1, Tdata(2,:)>=n-floor(traj.tfit/2) & Tdata(2,:)<=n ) ) ...
        & ismember(Tdata(1,:),Tdata(1, Tdata(2,:)>=n & Tdata(2,:)<=n+floor(traj.tfit/2) ) ) ) ... % allows interpolation
		& Tdata(2,:)>=n-floor(traj.tfit/2) & Tdata(2,:)<=n+floor(traj.tfit/2) ;
    
    
                
     
    % get indexing
    ind=Tdata(1:2,N); % index and frame
    
    % get quadrics
    Q=Tdata(3:12,N);
    
    %Xq=reshape(qvec2pnts(Q,10),3,[]);
    %figure; plot3(Xq(1,:),Xq(2,:),Xq(3,:),'.'); axis tight; axis equal
    
    % get position shape and orientation
    [ X, ax, ang ] = qvec2eshp( Q );
    
    %figure; scatter3(X(1,:),X(2,:),X(3,:),[],ind(2,:),'.'); axis tight; axis equal
    %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis tight; axis equal
    
    % get velocity
    V=Tdata(13:15,N);
    
    %figure; quiver3(X(1,:),X(2,:),X(3,:),V(1,:),V(2,:),V(3,:),0); axis tight; axis equal
    
	%-- #) Pad boundary start and end knots which are gradient preservind (Neumann condition)
	
	% current index mapping
	[ul,~,li]=unique(ind(1,:));
	G_sol=sparse(li,1:length(li),ones(size(li)),length(ul),length(li));
	
%     figure; spy(G_sol)
% 	figure; plot(sum(G_sol,1),'.'); % should be ones.
	
	% boundary indices
	[~,ni_min]=max( G_sol*spdiags( max(ind(2,:))+1-ind(2,:)' ,0,size(ind,2),size(ind,2)) ,[],2); % min index begin track
	[~,ni_max]=max( G_sol*spdiags(                 ind(2,:)' ,0,size(ind,2),size(ind,2)) ,[],2); % max index end track
	
    %figure; plot(ni_min,'.'); % should be ones.
	%hold on; plot(ni_max,'.'); % should be ones.
	
	% append boundary indexing
	N_min= ind(2,ni_min)*G_sol; % begin index
	N_max= ind(2,ni_max)*G_sol; % end index
	ind=cat(2,[ ind(1,:) ; 2*N_min-ind(2,:)], ind ,[ ind(1,:) ; 2*N_max-ind(2,:) ] );
	
	% append boundary points
	X_min= X(:,ni_min)*G_sol; % begin position
	X_max= X(:,ni_max)*G_sol; % end position
	X=cat(2, 2*X_min-X , X , 2*X_max-X );
    
   
    
%     figure; scatter3(X(1,:),X(2,:),X(3,:),[],ind(2,:),'.'); axis tight; axis equal; view(2);drawnow;
    %sel=round(median(ind(1,:)));
    %figure; plot(ind(2,ind(1,:)==sel),X(1,ind(1,:)==sel),'.'); 
    
	% append boundary velocity
	V_min= V(:,ni_min)*G_sol; % begin velocity
	V_max= V(:,ni_max)*G_sol; % end velocity
	V=cat(2,2*V_min-V,V,2*V_max-V);
	
	% append boundary axis
	ax_min= ax(:,ni_min)*G_sol; % begin axis
	ax_max= ax(:,ni_max)*G_sol; % end axis
	ax=cat(2,2*ax_min-ax,ax,2*ax_max-ax);
	
	% append boundary angles
	ang_min= ang(:,ni_min)*G_sol;% begin angles
	ang_max= ang(:,ni_max)*G_sol;% end axis
	ang=cat(2,2*ang_min-ang,ang,2*ang_max-ang);
	
	% keep valid stencil points
	g_ind=ind(2,:)>=n-floor(traj.tfit/2) & ind(2,:)<=n+floor(traj.tfit/2);
	ind=ind(:,g_ind);
    X=X(:,g_ind);
	V=V(:,g_ind);
	ax=ax(:,g_ind);
	ang=ang(:,g_ind);
    
	%hold on; plot(ind(2,ind(1,:)==sel),X(1,ind(1,:)==sel),'o'); 
    
	% remove double points at boundary
	[~,ai,~]=unique(ind','rows');
    
	ind=ind(:,ai);
	X=X(:,ai);
	V=V(:,ai);
	ax=ax(:,ai);
	ang=ang(:,ai);
	
    %hold on; plot(ind(2,ind(1,:)==sel),X(1,ind(1,:)==sel),'+'); 
    %figure; scatter3(X(1,:),X(2,:),X(3,:),[],ind(2,:),'.'); axis tight; axis equal
    
	%-- #) compute smoothing and differentiation
    
    % fit data
    Xdat=zeros(5,0);
    for t=-ctrl.tres:ctrl.tres
        Xdat=cat(2,Xdat,[ind+[0 t/(2*ctrl.tres+1)]'
            X+V*t/(2*ctrl.tres+1)]);
    end

	% compute trajectory
	[coef,residual]=polytraj(Xdat, traj.tord); % front back stencil velocity
    
    % compute angles
    coefang=polytraj([ind ; ang],traj.tord);
	
    % compute axis
    coefax=polytraj( [ind ; ax] ,traj.tord);
	
	%-- #) evaluate filtered data
	
	% remove stencil from indexing
    [~,ai,~]=unique(ind(1,:));
    
    res = residual(1,ai);
    
	ind=[ind(1,ai)
        n*ones(1,length(ai))]; % allows to interpolate
	
    % real time data
    t=ind(2,:)/prop.fps;
    
    % consistent position
    X=polyeval(coef,ind,0);
    
    %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis tight; axis equal;view(2)
    
	% velocity
    V=polyeval(coef,ind,1)*prop.fps;
    
%     figure; quiver3(X(1,:),X(2,:),X(3,:),V(1,:)/prop.fps,V(2,:)/prop.fps,V(3,:)/prop.fps,0); axis tight; axis equal
    
    % accelaration
    A=polyeval(coef,ind,2)*prop.fps*prop.fps;
    
    % filter axis
    ax=polyeval(coefax,ind,0);
    
    % filter angles
    ang=polyeval(coefang,ind,0);
    
    % evaluate smoothed conics
    Q=eshp2qvec(X,ax,ang); % consistency with position
    
%     Xq=reshape(qvec2pnts(Q,10),3,[]);
%     figure; plot3(Xq(1,:),Xq(2,:),Xq(3,:),'.'); axis tight; axis equal; view(2)
    
    % reshape coefficients
    coefx=reshape(coef(3,:),traj.tord+1,[]);
    coefy=reshape(coef(4,:),traj.tord+1,[]);
    coefz=reshape(coef(5,:),traj.tord+1,[]);
    
    coef=cat(1,coefx,coefy,coefz);
    
	%-- #) write output
	
    % write output : share headers
    Index=cat(2,Index,ind);
    Time=cat(2,Time,t);
    Quad=cat(2,Quad,Q);
    Curve=cat(2,Curve,coef);
    Pos=cat(2,Pos,X);
    Vel=cat(2,Vel,V); 
    Acc=cat(2,Acc,A);
    Axi=cat(2,Axi,ax);
    Ang=cat(2,Ang,ang);
    Accuracy = cat(2,Accuracy,Tdata(17,ai));
    
    FitRes = cat(2,FitRes,res);
    
    
	%-- #) Message
	
    % message
    disp(['post processed trajectories frame ',num2str(n)])
	
end

toc

%% save
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat'],Index,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat'],Time,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat'],Quad,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Curve.dat'],Curve,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat'],Pos,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat'],Vel,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat'],Acc,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'EllipsoidAxis.dat'],Axi,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'EllipsoidAngles.dat'],Ang,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accuracy.dat'],Accuracy,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'FitResidual.dat'],FitRes,'w');

end