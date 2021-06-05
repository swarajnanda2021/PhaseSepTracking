function postproc_frenetserret
%postproc_frenetserret Process frenet serret frame work
%
%   Input trajectory processing
%
%   [To Be Checked]

global folder date rec post traj

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Curve=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Curve.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret'])
end

%% Initiate files to write
Tvec=zeros(3,0);
Nvec=zeros(3,0);
Bvec=zeros(3,0);
Cur=zeros(1,0);
Tor=zeros(1,0);
Dvec=zeros(3,0);
Plncur=zeros(1,0);
Rad=zeros(1,0);
Pitch=zeros(1,0);
Vvec=zeros(3,0);
Avec=zeros(3,0);
Ome=zeros(1,0);

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for n=post.tproc(1):post.tproc(2)%unique(Index(2,:))
    % get trajectory of interest
    NN=Index(2,:)==n ; % rem nans 
	
	% indexing
	ind=Index(1:2,NN);
    
    % get coefficient that describe the curve along the trajectory
    coef=Curve(:,NN);
    coefi=reshape( repmat(ind(1,:),traj.tord+1,1          ) ,1,[]);
    coefo=reshape( repmat((1:4)'  ,1          ,size(ind,2)) ,1,[]);
    coefx=reshape(coef(                   1:traj.tord+1  ,:),1,[]);
    coefy=reshape(coef(  (traj.tord+1) + (1:traj.tord+1) ,:),1,[]);
    coefz=reshape(coef(2*(traj.tord+1) + (1:traj.tord+1) ,:),1,[]);
    coef=cat(1,coefi,coefo,coefx,coefy,coefz);
    
    % get physical velocity along the curve
    vel=Velocity(:,NN);
    
    % get physical accelaration along the curve
    acc=Accelaration(:,NN);
    
    % get coefficient x y z coordinates as they decouple for polynomials
    
    % define derivative to the path (similar but not same as vel and acc!)
    dR=polyeval(coef,ind,1);
    ddR=polyeval(coef,ind,2);
    dddR=polyeval(coef,ind,3);
    
    % derive tangent vector to curve of the frenet trihedron
    T=dR./sqrt( sum( dR.^2 , 1 ) );
    
    % derive the normal vector to curve of the frenet trihedron
    N=cross( dR , cross( ddR, dR ) ) ./ ...
        ( sqrt( sum( dR.^2 , 1 ) ) .* ...
        sqrt( sum(  cross( ddR, dR ).^2 , 1 ) ) ); 
    
    % derive the binormal vector to curve of the frenet trihedron
    B= cross( dR, ddR ) ./ sqrt( sum( cross( dR, ddR ).^2 , 1 ) );
    
    % derive curvature to curve
    cur=sqrt( sum( cross( dR, ddR ).^2 , 1) ) ./ ...
        sqrt( sum( dR.^2 , 1) ).^3 ;
    
    % derive torsion to curve
    tor= dot( dddR, cross( dR, ddR ) ) ./ ... % sqrt( sum( dot( dR, cross( ddR, dddR ) ).^2 , 1) ) ./ ...
        sqrt( sum( cross( dR, ddR ).^2 , 1) ).^2 ;
    
    % Darboux vector
    D=tor.*T+cur.*B;
    
    % plane curvature
    plncur=sqrt(cur.^2 + tor.^2);
    
    % radius
    rad=cur./(cur.^2+tor.^2);
    
    % pitch
    pitch=2*pi*tor./(cur.^2 + tor.^2);
    
    % (linear) velocity to the frenet frame
    V=[dot(T,vel)
        dot(N,vel)
        dot(B,vel)];
    
    % accelaration to the frenet frame
    A=[dot(T,acc)
        dot(N,acc)
        dot(B,acc)];
    
    % angular velocity    dot(T,vel).*D;%
    O=sign( tor ).*dot(T,vel).*sqrt( cur.^2 + tor.^2 ); 
    
    % write
    Tvec=cat(2,Tvec,T);
    Nvec=cat(2,Nvec,N);
    Bvec=cat(2,Bvec,B);
    Cur=cat(2,Cur,cur);
    Tor=cat(2,Tor,tor);
    Dvec=cat(2,Dvec,D);
    Plncur=cat(2,Plncur,plncur);
    Rad=cat(2,Rad,rad);
    Pitch=cat(2,Pitch,pitch);
    Vvec=cat(2,Vvec,V);
    Avec=cat(2,Avec,A);
    Ome=cat(2,Ome,O);
    
    % message
    disp(['post processed frenet serret framework ',num2str(n)])
end

toc

%% save
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Tangent.dat'],Tvec,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Normal.dat'],Nvec,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Binormal.dat'],Bvec,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Curvature.dat'],Cur,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Torsion.dat'],Tor,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'DarbouxVector.dat'],Dvec,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'PlaneCurvature.dat'],Plncur,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'HelixRadius.dat'],Rad,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'HelixPitch.dat'],Pitch,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'HelixOmega.dat'],Ome,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Velocity.dat'],Vvec,'w');

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'Frenetserret' vsl 'Accelaration.dat'],Avec,'w');

% add here
% angular velocity
% angular accelaration

end