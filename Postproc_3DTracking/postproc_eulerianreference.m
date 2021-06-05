function postproc_eulerianreference
%postproc_eulerianreference Process eulerian reference
%   Binning the data
%
%   Commute coherent group motion from fitting a displacement field
%
%   John: Integrate the outlier filtering in this script

%% get globals
global folder date rec post prop eulr

%% Get data
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);
load([folder date rec vsl 'Postproc3DTracking' vsl 'TrackQualityMetric.mat'],'t_spurious');

% Curvature=fload([folder date rec vsl 'Postproc' vsl 'Frenetserret' vsl 'Curvature.dat']);

%% make folder for data
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference'])
end

%% define indexed grid

% generate binning grid
[Y,X,Z]=meshgrid(post.dom(2,1):eulr.gres:post.dom(2,2),...
    post.dom(1,1):eulr.gres:post.dom(1,2),...
    post.dom(3,1):eulr.gres:post.dom(3,2));

% write grid as unstructured data
vecgrid.IND=1:numel(X);
vecgrid.X=[X(:)'
    Y(:)'
    Z(:)'];

% clear variables
clear X Y Z

%% initiate data sets

% griddata
Gdat_binindex=zeros(2,0); % index
Gdat_time=zeros(1,0); % time

% structure functions
Gdat_density=zeros(1,0); % number density
Gdat_position=zeros(3,0); % point position correction
Gdat_binnedvel=zeros(6,0); % binned velocity vector
Gdat_polarization=zeros(6,0); % polarization vector
Gdat_linearvel=zeros(2,0); % linear velocity vector magnitude
Gdat_binnedacc=zeros(6,0); % binned velocity vector

% motion function
Gdat_dispfield=zeros(7*5,0); % displacement field
% Gdat_fdat=zeros(7,0); % fitted trajectories on grid point memory
% intrusive
Gdat_fitresidual=zeros(2,0); % residual mean and std / summarizing any type of deviation

% kinematics
Gdat_velocity=zeros(3,0); % velocity vector
Gdat_accelaration=zeros(3,0); % accelaration vector
Gdat_velgradient=zeros(9,0); % velocity gradient
Gdat_strainrate=zeros(6,0); % symmetric strainrate voight vector
% principle strainrate [eigen vector and values]
Gdat_vorticity=zeros(3,0); % skew symmetric vorticity vector
Gdat_divergence=zeros(1,0); % divergence
Gdat_matderivative=zeros(3,0); % material derivative
% massflux

% % trajectory statistics
% Gdat_cur=zeros(2,0); % curvature field

%% loop frames in data, preserve time resolution
for n=post.tproc(1):post.tproc(2)%unique(Index(2,:))
    
    % select time frame
    N=Index(2,:)>=n-(eulr.tres-1)/2 & Index(2,:)<=n+(eulr.tres-1)/2 ...
        & ismember(Index(1,:),Index(1,Index(2,:)==n));
    
    
    % Long tracks
    [~,test,~] = unique(Index(1,:),'rows');
    l = unique(Index(1,:));
    L=histcounts(Index(1,:),[l-1/2,max(l)+1/2]);
    l=l(L>=eulr.lmin);
    L=ismember(Index(1,:),l);
    
    % select based on Lmin
    
    outl = ~ismember(Index(1,:),unique(t_spurious));
    % indexing
    ind=Index(:,N&L&outl);
    
    % select data
    pos=Position( :,N&L&outl);
    
    % trajectory velocity
    vel=Velocity( :,N&L&outl);
    
    % linear velocity
    linvel=sqrt(sum(vel.^2,1));
    
    % polarization
    pol=vel./linvel; % see how the above connects
    
    % trajectory acc
    acc=Accelaration( :,N&L&outl);
    
    %     cur=Curvature( :,N );
    
    % loop grid nodes
    for i=vecgrid.IND
        
        % data selection
        S=ismember(ind(1,:),ind(1,...
            sqrt(sum((pos-[vecgrid.X(1,i)
                            vecgrid.X(2,i)
                            vecgrid.X(3,i)]).^2,1))<=eulr.sres ...
            & ind(2,:)==n )) ;
        
        % local vector validation
        if ~strcmp('none',eulr.outl)
            
            %figure; quiver3(zeros(1,nnz(S)),zeros(1,nnz(S)),zeros(1,nnz(S)),vel(1,S),vel(2,S),vel(3,S),0); axis equal
            
            % track indiced
            tind=ind(1,S);
            
            % select test
            [~,val]=stattest(vel(:,S),eulr.outl,eulr.conf,eulr.iter);
            
            % make sure the full track is select
            val=ismember(tind,tind(val));
            
            % update segmentation
            S(S)=val;
            
            %hold on;quiver3(zeros(1,nnz(S)),zeros(1,nnz(S)),zeros(1,nnz(S)),vel(1,S),vel(2,S),vel(3,S),0);axis equal; view(2)
            
        end
        
        % process data
        if nnz(S)>0 % anything there?
            
            % write bin index
            Gdat_binindex=cat(2,Gdat_binindex,[i
                n]);
            
            % write time
            Gdat_time=cat(2,Gdat_time,n/prop.fps);
            
            % density per unit volume
            dens=nnz(S)/eulr.tres/(4/3*pi*eulr.sres^3); % average number objects
            
            % write object density within timespan
            Gdat_density=cat(2,Gdat_density,dens);
            
            % compute grid position
            gridpos=mean(pos(:,S),2)-vecgrid.X(:,i);
            
            % write averaged position
            Gdat_position=cat(2,Gdat_position,gridpos);
            
            % write binned velocity vector
            Gdat_binnedvel=cat(2,Gdat_binnedvel,[mean(vel(:,S),2)
                                                    std(vel(:,S),[],2)]);
            
            % binned velocity vector
            Gdat_linearvel=cat(2,Gdat_linearvel,[mean(linvel(:,S),2)
                                                    std(linvel(:,S),[],2)]);
            
            % binned polarization vector
            Gdat_polarization=cat(2,Gdat_polarization,[mean(pol(:,S),2)
                                                    std(pol(:,S),[],2)]);
            
            % write binned velocity vector
            Gdat_binnedacc=cat(2,Gdat_binnedacc,[mean(acc(:,S),2)
                                                    std(acc(:,S),[],2)]);
            
            % spatial derivative min number fish = 1+3+1 (3x first order, always at least 1 overfit)
            % temporal derivative min number frame = 1+2+1 (1x second order, always at least 1 overfit)
            % 50 [%] spatio temporal data completeness / not neccessarily needed
            if length(unique(ind(1,S)))>=5 && length(unique(ind(2,S)))>=4 % && ...
                   % size(unique(ind(:,S)','rows'),1)>=0.5*length(unique(ind(1,S)))*length(unique(ind(2,S)))
                
                % spatio temporally centered displacement data trajectory
                tdat=[ind(:,S)-[0
                    n]
                    pos(:,S)-[vecgrid.X(1,i)
                    vecgrid.X(2,i)
                    vecgrid.X(3,i)]];
                
                % figure; scatter3(tdat(3,:),tdat(4,:),tdat(5,:),[],tdat(2,:),'.'); xlim([-eulr.sres eulr.sres]) ; ylim([-eulr.sres eulr.sres]) ; zlim([-eulr.sres eulr.sres]); 
                
                % fit displacement field
                [coef,fdat]=polydisp(tdat,[1 1 1 2],4);
                
%                 % vector validation (based only on local time)
%                 if strcmp('median',eulr.outl)
%                     
%                     % filter median
%                     [~,val] = testmedian(fdat(5,:),2);
%                     
%                     % update segmentation
%                     S(S)=val;
%                     
%                 end
                
                % compute residual statistics
                res=[mean(fdat(5,:))
                    std(fdat(5,:))];
                
                % write displacement field
                Gdat_dispfield=cat(2,Gdat_dispfield,reshape(coef,[],1)); % displacement field
                
%                 % write trajectory fit data
%                 Gdat_fdat=cat(2,Gdat_fdat,[i*ones(1,size(fdat,2))
%                     n*ones(1,size(fdat,2))
%                     fdat]); % mem error
                
                % write fit residual
                Gdat_fitresidual=cat(2,Gdat_fitresidual,res); % metric units
                
                % evaluate displacement fit trajectory
                % dX= dispeval( coef, fdat(1:4,:) ,[0 0 0 0] );
                % compute fit positions
                % X=fdat(2:4,:)+dX;
                % dumerr=fdat(5,:);
                % dumerr( fdat(5,:)>(mean(fdat(5,:))+3*std(fdat(5,:))) )=mean(fdat(5,:))+3*std(fdat(5,:));
                % figure; plot3(tdat(3,:),tdat(4,:),tdat(5,:),'.','Color',[0.75 0.75 0.75],'MarkerSize',10)
                % hold on; scatter3(X(1,:),X(2,:),X(3,:),200,dumerr,'.'); colorbar;
                % quiver3(X(1,:),X(2,:),X(3,:),tdat(3,:)-X(1,:),tdat(4,:)-X(2,:),tdat(5,:)-X(3,:),0,'.r'); axis equal; xlim([-eulr.sres eulr.sres]) ; ylim([-eulr.sres eulr.sres]) ; zlim([-eulr.sres eulr.sres])
                
                % derive velocity from displacement field
                U=dispeval( coef, [0;gridpos] ,[0 0 0 1] )*prop.fps;
                
                % quiver3(gridpos(1),gridpos(2),gridpos(3),U(1)*prop.fps,U(2)*prop.fps,U(3)*prop.fps,0,'r'); xlim([-eulr.sres eulr.sres]) ; ylim([-eulr.sres eulr.sres]) ; zlim([-eulr.sres eulr.sres])
                
%                 % derive velocity from displacement field
%                 U=dispeval( coef, [0;0;0;0] ,[0 0 0 1] );
%                 
%                 % quiver3(0,0,0,U(1)*prop.fps,U(2)*prop.fps,U(3)*prop.fps,0,'r'); xlim([-eulr.sres eulr.sres]) ; ylim([-eulr.sres eulr.sres]) ; zlim([-eulr.sres eulr.sres])
                
                %%% if U spurious against mean vel + std vel
                
                % write velocity
                Gdat_velocity=cat(2,Gdat_velocity,U);
                
                % derive accelaration from displacement field
                dUdt=dispeval( coef, [0;gridpos] ,[0 0 0 2] )*prop.fps*prop.fps;
                
                % quiver3(0,0,0,dUdT(1),dUdt(2),dUdt(3),0); xlim([-eulr.sres eulr.sres]) ; ylim([-eulr.sres eulr.sres]) ; zlim([-eulr.sres eulr.sres])
                
                % write velocity
                Gdat_accelaration=cat(2,Gdat_accelaration,dUdt);
                
                % derive velocity gradient from displacement field
                velgrad=[dispeval( coef, [0;gridpos] ,[1 0 0 1] ),...
                    dispeval( coef, [0;gridpos] ,[0 1 0 1] ),...
                    dispeval( coef, [0;gridpos] ,[0 0 1 1] )]*prop.fps;
                
                % write velocity gradient
                Gdat_velgradient=cat(2,Gdat_velgradient,reshape(velgrad,[],1));
                
                % decompose symetric and skewsymetric part resp.
                [strain,rotation] = symskewdec(velgrad) ;
                
                % strain rate tensor
                strainvec = symmat2voightvec(strain); % voight vector 
                
                % write strain rate
                Gdat_strainrate=cat(2, Gdat_strainrate,strainvec) ; % symmetric strainrate
                
                % vorticity
                vor = skewmat2rotvec(rotation); % voricty vector from cross product matrix notation
                
                % quiver3(0,0,0,vor(1),vor(2),vor(3),0); xlim([-eulr.sres eulr.sres]) ; ylim([-eulr.sres eulr.sres]) ; zlim([-eulr.sres eulr.sres])
                
                % write vorticity vector
                Gdat_vorticity=cat(2,Gdat_vorticity,vor); % skew symmetric vorticity
                
                % compute divergence
                div=trace(velgrad); % sum dudx dudy..
                
                % write divergence
                Gdat_divergence=cat(2,Gdat_divergence,div); % divergence
                
                % compute material derivative
                matder=dUdt+velgrad'*U; % [to be checked]
                
                % write material derivative
                Gdat_matderivative=cat(2,Gdat_matderivative,matder); % material derivative
                
                % Massflux
                % Principle strains
                
                %             % compute average velocity
                %             Gdat_vel=cat(2,Gdat_vel,...
                %                 [mean(vel(:,S),2)
                %                 std(vel(:,S),[],2)]);
                %
                %             % compute average accelaration
                %             Gdat_acc=cat(2,Gdat_acc,...
                %                 [mean(acc(:,S),2)
                %                 std(acc(:,S),[],2)]);
                %
                %             % compute average director polarization
                %             Gdat_dirpol=cat(2,Gdat_dirpol,...
                %                 [mean(dir(:,S),2)
                %                 std(dir(:,S),[],2)]);
                %
                %             % compute average curvature
                %             Gdat_cur=cat(2,Gdat_cur,...
                %                 [mean(cur(:,S),2)
                %                 std(cur(:,S),[],2)]);
                
            else % write nan's s.t. we can see unevaluated cells
                
%                 % write object density within timespan
%                 Gdat_density=cat(2,Gdat_density,nan);
%                 
%                 % write object density within timespan
%                 Gdat_position=cat(2,Gdat_position,nan(3,1));
%                 
%                 % write binned velocity vector
%                 Gdat_binnedvel=cat(2,Gdat_binnedvel,nan(6,1));
%                 
%                 % write binned velocity vector
%                 Gdat_linearvel=cat(2,Gdat_linearvel,nan(2,1));
%                 
%                 % write binned polarization vector
%                 Gdat_polarization=cat(2,Gdat_polarization,nan(6,1));
%                 
%                 % write binned velocity vector
%                 Gdat_binnedacc=cat(2,Gdat_binnedacc,nan(6,1));
                
                % write displacement field
                Gdat_dispfield=cat(2,Gdat_dispfield,nan(5*7,1)); % displacement field
                
%                 % write trajectory fit data
%                 Gdat_fdat=cat(2,Gdat_fdat,nan(7,1)); % mem err
                
                % write fit residual
                Gdat_fitresidual=cat(2,Gdat_fitresidual,nan(2,1)); % metric units
                
                % write velocity
                Gdat_velocity=cat(2,Gdat_velocity,nan(3,1));
                
                % write velocity
                Gdat_accelaration=cat(2,Gdat_accelaration,nan(3,1));
                
                % write velocity gradient
                Gdat_velgradient=cat(2,Gdat_velgradient,nan(9,1));
                
                % write strain rate
                Gdat_strainrate=cat(2, Gdat_strainrate,nan(6,1)) ; % symmetric strainrate
                
                % write vorticity vector
                Gdat_vorticity=cat(2,Gdat_vorticity,nan(3,1)); % skew symmetric vorticity
                
                % write divergence
                Gdat_divergence=cat(2,Gdat_divergence,nan); % divergence
                
                % write material derivative
                Gdat_matderivative=cat(2,Gdat_matderivative,nan(3,1)); % material derivative
                
                
            end
            
        end
        
    end
    
    % message
    disp(['processed eulerian reference @ frame ',num2str(n)])
    
end

%% Save

% save eulerian reference grid
save([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'vecgrid'],'vecgrid'); % save vecgrid

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Index.dat'],Gdat_binindex,'w');

% indexing
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Position.dat'],Gdat_position,'w');

% time
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Time.dat'],Gdat_time,'w');

% number density
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'NumberDensity.dat'],Gdat_density,'w');

% binned velocity
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedVelocity.dat'],Gdat_binnedvel,'w');

% linear velocity 
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'LinearVelocity.dat'],Gdat_linearvel,'w');

% polarization vector
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Polarization.dat'],Gdat_polarization,'w');

% binned acc
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedAccelaration.dat'],Gdat_binnedacc,'w');

% displacement field
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'DisplacementField.dat'],Gdat_dispfield,'w');

% % fitted trajectories on grid point
% fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl
% 'TrajectoryFit.dat'],Gdat_fdat,'w'); % memory intrusive

% residual
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'FitResidual.dat'],Gdat_fitresidual,'w');

% velocity vector field
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Velocity.dat'],Gdat_velocity,'w');

% accelaration vectors
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Accelaration.dat'],Gdat_accelaration,'w');

% velocity gradient
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'VelocityGradient.dat'],Gdat_velgradient,'w');

% symmetric strainrate voight vector
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'StrainRate.dat'],Gdat_strainrate,'w');

% % principle strainrate
% fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'PrincipleStrain.dat'],Gdat_princstrain,'w');

% vorticity vector
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Vorticity.dat'],Gdat_vorticity,'w');

% divergence
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Divergence.dat'],Gdat_divergence,'w');

% material derivative
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'MaterialDerivative.dat'],Gdat_matderivative,'w');

% % massflux
% fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'MassFlux.dat'],Gdat_massflux,'w');

% % director polarization
% fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Polarization.dat'],Gdat_dirpol,'w');
%
% % trajectory curvature field
% fsave([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Curvature.dat'],Gdat_cur,'w');

end

