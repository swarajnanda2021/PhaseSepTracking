function postproc_lagrangianreference
%postproc_lagrangianreference
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global folder date rec lagr

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Pos=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Vel=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
!etc

%% make folder for data
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'LagrangianReference'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'LagrangianReference'])
end

%% change reference frame by groupframe


%% save

fsave([folder date rec vsl 'Postproc3DTracking' vsl 'LagrangianFrame' vsl 'Something.dat'],Something,'w');

end
