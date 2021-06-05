function [prop,ctrl] = proc_set
%proc_set Processing settings
% Output,
%   prop - Properties
%   ctrl - Control

%% global
global folder date rec

%% Specific
prop.ext='im7';
[vid,prop.res]=import_frames({folder date rec},prop.ext,1,1);
dum=class(vid);
prop.bit=str2double(dum(5:end));

ctrl.tspan=[1 prop.res(3)];

end