function [prop,ctrl] = proc_set
%proc_set Processing settings
% Output,
%   prop - Properties
%   ctrl - Control

%% global
global folder date cal

%% properties
% first cam
prop.ext='im7';
[vid,prop.res]=import_frames({folder date [cal vsl 'Calibration' vsl 'camera1']},prop.ext,1,1);
dum=class(vid);
prop.bit=str2double(dum(5:end));

% how many cams?
cams=dir([folder date cal vsl 'Calibration' vsl 'camera*.set']); % by set file
count=0;
for k=1:length(cams)
    dum=dir([folder date cal vsl 'Calibration' vsl cams(k).name(1:end-4)]);
    if size(dum,1)>2
        count=count+1;
    end
end
prop.res(4)=count;

%% controls
ctrl.fobj=strel('disk',3,0); % or better use fwin
ctrl.optl=10^-6;
ctrl.mord=3; % mapping calibration order>1 for usage else ==1 and <1

%% save proc set
disp('save properties and controls')
tic
save([folder date cal vsl 'prop'],'prop'); % save properties
save([folder date cal vsl 'ctrl'],'ctrl'); % save controls
toc

end