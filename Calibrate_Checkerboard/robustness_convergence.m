%% Info
%
% $Author: Koen Muller$
%

%% Matlab start
close all
clear all
clc

%% Additional Paths
dir= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(dir) %  back to original folder
clear dir

%% Location data
rcs=['D:\Desktop' '\2017_08_07_to_14' '\Oceanium_Calibration' '\Robustness_Convergence']; % robustness and convergence study

%% get data from different runs
% define folders
folder=dir([rcs vsl '*boards']); % folder

% make available markers and colors
smb={'o' 's' 'd' '<' '>' '^' 'v'}; % symbol
col={'r' 'g' 'b' 'k' 'm' }; % color 'y' 'c'
[i,j]=meshgrid(1:length(smb),1:length(col));
mrk=[cell2mat(smb(reshape(i,[],1)))' cell2mat(col(reshape(j,[],1)))'];

% intiate variables
count=0;
Nchk=cell(0);
Smrk=cell(0);
Qfoc=cell(0);
Rerr=cell(0);

% loop folders
for i=1:length(folder)
    subfolder=dir([rcs vsl folder(i).name vsl 'run*']); % subfolder
    
    for j=1:length(subfolder)
        
        % load calibration
        prop=importdata([rcs vsl folder(i).name vsl subfolder(j).name vsl 'prop.mat']);
        filesel=importdata([rcs vsl folder(i).name vsl subfolder(j).name vsl 'filesel.mat']);
        Kmat=importdata([rcs vsl folder(i).name vsl subfolder(j).name vsl 'Kmat.mat']);
        Dpar=importdata([rcs vsl folder(i).name vsl subfolder(j).name vsl 'Dpar.mat']);
        Dmap=importdata([rcs vsl folder(i).name vsl subfolder(j).name vsl 'Dmap.mat']);
        
        % load reprojection error data
        rpe=cell(2,length(filesel));
        for k=1:length(filesel) % new data set indicator
            
            if ~isempty(filesel(k).cal)
                n=unique(filesel(k).cal(:,1));
                c=unique(filesel(k).cal(:,2));
                
                res=fload([rcs vsl folder(i).name vsl subfolder(j).name vsl filesel(k).name '_rpe.dat']);
                res=res(:,ismember(res(1,:),n) & ismember(res(2,:),c));
                
                rpe{1,k}=k*ones(1,size(res,2));
                rpe{2,k}=res;
            end
            
        end
        rpe=cell2mat(rpe);
        
        % read properties experiments
        Npix=prod(prop.res(1:2));
        flens=prop.flens; % focal length lens [mm]
        nref=prop.nref;%0; % relative refractive index (salt)water/air [-]
        Ppix=Npix/prop.areachip; % pix / area chip [mm^2]
        
        % define jacobian of distortion mapping
%         Djac=matlabFunction(reshape(jacobian(sym(Dmap),sym('x',[2 1])),4,[]),'Optimize',false);
%         [y,x]=meshgrid(1:prop.res(2),1:prop.res(1));
%         x=[x(:) y(:)]';
        b=sort(symvar(sym(Dmap)));
        x=subfunv(Dmap,b(1),0); % will not be part of integration
        [~,J]=difobjf(x,sym('x',[2 1]));
        V=intobjf(J,sym('x',[1 2]),[1,prop.res(1);1,prop.res(2)]); % analytic integral
        Djac=divfunh(V,Npix); % gives same result as nummeric integration // tested

        % mark slection experiment
        msel=randi(length(mrk));
        
        % loop cameras
        for c=1:prop.res(4)
            
            % number boards
            C=rpe(3,:)==c;
            num=size(unique(rpe(1:2,C)','rows'),1);
            
            % distortion expansion
%             inp=num2cell([repmat(Dpar{c}',1,size(x,2)) ; x],2);
%             J=Djac(inp{:});
%             J=J(1,:).*J(4,:)-J(2,:).*J(3,:);
%             Jtilde=mean(J(:));
            inp=num2cell(Dpar{c}(2:end));
            Jtilde=Djac(inp{:});
            
            % camera parameters
            focx=Kmat{c}(1,1);
            focy=Kmat{c}(2,2);
            
            % expected focal length from setup
            Fexp=sqrt(focx*focy*Jtilde*nref^2/Ppix); % including magnifying effect of nref
            
            % write variables
            count=count+1;
            Smrk{count}=mrk(msel,:);
            Nchk{count}=num;
            Qfoc{count}=Fexp/flens;
            
            Rerr{1,count}=mean(rpe(12,C),2)*100;
            Rerr{2,count}=std(rpe(12,C),[],2)*100;
        end
        
        % remove that markers
        mrk(msel,:)=[];
        
    end
end
Smrk=cell2mat(Smrk');
Nchk=cell2mat(Nchk);
Qfoc=cell2mat(Qfoc);
Rerr=cell2mat(Rerr);

%% plot robustness
h=figure;
subplot(2,1,1)
h.Position=[460 190 400 325];
hold on
rectangle('Position',[0 0 2 2],'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','none')
plot([2 2],[0 2],'k--','LineWidth',1)
for i=1:length(Smrk)
    plot(Nchk(i),(Qfoc(i)),Smrk(i,:),'MarkerFaceColor',Smrk(i,2),'LineWidth',1)
end
plot([2 50],[1 1],'r--','LineWidth',1)
hold off
xlim([0 50])
ylim([0.0 2])
% xlabel('$N_{\rm{images}} \; [\rm{ \#}]$','interpreter','latex','FontSize',14)
ylabel('$f_{\rm{eff}} / f_{\rm{lens}} [-]$','interpreter','latex','FontSize',14)
h=gca;
h.YTick=[0 : 0.5 : 2];
h.XTick=[2,10:10:50];
box on

%% plot convergence
subplot(2,1,2)
hold on
rectangle('Position',[0 0 2 100],'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','none')
plot([2 50],[2 2],'r--','LineWidth',1)
plot([2 2],[0 100],'k--','LineWidth',1)
errorbar(Nchk,Rerr(1,:),Rerr(2,:),'bsq','LineWidth',1)
hold off
xlim([0 50])
ylim([0 40])
xlabel('$N \; [\rm{ \#}]$','interpreter','latex','FontSize',14)
ylabel('$\varepsilon_j^{n,c} \; [\%]$','interpreter','latex','FontSize',14)
h=gca;
h.YTick=[0:20:100.0];
h.XTick=[2,10:10:50];
box on

%% plot convergence inset
axes('Position',[.38 .25 .5 .18])
hold on
errorbar(Nchk(Nchk>=10),Rerr(1,Nchk>=10),Rerr(2,Nchk>=10),'bsq','LineWidth',1)
plot([15 50],[2 2],'r--','LineWidth',1)
hold off
xlim([10 50])
ylim([0 10])
% xlabel('$N_{\rm{images}} \; [\rm{ \#}]$','interpreter','latex','FontSize',14)
% ylabel('$\bar{\varepsilon} \; [\%]$','interpreter','latex','FontSize',14)
h=gca;
h.YTick=[0,2,5,10];
h.XTick=[10:10:50];
h.FontSize=8;
box on
