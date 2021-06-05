function [camprop]=exp_stats
%exp_stats Export statistics
%
%   Input:
%       Dmap.mat Distortion mapping 
%       Dpar.mat Distortion parameters
%       Kmat.mat Camera matrix 
%       rpe.dat Reprojection error
%       pld.dat Pointline distance
%
%   Output:
%       camprop A summary of the statistics

%% Get global variables
global folder date cal filesel prop

%% load data
pld=cell(2,length(filesel));
rpe=cell(2,length(filesel));
for k=1:length(filesel) % new data set indicator
    
    if ~isempty(filesel(k).dis)
        n=unique(filesel(k).dis(:,1));
        c=unique(filesel(k).dis(:,2));
        
        dis=fload([folder date cal vsl filesel(k).name '_pld.dat']);
        dis=dis(:,ismember(dis(1,:),n) & ismember(dis(2,:),c));
        
        pld{1,k}=k*ones(1,size(dis,2));
        pld{2,k}=dis;
    end
    
    if ~isempty(filesel(k).cal)
        n=unique(filesel(k).cal(:,1));
        c=unique(filesel(k).cal(:,2));
        
        res=fload([folder date cal vsl filesel(k).name '_rpe.dat']);
        res=res(:,ismember(res(1,:),n) & ismember(res(2,:),c));
        
        rpe{1,k}=k*ones(1,size(res,2));
        rpe{2,k}=res;
    end
    
end
rpe=cell2mat(rpe);
pld=cell2mat(pld);

%% load calibration files
% distortion
Dmap=importdata([folder date cal vsl 'Dmap.mat']);
Dpar=importdata([folder date cal vsl 'Dpar.mat']);

% calibration
Kmat=importdata([folder date cal vsl 'Kmat.mat']);
Rmat=importdata([folder date cal vsl 'Rmat.mat']);
tvec=importdata([folder date cal vsl 'tvec.mat']);

%% stats per camera

% define image grid
[y,x]=meshgrid(1:prop.res(2),1:prop.res(1));
x=[x(:)'
    y(:)'];

% initiate
Npix=prod(prop.res(1:2));
Ppix=Npix/prop.areachip;
camprop=struct([]);

% loop views
for c=1:prop.res(4)
    
%     % camera mapping for distortion characteristics
%     load([folder date cal vsl 'wrp_',num2str(c),'.mat'])
    
    % name
    camprop(c).name=['camera',num2str(c)];
    
    % camera parameters
    par=kmat2cpar(Kmat{c});
    camprop(c).focx=par(1);
    camprop(c).focy=par(4);
    camprop(c).skew=par(2);
    camprop(c).cenx=par(3);
    camprop(c).ceny=par(5);
    
    % dewarped coordinates
    inp=num2cell([repmat(Dpar{c}.map',1,size(x,2))
        x],2);
    X=Dmap(inp{:}); % distortion mapping
    X=homc2inhc(Dpar{c}.H*inhc2homc(X)); % back to image resolution
    
    % meshgrid
    Y=reshape(X(2,:),prop.res(1:2));
    X=reshape(X(1,:),prop.res(1:2));
    
    % compute mapping expansion
    Xcoef=SGfilt_coef(X,strel('disk',1,0),[1 1 0],1); % flat, disk, same as central scheme / sobel filter
    Ycoef=SGfilt_coef(Y,strel('disk',1,0),[1 1 0],1); 
    dXdx=SGfilt_eval(Xcoef,[1 0 0]);
    dXdy=SGfilt_eval(Xcoef,[0 1 0]);
    dYdx=SGfilt_eval(Ycoef,[1 0 0]);
    dYdy=SGfilt_eval(Ycoef,[0 1 0]);
    J=dXdx(:).*dYdy(:)-dXdy(:).*dYdx(:); % need to redefine because had been normalized avoiding small divisions
    camprop(c).Jtilde=mean(J(:)); % this is not neccessary the best statistic
    
    % expected focal length from setup
    camprop(c).Fexp=sqrt(camprop(c).focx*camprop(c).focy ... 
        *camprop(c).Jtilde*prop.nref^2/Ppix); % including magnifying effect of nref
    
    % quality focal
    camprop(c).Qfoc=camprop(c).Fexp/prop.flens; % to be checked!! new update
    
    % camera position
    pos=-Rmat{c}\tvec{c};%xcen{c};
    camprop(c).posx=pos(1);
    camprop(c).posy=pos(2);
    camprop(c).posz=pos(3);
    
    % camera orientation
    X=Rmat{c}\eye(3); % inverse are base-vec
    X=[X(:,1) X(:,3) X(:,2)]; % rmat as base vectors
    ang=rotm2eula(X,'ZXY')*180/pi; % modified tait-bryan angles for c-pos.
    camprop(c).alp=ang(1);
    camprop(c).bet=ang(2);
    camprop(c).gam=ang(3);
    
    % number boards
    C=rpe(3,:)==c;
    camprop(c).num=size(unique(rpe(1:2,C)','rows'),1);
    
    % point line distance
    C=pld(3,:)==c;
%     camprop(c).inideflavg=mean(pld(6,C),2);
%     camprop(c).inideflstd=std(pld(6,C),[],2);
    camprop(c).deflavg=mean(pld(8,C),2)*100; % fin
    camprop(c).deflstd=std(pld(8,C),[],2)*100; % in
    
    % reprojection error
    C=rpe(3,:)==c;
%     camprop(c).inirelreprjavg=mean(rpe(10,C),2) ;
%     camprop(c).inirelreprjstd=std(rpe(10,C),[],2);
    camprop(c).relreprjavg=mean(rpe(12,C),2)*100; % in percentage fin
    camprop(c).relreprjstd=std(rpe(12,C),[],2)*100;% fin
    camprop(c).reprjavg=mean(rpe(11,C),2) ;
    camprop(c).reprjstd=std(rpe(11,C),[],2);
    
end
struct2table(camprop);

%% save
save([folder date cal vsl 'camprop.mat'],'camprop')
camprop=[fieldnames(camprop),squeeze(struct2cell(camprop))];
camprop=cell2table([camprop(2:end,1),num2cell(cellfun(@(x)round(round(x,6,'significant'),6),camprop(2:end,2:end)))],'VariableNames',camprop(1,:));
writetable(camprop,[folder date cal vsl 'camprop.csv'],'Delimiter',',','QuoteStrings',true);

%% print
disp(camprop)

end