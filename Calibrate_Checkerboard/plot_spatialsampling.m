function plot_spatialsampling
%plot_spatialsampling Plot the spatial sampling of checkerboard over the
%   object domain
%
%   Input:
%       chkboard.mat Checkerboard
%       qcb.dat & tcb.dat Checkerboard positioning
%
%   Output:
%       Sampling over depth of field
%
%   Note: The sampling variables can be tuned internally in this script.

%% get global variables
global folder date cal filesel chkboard 

%% load image and calibration data
IMcnod=cell(2,length(filesel));
qcb=cell(2,length(filesel));
tcb=cell(2,length(filesel));
for k=1:length(filesel) % loop file selection
    
    if ~isempty(filesel(k).cal) % if selected
        % get nodes
        imcnod=fload([folder date cal vsl filesel(k).name '_IMcnod.dat']);
        imcnod=imcnod(:, ismember(imcnod([1,2],:)',filesel(k).cal,'rows') );
        
        % write nodes
        IMcnod{2,k}=imcnod;
        IMcnod{1,k}=k*ones(1,size(IMcnod{2,k},2));
        
        % get unit quaternions
        q=fload([folder date cal vsl filesel(k).name '_qcb.dat']);
        q=q(:, ismember(q([1,2],:)',filesel(k).cal,'rows') );
        
        % write unit quaternions
        qcb{1,k}=k*ones(1,size(q,2));
        qcb{2,k}=q;
        
        % get translations
        t=fload([folder date cal vsl filesel(k).name '_tcb.dat']);
        t=t(:, ismember(t([1,2],:)',filesel(k).cal,'rows') );
        
        % write translations
        tcb{1,k}=k*ones(1,size(t,2));
        tcb{2,k}=t;
        
    end
end
IMcnod=cell2mat(IMcnod);
qcb=cell2mat(qcb);
tcb=cell2mat(tcb);

%% load calibration files
Rmat=importdata([folder date cal vsl 'Rmat.mat'],'Rmat'); % rotation matrices
tvec=importdata([folder date cal vsl 'tvec.mat'],'tvec'); % translations vectors
Pmat=cell(1,4); % initiate projection matrices
for i=1:size(Pmat,2) % construct projection matrices
    Pmat{i}=krtm2pmat(eye(3),Rmat{i},tvec{i});
end

%% Position checkerboards in object space
% triangulate same order as IMcnod such that they correspond
R=cellfun(@(x,y)Rmat{y}\qcom2rotm(x),...
    num2cell(qcb(4:end,:),1),num2cell(tcb(3,:),1),'UniformOutput',false); % rotations
t=cellfun(@(x,y)Rmat{y}\(x-tvec{y}),...
    num2cell(tcb(4:end,:),1),num2cell(tcb(3,:),1),'UniformOutput',false); % translations
Xchk=cell2mat(cellfun(@(x,y)x*inhc2homc(chkboard.points,0)+y,R,t,'UniformOutput',false)); % place

% unique checkerboards
[~,ia,~]=unique(IMcnod([1,2,4],:)','rows','stable');
ia=ia'; % use to get unique checkerboards points
Xchk=Xchk(:,ia); % keep unique

%% Sampling over depth of tank
 % sample depth of field
dep=linspace(0,25,26);

% loop
I = discretize(Xchk(2,:),dep);
num_chk=ones(size(dep));
for i=1:length(dep)
    sel=I==i;
    
    num_chk(i)=nnz(sel)/size(chkboard.points,2); % yet correct for num view
end

%make figure
h=figure;
% h.Color=[1 1 1];
% h.Position=[300 300 500 200];
plot(dep(num_chk==0),num_chk(num_chk==0),'r.')
hold on
plot(dep(num_chk~=0),num_chk(num_chk~=0),'b--')
plot(dep(num_chk~=0),num_chk(num_chk~=0),'b.')
hold off
xlim([0 25])
% ylim([0 30])
% ax = gca;
% ax.YTick =0:5:30;
% ax.XTick =0:5:25;
% ax.LineWidth=1;
box on
xlabel('depth $Y \ [\rm{m}]$','interpreter','latex','FontSize',14)
ylabel('sampling $N \ [\rm{\#}]$','interpreter','latex','FontSize',14)

%% Slices over depth of field
% Sampling depth of field and XZ-plane
dep=[4 17 25];
D = discretize(Xchk(2,:),dep);

% limits
Xmin=-4;
Xmax=3;
Xdel=0.5;
Zmin=-5;
Zmax=0;
Zdel=0.5;

% make edges
Xedg=Xmin:Xdel:Xmax;
Zedg=Zmin:Zdel:Zmax;

% make grid
Xmap=Xmin+Xdel/2:Xdel:Xmax-Xdel/2;
Zmap=Zmin+Zdel/2:Zdel:Zmax-Zdel/2;
[Zmap,Xmap]=meshgrid(Zmap,Xmap);

% loop over depth
for k=unique(D(~isnan(D)))
    % selection
    K=D==k;
    
    % use histcounts2 to find binning for data
    [Smap]=histcounts2(Xchk(1,K),Xchk(3,K),Xedg,Zedg);
    
    % figure sampling
    h=figure;
%     h.Color=[1 1 1];
%     h.Position=[300 300 250 200];
    contourf(Xmap',Zmap',Smap',10)
    view(2)
    box on
    shading flat
    axis equal
%     caxis([0 90])
    colorbar
    colormap parula
    xlim([Xmin Xmax])
    ylim([Zmin Zmax])
    xlabel('$X \ [\rm{m}]$','interpreter','latex')
    ylabel('$Z \ [\rm{m}]$','interpreter','latex')
    ax = gca;
    ax.XTick =-4:3;
    ax.YTick =-5:1;
    ax.LineWidth=0.5;
    drawnow
    
end

% % generate separate colorbar
% h=figure;
% h.Color=[1 1 1];
% h.Position=[300 300 250 200];
% axis off
% caxis([0 90])
% c=colorbar('northoutside');
% c.Ticks=[0 90];
% ylabel(c,'sampling $N \ [\rm{\#}]$','interpreter','latex','FontSize',12)

end
