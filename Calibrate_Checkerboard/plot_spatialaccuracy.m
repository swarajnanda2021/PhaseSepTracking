function plot_spatialaccuracy
%plot_spatialaccuracy Plot spatialaccuracy of the checkerboard and point
%   triangulations over the object domain
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
global folder date cal filesel chkboard Dmap

%% load image and calibration data
IMcnod=cell(2,length(filesel));
qcb=cell(2,length(filesel));
tcb=cell(2,length(filesel));
for k=1:length(filesel) % loop file selection
    
    if ~isempty(filesel(k).cal) % if selected
        % select only nodes imaged from minimal two views
        [un,~,in]=unique(filesel(k).cal(:,1));
        sel=histcounts(filesel(k).cal(:,1),[un-1/2;max(un)+1/2])'>1;
        sel=sel(in);
        
        % get nodes
        imcnod=fload([folder date cal vsl filesel(k).name '_IMcnod.dat']);
        imcnod=imcnod(:, ismember(imcnod([1,2],:)',filesel(k).cal(sel,:),'rows') );
        
        % write nodes
        IMcnod{2,k}=imcnod;
        IMcnod{1,k}=k*ones(1,size(IMcnod{2,k},2));
        
        % get unit quaternions
        q=fload([folder date cal vsl filesel(k).name '_qcb.dat']);
        q=q(:, ismember(q([1,2],:)',filesel(k).cal(sel,:),'rows') );
        
        % write unit quaternions
        qcb{1,k}=k*ones(1,size(q,2));
        qcb{2,k}=q;
        
        % get translations
        t=fload([folder date cal vsl filesel(k).name '_tcb.dat']);
        t=t(:, ismember(t([1,2],:)',filesel(k).cal(sel,:),'rows') );
        
        % write translations
        tcb{1,k}=k*ones(1,size(t,2));
        tcb{2,k}=t;
        
    end
end
IMcnod=cell2mat(IMcnod);
qcb=cell2mat(qcb);
tcb=cell2mat(tcb);

%% load calibration files
Dpar=importdata([folder date cal vsl 'Dpar.mat'],'Dpar'); % distortion parameters
Kmat=importdata([folder date cal vsl 'Kmat.mat'],'Kmat'); % calibration matrix
Rmat=importdata([folder date cal vsl 'Rmat.mat'],'Rmat'); % rotation matrices
tvec=importdata([folder date cal vsl 'tvec.mat'],'tvec'); % translations vectors
Pmat=cell(1,4); % initiate projection matrices
for i=1:size(Pmat,2) % construct projection matrices
    Pmat{i}=krtm2pmat(eye(3),Rmat{i},tvec{i});
end

%% initial triangulation checkerboards
% triangulate same order as IMcnod such that they correspond
R=cellfun(@(x,y)Rmat{y}\qcom2rotm(x),...
    num2cell(qcb(4:end,:),1),num2cell(tcb(3,:),1),'UniformOutput',false); % rotations
t=cellfun(@(x,y)Rmat{y}\(x-tvec{y}),...
    num2cell(tcb(4:end,:),1),num2cell(tcb(3,:),1),'UniformOutput',false); % translations
Xchk=cell2mat(cellfun(@(x,y)x*inhc2homc(chkboard.points,0)+y,R,t,'UniformOutput',false)); % place

% make figure of checkerboard positioning
% figure
% plot3(Xchk(1,:),Xchk(2,:),Xchk(3,:),'.') % aprox
% axis equal
% grid minor
% colorbar
% box on
% ax=gca;
% ax.XLim=[-5 5];
% ax.YLim=[-1 25];
% ax.ZLim=[-5 2];
% ax.XTick =-5:5:5;
% ax.YTick =-0:5:25;
% ax.ZTick =-5:5:0;
% ax.Clipping= 'off';
% xlabel('x $[\rm{m}]$','interpreter','latex','FontSize',14)
% ylabel('y $[\rm{m}]$','interpreter','latex','FontSize',14)
% zlabel('z $[\rm{m}]$','interpreter','latex','FontSize',14)
% camproj('perspective')
% ax.CameraPosition = [20 -10 2];
% drawnow

%% Object triangulation checkerboard points
% unique checkerboards
[~,ia,ic]=unique(IMcnod([1,2,4],:)','rows','stable');%
ia=ia'; % use to get unique checkerboards points
ic=ic'; % use to enumerate unique checkerboard points

% dewarp
Ddat=cell2mat(cellfun(@(x)reshape(Dpar{x}.map,[],1),num2cell(IMcnod(3,:),1),'UniformOutput',false));
x=IMcnod(5:6,:);
inp=num2cell([Ddat;x],2);
x=Dmap(inp{:});
x=cell2mat(cellfun(@(x,y)homc2inhc(Dpar{y}.H*inhc2homc(x)),num2cell(x,1),num2cell(IMcnod(3,:),1),'UniformOutput',false));

% Assemble triangulation data
xdat=cell(1,length(Pmat));
for c=1:length(Pmat)
    xdat{c}=nan*zeros(3,length(unique(ic))); %0
    xdat{c}(:,ic(IMcnod(3,:)==c))=Kmat{c}\inhc2homc(x(:,IMcnod(3,:)==c));
end

% compute triangulation
[Xt,~,r]=objtriang(xdat,Pmat);

% accuracy in centimeters
acc=r*100;

% make figure point triangulation
figure
scatter3(Xt(1,:),Xt(2,:),Xt(3,:),[],acc,'.') % aprox
axis equal
grid minor
colorbar
box on
ax=gca;
ax.XLim=[-5 5];
ax.YLim=[-1 25];
ax.ZLim=[-5 2];
ax.XTick =-5:5:5;
ax.YTick =-0:5:25;
ax.ZTick =-5:5:0;
ax.Clipping= 'off';
xlabel('x $[\rm{m}]$','interpreter','latex','FontSize',14)
ylabel('y $[\rm{m}]$','interpreter','latex','FontSize',14)
zlabel('z $[\rm{m}]$','interpreter','latex','FontSize',14)
camproj('perspective')
ax.CameraPosition = [20 -10 2];
drawnow

%% accuracy over depth of tank
% Sample depth of field
dep=linspace(0,25,26);
I = discretize(Xt(2,:),dep);
acc_avg=ones(size(dep));
acc_std=ones(size(dep));
for i=1:length(dep)
    sel=I==i;
    
    acc_avg(i)=mean(acc(sel));
    acc_std(i)=std(acc(sel));
    
end

% make figure
h=figure;
% h.Color=[1 1 1];
% h.Position=[300 300 500 200];
hold on
% rectangle('Position',[0 0 4 1],'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','none')
% plot([4 4],[0 1],'k--','LineWidth',1)
% rectangle('Position',[23 0 2 1],'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','none')
% plot([23 23],[0 1],'k--','LineWidth',1)
errorbar(dep,acc_avg,acc_std,'-s','Linewidth',1,'Color','blue','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
hold off
xlim([0 25])
% ylim([0 1])
% ax = gca;
% ax.YTick =0:0.25:1;
% ax.XTick =0:5:25;
% ax.LineWidth=1;
box on
xlabel('depth $Y \ [\rm{m}]$','interpreter','latex','FontSize',14)
ylabel('skew $s_j^n \ [\rm{cm}]$','interpreter','latex','FontSize',14)

%% slices over depth of field
% Sample depth of field and XZ-plane
dep=[4 17 25];
D = discretize(Xt(2,:),dep);

% limits grid
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

% loop depth 
for k=unique(D(~isnan(D)))
    % selection
    K=D==k;
    
    % position data and accuracy
    xpos=Xt([1,3],K);
    apos=acc(K);
    
    % use histcounts2 to find binning for data
    [Smap,~,~,binX,binY]=histcounts2(xpos(1,:),xpos(2,:),Xedg,Zedg);
    
    % bin data
    AVGmap=nan*zeros(size(Smap));
    STDmap=nan*zeros(size(Smap));
    for i=unique(binX)
        I=binX==i;
        
        for j=unique(binY)
            J=binY==j;
            
            if i>0 && j>0 && nnz(I&J)>0
                AVGmap(i,j)=nanmean(apos(I&J));
                STDmap(i,j)=nanstd(apos(I&J));
            end
        end
    end
    
    % figure  error
    h=figure;
%     h.Color=[1 1 1];
%     h.Position=[300 300 250 200];
    contourf(Xmap',Zmap',AVGmap',10)
    view(2)
    box on
    shading flat
    axis equal
%     caxis([0 0.75])
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

% % generate colorbars
% h=figure;
% h.Color=[1 1 1];
% h.Position=[300 300 250 200];
% axis off
% c=colorbar('northoutside');
% c.Ticks=[0 1];
% ylabel(c,'skew $s_j^n \ [\rm{cm}]$','interpreter','latex','FontSize',12)

end
