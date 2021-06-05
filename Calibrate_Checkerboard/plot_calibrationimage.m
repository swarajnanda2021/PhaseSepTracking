function plot_calibrationimage(f,c,n)
%plot_calimage Plot calibration image to interest and visualize the
%   different calibration steps
%
%   Input:
%       f file index
%       c camera index
%       n frame index
%
%       IMblin.dat Image baselines checkerboard
%       DWblin.dat Dewarped baselines checkerboard
%       IMcfit.dat Checkerboard curvefit gridlines
%       IMcnod.dat Checkerboard nodes
%       Kmat.mat wrp.mat
%
%   Output:
%       Plot desired image

%% get globals
global folder date cal prop Dmap Cfit filesel ctrl

%% load data
IMblin=fload([folder date cal vsl filesel(f).name '_IMblin.dat']);
DWblin=fload([folder date cal vsl filesel(f).name '_DWblin.dat']);
IMcfit=fload([folder date cal vsl filesel(f).name '_IMcfit.dat']);
IMcnod=fload([folder date cal vsl filesel(f).name '_IMcnod.dat']);

%% load calibration
Dpar=importdata([folder date cal vsl 'Dpar.mat']);
Kmat=importdata([folder date cal vsl 'Kmat.mat']);
wrp=importdata([folder date cal vsl 'wrp_',num2str(c),'.mat']);

%% image data
imd=double(import_frames({folder date cal [vsl filesel(f).name]},prop.ext,n,c)); % board scatters

%% plot camera image
h=figure;
h.Color=[0 0 0];
% h.Position=[300 300 300 300];
% axes('Units', 'normalized', 'Position', [0 0 1 1])
plt=log(imd+1);
surf(wrp.x',wrp.y',zeros(size(plt')),plt')
    caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
view(2)
shading flat
axis off
axis equal
colormap gray
drawnow

%% get lines data
C=IMblin(2,:)==c;
N=IMblin(1,:)==n;
lcb=IMblin(4:end,C&N);

% get outer perimeter board
xboun=[1 prop.res(1) prop.res(1) 1
    1 1 prop.res(2) prop.res(2)];
pcb=lvec2pnts(lcb,xboun);

%% plot camera image
h=figure;
h.Color=[0 0 0];
% h.Position=[300 300 300 300];
% axes('Units', 'normalized', 'Position', [0 0 1 1])
surf(wrp.x',wrp.y',zeros(size(imd')),log(imd'+1))
hold on
for i=1:sum(ctrl.nnod)
    plot([pcb(1,i) pcb(3,i)],[pcb(2,i) pcb(4,i)],'r','Linewidth',1)
    if i<=ctrl.nnod(1)
        text(pcb(1,i)+20,pcb(2,i)+70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
    else
        text(pcb(1,i)-70,pcb(2,i)+70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
    end
end
hold off
view(2)
shading flat
axis off
axis equal
colormap gray
drawnow

%% get curve and points data
C=IMcfit(2,:)==c;
N=IMcfit(1,:)==n;
acb=IMcfit(4:end,C&N);
t=linspace(0,1,100);

C=IMcnod(2,:)==c;
N=IMcnod(1,:)==n;
xcb=IMcnod(4:end,C&N);

%% plot camera image
h=figure;
h.Color=[0 0 0];
% h.Position=[300 300 300 300];
% axes('Units', 'normalized', 'Position', [0 0 1 1])
plt=log(imd+1);
surf(wrp.x',wrp.y',zeros(size(plt')),plt')
    caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
hold on
% for i=1:sum(ctrl.nnod)
%     plot([pcb(1,i) pcb(3,i)],[pcb(2,i) pcb(4,i)],'r','Linewidth',1)
%     if i<=ctrl.nnod(1)
%         text(pcb(1,i)+20,pcb(2,i)+70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     else
%         text(pcb(1,i)-70,pcb(2,i)+70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     end
% end
for i=1:size(acb,2)
    inp=num2cell([repmat(acb(:,i),1,length(t));t],2);
    xlin=Cfit(inp{:});
    plot(xlin(1,:),xlin(2,:),'g-','LineWidth',1)
    plot(xlin(1,:),xlin(2,:),'g-','LineWidth',1)
%     if i<=ctrl.nnod(1)
%         text(xlin(1,1)-20,xlin(2,1)-70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     else
%         text(xlin(1,1)-100,xlin(2,1)+10,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     end
end
plot(xcb(1,:),xcb(2,:),'rsq','LineWidth',1)
for i=1:size(xcb,2)
    text(xcb(1,i)+20,xcb(2,i)+100,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)%,'LineWidth',1)
end
hold off
view(2)
shading flat
axis off
axis equal
colormap gray
drawnow

%% dewarp distortion
[xw,yw]=imtrans(wrp.X,wrp.Y,[],inv(Kmat{c}));

% dew lines
C=DWblin(2,:)==c;
N=DWblin(1,:)==n;
ldw=DWblin(4:end,C&N); % lines dewarped

% get outer perimeter board
pcb=zeros(4,size(ldw,2));
for i=1:sum(ctrl.nnod)
    if i<=ctrl.nnod(1)
        xmin=homc2inhc(cross(ldw(:,i),cross([min(xw(:));min(yw(:));1],[max(xw(:));min(yw(:));1]))); % xmin
        xmax=homc2inhc(cross(ldw(:,i),cross([min(xw(:));max(yw(:));1],[max(xw(:));max(yw(:));1]))); %xmax
    else
        xmax=homc2inhc(cross(ldw(:,i),cross([min(xw(:));min(yw(:));1],[min(xw(:));max(yw(:));1]))); % xmin
        xmin=homc2inhc(cross(ldw(:,i),cross([max(xw(:));min(yw(:));1],[max(xw(:));max(yw(:));1]))); %xmax
    end
%     xcen=(xmin+xmax)/2;
%     xdel=xmax-xmin;
    pcb(1:2,i)=xmin;%xcen-xdel/2;
    pcb(3:4,i)=xmax;%xcen+xdel/2;
end

%% plot dewarped image
h=figure;
h.Color=[0 0 0];
% h.Position=[300 300 300 300];
% axes('Units', 'normalized', 'Position', [0 0 1 1])
plt=log((imd./wrp.J)+1);
surf(xw',yw',zeros(size(plt')),plt')
    caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
hold on
% for i=1:sum(ctrl.nnod)
%     plot([pcb(1,i) pcb(3,i)],[pcb(2,i) pcb(4,i)],'r','Linewidth',1)
%     if i<=ctrl.nnod(1)
%         text(pcb(1,i)+20,pcb(2,i)+70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     else
%         text(pcb(1,i)-70,pcb(2,i)+70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     end
% end
for i=1:size(acb,2)
    inp=num2cell([repmat(acb(:,i),1,length(t));t],2);
    xlin=Cfit(inp{:});
    inp=num2cell([repmat(Dpar{c}.map',1,length(t));xlin],2);
    xlin=Dmap(inp{:});
    xlin=homc2inhc(Dpar{c}.H*inhc2homc(xlin));
    plot(xlin(1,:),xlin(2,:),'g-','LineWidth',1)
%     if i<=ctrl.nnod(1)
%         text(xlin(1,1)-20,xlin(2,1)-70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     else
%         text(xlin(1,1)-100,xlin(2,1)+10,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     end
end
inp=num2cell([repmat(Dpar{c}.map',1,size(xcb,2));xcb],2);
xpnt=Dmap(inp{:});
xpnt=homc2inhc(Dpar{c}.H*inhc2homc(xpnt));
plot(xpnt(1,:),xpnt(2,:),'rsq','LineWidth',1)
for i=1:size(xcb,2)
    text(xpnt(1,i)-20,xpnt(2,i)+100,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)%,'LineWidth',1)
end
hold off
view(2)
shading flat
axis off
axis equal
colormap gray
drawnow

%% plot map to ref pln
h=figure;
h.Color=[0 0 0];
% h.Position=[300 300 300 300];
% axes('Units', 'normalized', 'Position', [0 0 1 1])
plt=log((imd./wrp.J)+1);
surf(wrp.X',wrp.Y',zeros(size(plt')),plt')
    caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
    hold on
for i=1:size(acb,2)
    inp=num2cell([repmat(acb(:,i),1,length(t));t],2);
    xlin=Cfit(inp{:});
    inp=num2cell([repmat(Dpar{c}.map',1,length(t));xlin],2);
    xlin=Dmap(inp{:});
    xlin=homc2inhc(Dpar{c}.H*inhc2homc(xlin));
    xlin=homc2inhc(Kmat{c}\inhc2homc(xlin));
    plot(xlin(1,:),xlin(2,:),'g-','LineWidth',1)
%     if i<=ctrl.nnod(1)
%         text(xlin(1,1)-20,xlin(2,1)-70,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     else
%         text(xlin(1,1)-100,xlin(2,1)+10,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)
%     end
end
inp=num2cell([repmat(Dpar{c}.map',1,size(xcb,2));xcb],2);
xpnt=Dmap(inp{:});
xpnt=homc2inhc(Dpar{c}.H*inhc2homc(xpnt));
xpnt=homc2inhc(Kmat{c}\inhc2homc(xpnt));
plot(xpnt(1,:),xpnt(2,:),'rsq','LineWidth',1)
for i=1:size(xcb,2)
    text(xcb(1,i)-20,xcb(2,i)+100,num2str(i),'Color',[1 1 1]-eps,'FontSize',14)%,'LineWidth',1)
end
    
view(2)
shading flat
axis equal
colormap gray
ax = gca;
ax.Color='none';
ax.Layer='top';
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.XColor='white';
ax.YColor='white';
ax.GridLineStyle='none';
ax.MinorGridLineStyle='none';
ax.YTick =-1:0.05:1;
ax.XTick =-1:0.05:1; % 0.05
ax.LineWidth=1;

xlabel('x $[\rm{m}]$','interpreter','latex','Position',[0.15 -0.0125])
ylabel('y $[\rm{m}]$','interpreter','latex','Position',[-0.0125 0.15])
drawnow

end