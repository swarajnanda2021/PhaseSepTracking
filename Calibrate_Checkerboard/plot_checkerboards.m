function plot_checkerboards
%plot_chktriang Plot the checkerboard traingulation and camera pose for
%   inspection to the obtained camera calibration.
%
%   Input:
%       Kmat.mat Rmat.mat tvec.mat calibration files
%       
%   Output:
%       Checkerboard positions per camera
%       Checkerboard positions in world coordinates

%% Get global variables
global folder date cal prop chkboard filesel

%% load data
IMcnod=cell(2,length(filesel));
qcb=cell(2,length(filesel));
tcb=cell(2,length(filesel));
rpe=cell(2,length(filesel));
for k=1:length(filesel) % new data set indicator
    if ~isempty(filesel(k).cal)
        n=unique(filesel(k).cal(:,1));
        c=unique(filesel(k).cal(:,2));
        
        imcnod=fload([folder date cal vsl filesel(k).name '_IMcnod.dat']);
        imcnod=imcnod(:,ismember(imcnod(1,:),n) & ismember(imcnod(2,:),c));
        
        IMcnod{1,k}=k*ones(1,size(imcnod,2));
        IMcnod{2,k}=imcnod;
        
        q=fload([folder date cal vsl filesel(k).name '_qcb.dat']);
        q=q(:,ismember(q(1,:),n) & ismember(q(2,:),c));
        
        qcb{1,k}=k*ones(1,size(q,2));
        qcb{2,k}=q;
        
        t=fload([folder date cal vsl filesel(k).name '_tcb.dat']);
        t=t(:,ismember(t(1,:),n) & ismember(t(2,:),c));
        
        tcb{1,k}=k*ones(1,size(t,2));
        tcb{2,k}=t;
        
        r=fload([folder date cal vsl filesel(k).name '_rpe.dat']);
        r=r(:,ismember(r(1,:),n) & ismember(r(2,:),c));
        
        rpe{1,k}=k*ones(1,size(r,2));
        rpe{2,k}=r;
    end
end
IMcnod=cell2mat(IMcnod);
qcb=cell2mat(qcb);
tcb=cell2mat(tcb);
rpe=cell2mat(rpe);

%% load calibration files
load([folder date cal vsl 'Kmat.mat'],'Kmat')
load([folder date cal vsl 'Rmat.mat'],'Rmat')
load([folder date cal vsl 'tvec.mat'],'tvec')
Pmat=krtm2pmat(Kmat,Rmat,tvec);

%% plotting
figure(prop.res(4)+1)
hold on
quiver3(0,0,0,1,0,0,0,'k','LineWidth',1,'MaxHeadSize',1)
quiver3(0,0,0,0,1,0,0,'k','LineWidth',1,'MaxHeadSize',1)
quiver3(0,0,0,0,0,1,0,'k','LineWidth',1,'MaxHeadSize',1)
hold off

% get points
Xo=inhc2homc(chkboard.outline,0); %outline

for c=1:prop.res(4)
    
    figure(prop.res(4)+1)
    hold on
    xcen=-Rmat{c}\tvec{c};
    text(xcen(1)+0.5,xcen(2)+0.5,xcen(3)+0.5,num2str(c))
    Xcen=Rmat{c}\eye(3);
    quiver3(xcen(1),xcen(2),xcen(3),Xcen(1,1),Xcen(2,1),Xcen(3,1),0,'r','LineWidth',1,'MaxHeadSize',1)
    quiver3(xcen(1),xcen(2),xcen(3),Xcen(1,2),Xcen(2,2),Xcen(3,2),0,'g','LineWidth',1,'MaxHeadSize',1)
    quiver3(xcen(1),xcen(2),xcen(3),Xcen(1,3),Xcen(2,3)*2,Xcen(3,3),0,'b','LineWidth',1,'MaxHeadSize',1)
    hold off
    
    figure(c)
    hold on
    text(0.5,0.5,0.5,num2str(c))
    quiver3(0,0,0,1,0,0,0,'r','LineWidth',1,'MaxHeadSize',1)
    quiver3(0,0,0,0,2,0,0,'b','LineWidth',1,'MaxHeadSize',1)
    quiver3(0,0,0,0,0,1,0,'g','LineWidth',1,'MaxHeadSize',1)
    hold off
    
    C=qcb(3,:)==c;
    % find rigid body motion checkerboard, per file per board
    for f=unique(qcb(1,C))
        % file set indicators
        F=qcb(1,:)==f;
        
        for n=unique(qcb(2,F&C)) %min
            % frame set inidicators
            N=qcb(2,:)==n;
            
            % get transformation
            R=qcom2rotm(qcb(4:end,F&N&C));
            t=tcb(4:end,F&N&C);
            
            % triangulate board
            X=R*Xo+t;
            
            % plotting
            figure(c)
            hold on
            patch(X(1,:),X(3,:),X(2,:),zeros(size(X(3,:))),'FaceColor',[0.5 0.5 1],'FaceAlpha',0.2,'EdgeColor',[0.7 0.7 1])%,'LineWidth',0.5)
            hold off
            
            R=Rmat{c}\R;
            t=Rmat{c}\(t-tvec{c});
            
            % triangulate board
            X=R*Xo+t;
            
            % plotting
            figure(prop.res(4)+1)
            hold on
            patch(X(1,:),X(2,:),X(3,:),zeros(size(X(3,:))),'FaceColor',[0.5 0.5 1],'FaceAlpha',0.2,'EdgeColor',[0.7 0.7 1])%,'LineWidth',0.5)
            hold off
            
%             % get points
%             X=inhc2homc(chkboard.points,0);
            
%             % triangulate board
%             X=R*X+t;
%             
%             %plotting
%             plot3(X(1,:),X(2,:),X(3,:),'.','Color',[1 1 1]/2,'LineWidth',0.25)
        end
    end
    figure(c)
    axis tight
    axis equal
    grid minor
    box on
    ax=gca;
    ax.XLim=[-4 4];
    ax.YLim=[-1 25];
    ax.ZLim=[-4 4];
    ax.XTick =-4:4:4;
    ax.YTick =-0:5:25;
    ax.ZTick =-4:4:4;
    ax.Clipping= 'off';
    % ax.LineWidth=1;
    xlabel('x $[\rm{m}]$','interpreter','latex','FontSize',14)
    ylabel('k $[\rm{m}]$','interpreter','latex','FontSize',14)
    zlabel('y $[\rm{m}]$','interpreter','latex','FontSize',14)
    camproj('perspective')
    ax.CameraPosition = [20 -10 4];
    title(['Checkerboards camera ',num2str(c)],'interpreter','latex')
    
%     axis off
    
end
figure(prop.res(4)+1)
axis tight
axis equal
grid minor
box on
ax=gca;
ax.XLim=[-5 5];
ax.YLim=[-1 25];
ax.ZLim=[-5 2];
ax.XTick =-5:5:5;
ax.YTick =-0:5:25;
ax.ZTick =-5:5:0;
ax.Clipping= 'off';
% ax.LineWidth=1;
xlabel('x $[\rm{m}]$','interpreter','latex','FontSize',14)
ylabel('y $[\rm{m}]$','interpreter','latex','FontSize',14)
zlabel('z $[\rm{m}]$','interpreter','latex','FontSize',14)
camproj('perspective')
ax.CameraPosition = [20 -10 2];
title('Checkerboards and camera pose object domain','interpreter','latex')
drawnow

% view(2)
% camproj('orthographic')

end