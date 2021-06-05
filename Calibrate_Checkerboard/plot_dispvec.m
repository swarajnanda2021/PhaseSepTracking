function plot_dispvec
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% Get global variables
global folder date cal filesel chkboard Dmap prop % prop rec

%% load data
IMcnod=cell(2,length(filesel));
qcb=cell(2,length(filesel));
tcb=cell(2,length(filesel));
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
        
    end
end
IMcnod=cell2mat(IMcnod);
qcb=cell2mat(qcb);
tcb=cell2mat(tcb);

%% load calibration files
% distortion
Dmap=importdata([folder date cal vsl 'Dmap.mat']);
Dpar=importdata([folder date cal vsl 'Dpar.mat']);

% calibration
Kmat=importdata([folder date cal vsl 'Kmat.mat']);
Rmat=importdata([folder date cal vsl 'Rmat.mat']);
tvec=importdata([folder date cal vsl 'tvec.mat']);

Pmat=krtm2pmat(Kmat,Rmat,tvec);

%% triangulate and project on views
Xcb=[IMcnod(1:4,:);repmat(chkboard.points,1,size(tcb,2))];
Xt=zeros(7,size(tcb,2)*size(chkboard.points,2));
xp=zeros(6,size(tcb,2)*size(chkboard.points,2));
Xt(1:4,:)=IMcnod(1:4,:);
xp(1:4,:)=IMcnod(1:4,:);
for c=1:prop.res(4)
    % camera set
    C1=qcb(3,:)==c;
    C2=Xt(3,:)==c;
    
    % find rigid body motion checkerboard, per file per board
    for f=unique(qcb(1,C1))
        % file set indicators
        F1=qcb(1,:)==f;
        F2=Xt(1,:)==f;
        
        for n=unique(qcb(2,F1&C1)) %min
            % frame set inidicators
            N1=qcb(2,:)==n;
            N2=Xt(2,:)==n;
            
            % get transformation
            R=qcom2rotm(qcb(4:end,F1&N1&C1));
            t=tcb(4:end,F1&N1&C1);
            
            % triangulate board
            Xt(5:7,F2&N2&C2)=Rmat{c}\(R*inhc2homc(chkboard.points,0)...
                +t-tvec{c});
            
        end
        
    end
    
    % project
    xp(5:6,C2)=homc2inhc(Pmat{c}*inhc2homc(Xt(5:7,C2)));
end

%% plot camera and checkboard disparity
for c=1:prop.res(4)
    % camera set
    C1=IMcnod(3,:)==c;
    C2=tcb(3,:)==c;
    
    % dewarp
    inp=num2cell([repmat(Dpar{c}.map',1,nnz(C1));IMcnod(5:6,C1)],2);
    xw=Dmap(inp{:});
    xw=homc2inhc(Dpar{c}.H*inhc2homc(xw));
    
    % associated disparity
    d=xp(5:6,C1)-xw;
    
    % plot
    figure
    subplot(1,2,1)
    plot(xw(1,:),xw(2,:),'.')
    hold on
    quiver(xw(1,:),xw(2,:),d(1,:),d(2,:),0)
    hold off
    xlabel('x [px]','interpreter','latex')
    ylabel('y [px]','interpreter','latex')
    title(['disparity vectors camera ',num2str(c)] ,'interpreter','latex')
    axis equal
    
    % dewarp camera matrix
    Xw=Kmat{c}\inhc2homc(xw);
    
    % projection homography
    R=cellfun(@(x)qcom2rotm(x),num2cell(qcb(4:end,C2),1),'UniformOutput',false);
    H=cellfun(@(R,t)[R(:,1:2),t],R,num2cell(tcb(4:end,C2),1),'UniformOutput',false);
    Xw=cellfun(@(H,x)H\x,H',squeeze(num2cell(reshape(Xw,3,size(chkboard.points,2),[]),[1 2])),'UniformOutput',false);
    Xw=homc2inhc(cell2mat(Xw'));
    
    % associated disparity
    D=Xw-Xcb(5:6,C1);
    
    % plot
    subplot(1,2,2)
    quiver(Xcb(5,C1),Xcb(6,C1),D(1,:),D(2,:),0)
    hold on
    plot(Xcb(5,C1),Xcb(6,C1),'.')
    hold off
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    xlim([min(chkboard.outline(1,:)) max(chkboard.outline(1,:))])
    ylim([min(chkboard.outline(2,:)) max(chkboard.outline(2,:))])
    title(['disparity vectors checkerboard ',num2str(c)] ,'interpreter','latex')
    axis equal
    
%     figure
%     hold on
%     scatter(xw(1,:),xw(2,:),[],sqrt(sum(d.^2,1)))
%     hold off
    
end

end