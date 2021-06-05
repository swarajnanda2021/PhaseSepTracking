function plot_errordistribution
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
load([folder date cal vsl 'Kmat.mat'],'Kmat')
load([folder date cal vsl 'Rmat.mat'],'Rmat')
load([folder date cal vsl 'tvec.mat'],'tvec')

%% figures distribution 
n=20; % number bins

%% histogram relative reprojection error
edges=linspace(0,max(rpe(end,rpe(end,:)<0.10))*100,n);
bins=edges(2:end)-diff(edges)/2;
counts=zeros(prop.res(4),length(bins));
for c=1:prop.res(4)
    counts(c,:)=histcounts(rpe(end,rpe(3,:)==3 & rpe(end,:)<0.10)*100,edges);
end

h=figure;
h.Position=[460 190 325 225];
h=bar(bins',counts','stacked');
col=[1 0 0
    0 1 0
    0 0 1
    1 0 1];
for c=1:prop.res(4)
    h(c).FaceColor=col(c,:);
    h(c).EdgeColor=col(c,:);
end
xlim([0 0.1*100])
ylim([0 4000])
h=gca;
h.XTick=0:0.02*100:0.1*100;
h.YTick=0:1000:4000;
xlabel('$\bar{\varepsilon} \; [\mathrm{\%}]$','interpreter','latex','FontSize',14)
ylabel('$N \times J \; [\#]$','interpreter','latex','FontSize',14)
box on
legend('cam1','cam2','cam3','cam4')

%% histogram relative percentage of distortion
edges=linspace(0,max(pld(end,pld(end,:)<0.005))*100,n);
bins=edges(2:end)-diff(edges)/2;
counts=zeros(prop.res(4),length(bins));
for c=1:prop.res(4)
    counts(c,:)=histcounts(pld(end,pld(3,:)==c & pld(end,:)<0.005)*100,edges);
end

h=figure;
h.Position=[460 190 325 225];
h=bar(bins',counts','stacked');
col=[1 0 0
    0 1 0
    0 0 1
    1 0 1];
for c=1:prop.res(4)
    h(c).FaceColor=col(c,:);
    h(c).EdgeColor=col(c,:);
end
xlim([0 0.005*100])
ylim([0 3000])
h=gca;
h.XTick=0:0.001*100:0.005*100;
h.YTick=0:750:3000;
xlabel('$\partial \bar{L} \; [\mathrm{\%}]$','interpreter','latex','FontSize',14)
ylabel('$N \times I \; [\#]$','interpreter','latex','FontSize',14)
box on
legend('cam1','cam2','cam3','cam4')


end

