function [ chkboard ] = def_chkboard
%def_chkboard Define checkerboard
%   
%   Input: processing setting
%
%	Output:
%       chkboard.mat Checkerboard parametrization

%% get globals
global folder date cal ctrl

%% define checkerboard
[y,x]=meshgrid(0:ctrl.nnod(2)-1,0:ctrl.nnod(1)-1);
chkboard.points=([reshape(x,[],1),reshape(y,[],1)]'*ctrl.tsiz); % consistent with id number
chkboard.outline=[min(x(:))-1 min(y(:))-1
    max(x(:))+1 min(y(:))-1
    max(x(:))+1 max(y(:))+1
    min(x(:))-1 max(y(:))+1
    min(x(:))-1 min(y(:))-1]'*ctrl.tsiz;

%% Make figure
figure
plot(chkboard.points(1,:),chkboard.points(2,:),'.r')
hold on
for i=1:size(chkboard.points,2)
    text(chkboard.points(1,i)+ctrl.tsiz/8,chkboard.points(2,i)+ctrl.tsiz/8,num2str(i))
end
plot(chkboard.outline(1,:),chkboard.outline(2,:),'-b','LineWidth',1)
hold off
axis equal
xlim([min(chkboard.outline(1,:))-ctrl.tsiz/2 max(chkboard.outline(1,:))+ctrl.tsiz/2])
ylim([min(chkboard.outline(2,:))-ctrl.tsiz/2 max(chkboard.outline(2,:))+ctrl.tsiz/2])
xlabel('$X^o$','interpreter','latex')
ylabel('$Y^o$','interpreter','latex')
title('Checkerboard Nodes and Outline','interpreter','latex')

%% save
save([folder date cal vsl 'chkboard.mat'],'chkboard')

end