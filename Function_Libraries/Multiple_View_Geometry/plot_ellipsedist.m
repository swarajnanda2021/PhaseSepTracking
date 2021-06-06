clear all
close all
clc

% Take conic vector C

C = [1 0 1 0 2 2]';
% ellipse centre
[Xc,ax,ang] = cvec2eshp(C);

ellipsepts = cvec2pnts(C);

[xg,yg] = meshgrid(linspace(min(ellipsepts(1,:)),max(ellipsepts(1,:)),100),linspace(min(ellipsepts(2,:)),max(ellipsepts(2,:)),100));
reproj_err = zeros(size(xg));
for i=1:size(xg,1)
    for j=1:size(xg,2)
        
        reproj_err(i,j) = ellipsedist(C,[xg(i,j) yg(i,j)]');
        
        
    end
end

figure(1)
hold all
[C,h] = contour(xg,yg,reproj_err);
v = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
clabel(C,h,v)
xlabel('x','interpreter','latex','FontSize',20)
ylabel('y','interpreter','latex','FontSize',20)
plot(ellipsepts(1,:),ellipsepts(2,:),'r--')
legend('$\epsilon^*_{reproj}$','ellipse','interpreter','latex','FontSize',20)