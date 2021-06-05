%% Notes
%
%   Author: Koen Muller, 27 May 2020
%
% example script for spatial 2-point corrolation

%% Start
close all
clear all
clc

%% Controls
% add here, skiped for now

%% create dummy data

% some grid
[X,Y,Z]=ndgrid(-3*pi:pi/10:3*pi,-3*pi:pi/10:3*pi,-3*pi:pi/10:3*pi);%1:0.5:10 % aliasing will make a difference in the result

% indexed grid
I=1:numel(X);
grd=[I
    X(:)'
    Y(:)'
    Z(:)'];

% some mask for reference data (incompleteness)
M=(X.^2+Y.^2+Z.^2+2^2)<=5.^2;

% some peaks
P=reshape(cos(X).*cos(Y).*cos(Z/2).*M,1,[]);

% write dat reference
M=reshape(M,1,[]);
dat1=[I(M)
    P(M)];

% some mask for e.g. displaced or another data to cross corrolate
M=(X.^2+Y.^2+Z.^2)<=5.^2;

% some peaks
P=reshape(sin(X).*cos(Y).*cos(Z/2).*M,1,[]);

% write dat second
M=reshape(M,1,[]);
dat2=[I(M)
    P(M)];

%% plot data

figure
scatter3(grd(2,dat1(1,:)),grd(3,dat1(1,:)),grd(4,dat1(1,:)),100,squeeze(dat1(2,:)),'.')
colorbar

figure
scatter3(grd(2,dat2(1,:)),grd(3,dat2(1,:)),grd(4,dat2(1,:)),100,squeeze(dat2(2,:)),'.')
colorbar

%% evaluate direct spatial 2-point corrolation

% evaluate
[dis,cor,wav,spe,nor] = spatial2pointcor(grd,dat1,dat2,'direct');

% remove empty corrolations
sel=~isinf(nor);
dis=dis(:,sel);
cor=cor(:,sel);
wav=wav(:,sel);
spe=spe(:,sel);
nor=nor(:,sel);

%% Several quick plots

figure
scatter3(dis(1,:),dis(2,:),dis(3,:),100,squeeze(cor(1,:)),'.')
colorbar

figure
scatter3(dis(1, dis(3,:)==0 ),dis(2, dis(3,:)==0 ),dis(3, dis(3,:)==0 ),100,squeeze(cor(1, dis(3,:)==0 ) ),'.')
view(2)
colorbar

figure
plot(dis(1, dis(2,:)==0 & dis(3,:)==0),squeeze(cor(1, dis(2,:)==0 & dis(3,:)==0)),'.')
colorbar

%% evaluate normalized spatial 2-point corrolation

% evaluate
[dis,cor,wav,spe,nor] = spatial2pointcor(grd,dat1,dat2,'normalized');

% remove empty corrolations
sel=~isinf(nor);
dis=dis(:,sel);
cor=cor(:,sel);
wav=wav(:,sel);
spe=spe(:,sel);
nor=nor(:,sel);

%% Several quick plots

figure
scatter3(dis(1,:),dis(2,:),dis(3,:),100,squeeze(cor(1,:)),'.')
colorbar

figure
scatter3(dis(1, dis(3,:)==0 ),dis(2, dis(3,:)==0 ),dis(3, dis(3,:)==0 ),100,squeeze(cor(1, dis(3,:)==0 ) ),'.')
view(2)
colorbar

figure
plot(dis(1, dis(2,:)==0 & dis(3,:)==0),squeeze(cor(1, dis(2,:)==0 & dis(3,:)==0)),'.')
colorbar

