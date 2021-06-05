function plot_calpoints(Dpar,Kmat,Rmat,tvec)
%plot_calpoints Plot calibration
%   Plots the resulting camera calibration in a visualization of
%       - world point of the mark data
%       - (virtual-)camera centers
%       - triangulated calibration point from multiview setup (todo)
%

%% Get global variables
global folder date cal prop ctrl Pmap Dmap

%% Load or create data (takes some time)
Cmark=fload([folder date cal vsl 'Cmark.dat']);

%% projection matrix
Pmat=krtm2pmat(Kmat,Rmat,tvec);

%% define objective
%world coordinates
X=sym('X',[3 1]);

% image coordinates
x=sym('x',[2 1]);

% projected points
xp=subs(sym(Pmap),sym('X',[4 1]),[X;1]);

% objective argument
ObjA=x-xp;

% define objective function
ObjF=defobjf(ObjA,[]);

% variabes
vars=sort(symvar(sym(ObjF)));

objf=matlabFunction(ObjF,'Optimize',false,'Vars',vars);

%% plot cams
figure
hold on
% global frame at origin
quiver3(0,0,0,1,0,0,10,'k.') % scale here
quiver3(0,0,0,0,1,0,10,'k.')
quiver3(0,0,0,0,0,1,10,'k.')
% camera centers
for c=1:length(tvec)
    xcen=-Rmat{c}\tvec{c};
    Xor=Rmat{c}\eye(3);
    
    plot3(xcen(1),xcen(2),xcen(3),'k.')
    quiver3(xcen(1),xcen(2),xcen(3),Xor(1,1),Xor(2,1),Xor(3,1),10,'r.') % scale here
    quiver3(xcen(1),xcen(2),xcen(3),Xor(1,2),Xor(2,2),Xor(3,2),10,'g.')
    quiver3(xcen(1),xcen(2),xcen(3),Xor(1,3),Xor(2,3),Xor(3,3),10,'b.')
    text(xcen(1),xcen(2),xcen(3),num2str(c))
end
hold off
axis square

%% plot imposed world points
X=unique(Cmark(4:6,:)','rows')'; 
hold on
plot3(X(1,:),X(2,:),X(3,:),'.','Color',[0.75 0.75 0.75])
hold off
view(3)
pause(1)

% xlim([-2 2])
% ylim([-2 2])
% zlim([-1 3])

%% plot triangulated world coordinates

if 0
    
hold on
% numbering knowns
nk=[1:12 16,17]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% find grad and hess w.r.t. unknowns
[gObjF,hObjF]=difobjf(ObjF,vars(nu));

gobjf=matlabFunction(gObjF,'Optimize',false,'Vars',vars);

hobjf=matlabFunction(hObjF,'Optimize',false,'Vars',vars);

% loop 2 plot
for i=unique(Cmark(4,:))
    I=Cmark(4,:)==i;
    for j=unique(Cmark(5,:))
        J=Cmark(5,:)==j;
        for k=unique(Cmark(6,:))
            K=Cmark(6,:)==k;
            
            % get cameras
            c=Cmark(1,I&J&K);
            x=Cmark(2:3,I&J&K);
            
            if length(c)>=2
                % distortion dewarp
                b=cell2mat(cellfun(@(x)x',Dpar(c),'UniformOutput',false));
                
                % distortion dewarp points
                inp=num2cell([b;x],2);
                xd=Dmap(inp{:});
                
                % linear triangulation points
                xh=num2cell(inhc2homc(xd),1);
                P=Pmat(c);
                X=homc2inhc(lintriang(xh,P));
                
                % calibration parameters
                p=cell2mat(cellfun(@(x)reshape(x',[],1),P,'UniformOutput',false));
                
                % solve objective function
                asn=repmat((1:length(nu))',1,size(xd,2)); % assembly numbering
                
                % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
                dat=[zeros(size(nk)) nk [p;xd]
                    ones(size(nu)) nu asn];
                
                % solution vec
                x0=X;
                
                % solution
                [xsol,res] = optobjf( [], objf, gobjf, hobjf , x0 , dat , 10^(-3) ,'newton');
                
                % plot scatter with color residual
                scatter3(xsol(1),xsol(2),xsol(3),[],res,'.')
                
                drawnow
            end
        end
    end
end

end
hold off

%% figure lay-out
axis tight
axis equal
grid minor
box on
xlabel('x $[\rm{mm}]$','interpreter','latex')
ylabel('y $[\rm{mm}]$','interpreter','latex')
zlabel('z $[\rm{mm}]$','interpreter','latex')
view(3)
title('Reconstruction of the (virtual-)camera centers and calibration points','interpreter','latex')
hold off
drawnow

end