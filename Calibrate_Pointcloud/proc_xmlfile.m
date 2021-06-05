function proc_xmlfile
%proc_xmlfile Process and check the calibration data from a Davis xml file
%   for consitency with the Davis image file, in extension to preference.
%   This can be extended to any similar type calibration format.
%
%   $Author: Koen Muller$
%
%   Notes,
%       - Input: .xml file from Davis
%       - Output: .dat&bin file containing binary data format and size
%       - This can be extended to any similar type calibration format.
%

%% Get global variables
global folder date cal prop ext ctrl

%% Load or create data
disp('Load calibration marks from Davis .xml')
tic

% read xml file using custom function from file exchange
obj=xmlread([folder date cal vsl 'Calibration' vsl 'MarkPositionTable.xml']);
obj=xml2struct(obj);

% initiate
count=0; % count variable
Cmark=cell(0); % data cell
for c=1:length(obj.MarkTable.Camera) % loop cameras
    
    % figure per camera
    figure(c) % figure for plotting
    
    % loop targets
    for v=1:length(obj.MarkTable.Camera{c}.View) % loop planes/view calibration taget
        
        % image data
        imd=log(double(import_frames({folder date [cal vsl 'Calibration' vsl 'camera' num2str(c)]},prop.ext,v,1))+1); % board scatters
        
        % local grid
        [y,x]=meshgrid(1:size(imd,2),1:size(imd,1)); % meshgrid image
        
        % subplot image
        subplot(1,2,1)
        cla
        surf(x'-1/2,y'-1/2,zeros(size(imd')),imd') % shift to mid pixel 
        view(2)
        shading flat
        axis tight
        axis equal
        colormap gray
        xlabel('x [px]','interpreter','latex')
        ylabel('y [px]','interpreter','latex')
        title(['Calibration points check: camera = ',num2str(c),' view = ',num2str(v)],'interpreter','latex')
        
% 		% detect points
%         [scls,scll]=imrelscl(imd,ctrl.fobj);
%         obj=(scls - scll); % use abs to include minima 
%         [xpnts]=featpoints(obj,ctrl.fobj,'minima');
% %         I=I(:,1)<mean(I(:,1))-std(I(:,1));
%         ypnts=xpnts(:,2);
%         xpnts=xpnts(:,1); % I
%         hold on;
%         plot(xpnts,ypnts,'.')
%         hold off;
        
        % loop over calibration points
        for m=1:length(obj.MarkTable.Camera{c}.View{v}.Mark)
            % camera points
            xc=str2double(obj.MarkTable.Camera{c}.View{v}.Mark{m}.RawPos.Attributes.x);
            yc=str2double(obj.MarkTable.Camera{c}.View{v}.Mark{m}.RawPos.Attributes.y);
            
            % world points
            X=str2double(obj.MarkTable.Camera{c}.View{v}.Mark{m}.FittedWorldPos.Attributes.x);
            Y=str2double(obj.MarkTable.Camera{c}.View{v}.Mark{m}.FittedWorldPos.Attributes.y);
            Z=str2double(obj.MarkTable.Camera{c}.View{v}.Mark{m}.FittedWorldPos.Attributes.z);
            
            % correct world points: modify here,
            xc=1+xc; % to flip axis prop.res(1)
            %yc=1+(prop.res(2)-1-yc); % to flip axis
            yc=1+(size(imd,2)-1-yc);%(prop.res(2)-1-yc); % to flip axis
            
%             % nearest neighbor in featpoint list (ref)
%             [~,i] = min(sum([xpnts-xc ypnts-yc].^2,2)); 
%             xc=xpnts(i);
%             yc=ypnts(i); % consistency image processing

            % manual hacks
          %  if c==5
          %      yc=yc-348; % correct manual
          %  end
            
            % write
            count=count+1;
            Cmark{count}=[c xc yc X Y Z]';
            
            % append plot
            subplot(1,2,1)
            hold on
            if round(X,1)==0 && round(Y,1)==0
                plot(xc,yc,'g+')
            elseif round(X,1)==0.1 && round(Y,1)==0
                plot(xc,yc,'b*')
            elseif round(X,1)==0 && round(Y,1)==0.1
                plot(xc,yc,'r*')
            else
                plot(xc,yc,'m.')
            end
            hold off
            
            % marker object positions            
            subplot(1,2,2)
            hold on
            if round(X,1)==0 && round(Y,1)==0
                plot3(X,Y,Z,'g+')
            elseif round(X,1)==0.1 && round(Y,1)==0
                plot3(X,Y,Z,'b*')
            elseif round(X,1)==0 && round(Y,1)==0.1
                plot3(X,Y,Z,'r*')
            else
                plot3(X,Y,Z,'m.')
            end
            hold off
            
        end
        
        % subplot target
        subplot(1,2,2)
        view(3)
        axis tight
        axis equal
        xlabel('X [m]','interpreter','latex')
        ylabel('Y [m]','interpreter','latex')
        zlabel('Z [m]','interpreter','latex')
        title(['Calibration Target Position Camera ',num2str(c)],'interpreter','latex')
        
        % draw figure
        drawnow
        
    end
    
end

%% save calibration data in binary format
Cmark=cell2mat(Cmark); % needed for next proc
fsave([folder date cal vsl 'Cmark.dat'],Cmark,'w')

end