function plot_imageprocessing
%plot_improc Plot image processing to check validity checkerboard and
%   exclude wrong ones

%% get globals
global folder date cal prop ctrl Cfit 

%% files
fil=dir([folder date cal vsl '*_IMcnod.dat']); % keep

%% create path for figure
if exist([folder date cal vsl 'Figures' vsl 'ImageProcessing'],'dir')~=7
    mkdir([folder date cal vsl 'Figures' vsl 'ImageProcessing'])
end

%% loop files

% curve sampling
t=linspace(0,1,100);

% image grid
[x,y]=ndgrid(1:prop.res(1),1:prop.res(2));

% down sample
[x,y]=boxgrid(x,y,[],[],[3 3]);

% initiate figure
h=figure(1);
h.Color=[1 1 1];% window
h.Position=[50 50 500 500];

% loop
for k=1:length(fil)
    
    %name
    name=fil(k).name(1:end-11);
    
    % curve data
    IMcfit=fload([folder date cal vsl name '_IMcfit.dat']);
    
    %data    
    IMcnod=fload([folder date cal vsl name '_IMcnod.dat']);
    
    for n=unique(IMcnod(1,:))
        
        % frame set
        N1=IMcnod(1,:)==n;
        N2=IMcfit(1,:)==n;
        
        for c=unique(IMcnod(2,N1))
            
            % clear image
            cla
            
            % camera selection
            C1=IMcnod(2,:)==c;
            C2=IMcfit(2,:)==c;
            
            % nodes
            X=IMcnod(4:5,C1&N1);
            
            % curves
            acb=IMcfit(4:end,C2&N2);
            
            % image
            imd=double(import_frames({folder date cal [vsl name]},prop.ext,n,c)); % board scatters
            
            % downgrade
            imd=imbox(imd,[3 3],'avg');
            
            % plot image
            plt=log(imd+1);
            surf(x'-1/2,y'-1/2,zeros(size(plt')),plt')
            caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+3*std(plt(:))])
            

            % plot curves
            hold on
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
            
            % plot nodes
            plot(X(1,:),X(2,:),'rsq')
            
            % define title
            if strcmp(name(23:39),'Calibration_Route')
                    figname=[name(12:42) '_cam_',num2str(c),'_frm_',num2str(n)];
            elseif strcmp(name(23:43),'Distortion_Correction') 
                    figname=[name(12:46) '_cam_',num2str(c),'_frm_',num2str(n)];
            end
            
            % layout
            view(2)
            shading flat
            axis equal
            colormap gray
            xlim([min(X(1,:))-range(X(1,:))/2 max(X(1,:))+range(X(1,:))/2])
            ylim([min(X(2,:))-range(X(2,:))/2 max(X(2,:))+range(X(2,:))/2])
            title(figname,'Interpreter','none')
            drawnow
            
            % print figure
            print(h,[folder date cal vsl 'Figures' vsl 'ImageProcessing' vsl name '_cam_',num2str(c),'_frm_',num2str(n)],'-opengl','-r100','-dbmp')
            
        end
        
    end
    
end


end

