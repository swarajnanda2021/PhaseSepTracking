function improc_chkboard
%improc_chkboard Semi automated script to process the calibration image in
%   the desired data-format. Semi-automation requires inspection of every
%   calibration image and can be further automated if desired.
%
%   Input:
%       Image data from processing folder
%
%   Output:
%       IMblin.dat baseline from radon transform to fit curve
%       IMcfit.dat curve fitted from baseline to checkerboard edges
%       IMcnod.dat checkerboard nodes from intersection between curves

%% Get global variables
global folder date cal rec chkboard ctrl prop Cfit

%% check for available data create tspan and indexing
count=0; % initiate counting for data cell
IMblin=cell(0); % initiate data cell base lines
IMcfit=cell(0); % initiate data cell curve fit
IMcnod=cell(0); % initiate data cell nodes (intersections)

%% define objective function to find nodes
% substitute path in mapping
vars=sort(symvar(sym(Cfit)));
xh=subs(sym(Cfit),vars,[sym('ah',[1 length(vars)-1]),sym('th')]);
xv=subs(sym(Cfit),vars,[sym('av',[1 length(vars)-1]),sym('tv')]);

% substract lines
ObjA=xh-xv;

% objective function
ObjF=defobjf(ObjA,[]);

% gradient and hessian iObjF
vars=sort(symvar(sym(ObjF)));

% write function
objf=matlabFunction(ObjF,'Vars',vars);%,'Optimize',false);

% numbering knowns
nk=(1:length(vars)-2)';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% grad hess
[gObjF,hObjF]=difobjf(ObjF,vars(nu));
gobjf=matlabFunction(gObjF,'Vars',vars);%,'Optimize',false);
hobjf=matlabFunction(hObjF,'Vars',vars);

% evaluate point on line v
xv=matlabFunction(xv);

%% Get data from image
disp('Get checkerboard image data')

% image grid
[y,x]=meshgrid(1:prop.res(2),1:prop.res(1));

% chkboard grid for radon transform
[YY,XX]=meshgrid(linspace(min(chkboard.outline(2,1:4)),max(chkboard.outline(2,1:4)),25*(ctrl.nnod(2)+1)+1),...
    linspace(min(chkboard.outline(1,1:4)),max(chkboard.outline(1,1:4)),25*(ctrl.nnod(1)+1)+1));%

% presets radon transform @ fixed resolution
cen=[min(XX(:))+max(XX(:)) min(YY(:))+max(YY(:))]/2; % center
ran=[-45 135]; % range in degrees, XX and YY are approximately matched in sampling
thet=(ran(1):ran(2))'; % theta

% make figure
figure

% loop camera
for c=1:prop.res(4)
    
    % loop frames to process
    for n=ctrl.fproc
        % message
        disp(['Process image data camera ',num2str(c),' frame ',num2str(n)])
        
        % image data
        imd=(double(import_frames({folder date cal rec},prop.ext,n,c))+0); % board scatters
        
        % show calibration image
        plt=log(abs(imd)+1);
        surf(x',y',zeros(size(plt')),plt')
        view(2)
        shading flat
        axis tight
        axis equal
        colormap gray
        xlabel('x [px]','interpreter','latex')
        ylabel('y [px]','interpreter','latex')
        title(['Image camera ',num2str(c),' frame ',num2str(n)],'interpreter','latex')
        drawnow
        
        % question dialog whether to process image
        answ=questdlg('Do you wish to process this image?','Process Image','Yes','No','Yes');
        
        % continue
        if strcmp(answ,'Yes')
            
            title('\textbf{ Select 4-polygon with long-side first from checkerboard origin }','interpreter','latex')
            
            drawnow
            
            answ='No'; %initiate
            while strcmp(answ,'No')
                h=impoly;
                
                % get rectangle
                pol=round(getPosition(h))';
                
                % plot
                answ=questdlg('Is selection correct?','Select checkerboard','Yes','No','Yes');
                if strcmp(answ,'No')
                    delete(h)
                end
            end
            
            % bounding box
            bbx=[max([min(pol(1,:)),1]) max([min(pol(2,:)),1])...
                min([max(pol(1,:)),prop.res(1)]) min([max(pol(2,:)),prop.res(2)])];
            
            % crop imd
            crp=imd(bbx(1):bbx(3),bbx(2):bbx(4));
            xc=x(bbx(1):bbx(3),bbx(2):bbx(4));
            yc=y(bbx(1):bbx(3),bbx(2):bbx(4));
            
            % local quad intensity normalization
            %             imn=imminmax(crp,strel('rectangle',ceil(1*(size(crp)-1)./(ctrl.nnod+1))),0.1);
            
            % image coefficient local pol fit
            coef=SGfilt_coef(crp,ctrl.fobj,ctrl.ford,1); % note maximum second ord.
            
            % zero derivative filter image
            d0=SGfilt_eval(coef,[0 0 0]);
            
            % plot selected filter board image
            plt=log(abs(d0)+1);
            surf(xc',yc',plt')
            view(2)
            shading flat
            axis tight
            axis equal
            colormap gray
            xlabel('x [-]','interpreter','latex')
            ylabel('y [-]','interpreter','latex')
            title('Checkerboard image','interpreter','latex')
            drawnow
            pause(0.5)
            
            % compute gradient image
            dx=SGfilt_eval(coef,[1 0 0]);
            dy=SGfilt_eval(coef,[0 1 0]);
            
            % image gradient
            img=sqrt(dx.^2+dy.^2);
            
            % plot
            plt=log(abs(img)+1);
            surf(xc',yc',plt')
            view(2)
            shading flat
            axis tight
            axis equal
            colormap gray
            xlabel('x [-]','interpreter','latex')
            ylabel('y [-]','interpreter','latex')
            title('Image gradient','interpreter','latex')
            drawnow
            pause(0.5)
            
            % image mask
            msk=img>imfilter(img,...
                fspecial('average',...
                ceil(0.5*(size(img)-1)./(ctrl.nnod+1))));%mean(img(:)); % segment by surface
            %             msk=bwmorph(msk,'open',1);%,ceil(norm(size(msk),2)/500));
            %             msk=bwareafilt(msk,1,'largest'); % lines +holes
            msk=bwmorph(msk,'spur'); % remove spurious stuff
            
            plt=double(msk);
            surf(xc',yc',zeros(size(plt')),plt')
            view(2)
            shading flat
            axis tight
            axis equal
            colormap gray
            xlabel('x [-]','interpreter','latex')
            ylabel('y [-]','interpreter','latex')
            title('Binarized gradient','interpreter','latex')
            drawnow
            pause(0.5)
            
            % get intensity normalized edge data checkerboard
            ecb=img.*msk;
            ecb=imminmax(ecb,strel('rectangle',ceil(0.5*(size(crp)-1)./(ctrl.nnod+1))),0.1);
            ecb=log(abs(ecb)+1); % shape: as a gaus distance norm log(Aexp(mux^2)+1)
            
            % plot
            plt=ecb;
            surf(xc',yc',plt')
            view(2)
            shading interp
            axis tight
            axis equal
            colormap gray
            xlabel('x [-]','interpreter','latex')
            ylabel('y [-]','interpreter','latex')
            title('Normalized edge data','interpreter','latex')
            drawnow
            pause(0.5)
            
%                         % downgrade resolution to match XX
%                         bres=max([round(size(ecb)./size(XX)); 1 1],[],1);
%                         lwr=imbox(ecb,bres);
%                         [xb,yb]=boxgrid(xc,yc,[],[],bres);
%             
%                         % plot
%                         plt=lwr;
%                         surf(xb',yb',plt')
%                         view(2)
%                         shading interp
%                         axis tight
%                         axis square
%                         colormap gray
%                         xlabel('x [-]','interpreter','latex')
%                         ylabel('y [-]','interpreter','latex')
%                         title('Normalized edge data','interpreter','latex')
%                         drawnow
%                         pause(0.5)
            
            xb=xc;
            yb=yc;
            lwr=ecb;
            
            % linear image trans - warp polygon to chkboard reference
            H=dirlintrans(inhc2homc(pol),inhc2homc(chkboard.outline(:,1:4)));
            
            % warp coordinates
            [xx,yy]=imtrans(XX,YY,[],H); % deformation has uniform expansion (contant determinant)
            
            % warp image on fixed - lower resolution
            lwr=imdewarp(xb,yb,lwr,xx,yy,'cubic',0);
            
            % plot
            plt=lwr;
            surf(XX',YY',plt')
            view(2)
            shading interp
            axis tight
            axis equal
            colormap gray
            xlabel('x [-]','interpreter','latex')
            ylabel('y [-]','interpreter','latex')
            title('Normalized edge data','interpreter','latex')
            drawnow
            pause(0.5)
            
            answ=questdlg(['Does the edge data mark the ',num2str(sum(ctrl.nnod)),' gridlines on the checkerboard?'],'Process Image','Yes','No','Yes');
            
            if strcmp(answ,'Yes') % continue
                
                % compute radon transform image
                [rad,rho]=radon(flipud(lwr'),thet); % transpose to have theta rho
                rho=rho*abs(XX(2,1)-XX(1,1));
                rad=rad.*(ones(size(rad,1),1)*(sum(rad>0,1)/size(rad,1))); % normalize width
                
                % plot radon transform
                cla
                plt=rad;
                surf(rho',thet',plt')
                shading interp
                view(2)
                axis tight
                axis square
                colormap parula
                xlabel('xp [-]','interpreter','latex')
                ylabel('$\theta [^o]$','interpreter','latex')
                title('Radon transform edge data','interpreter','latex')% est.filt.siz' num2str(sqrt(numel(plt))/50)
                drawnow
                
                % get start estimate peaks in radon transform
                dumf=strel('disk',round(max([4 sqrt(numel(rad))/50])),0);
                [dums,duml]=imrelscl(rad,dumf);
                dumobj=dums-duml;
                [X,~,C,~,W]=featpoints(dumobj,dumf,'maxima'); % peaks
                lv_p=W.*abs(sum(C.*[X(:,1).^2 X(:,1).*X(:,2) X(:,2).^2 X(:,1) X(:,2) ones(size(X(:,1)))],2)); % level
                rho_p=interp1(rho,X(:,1),'linear','extrap'); % rho-axis
                thet_p=interp1(thet,X(:,2),'linear','extrap'); % theta-axis
                
                % get peaks long edge
                sel=(thet_p>-10 & thet_p<10 );
                thet_p1=thet_p(sel);
                rho_p1=rho_p(sel);
                lv_p1=lv_p(sel);
                [lv_p1,i]=sort(lv_p1,'descend');
                rho_p1=rho_p1(i);
                thet_p1=thet_p1(i);
                
                % get peaks short edge
                sel=(thet_p>80 & thet_p<100 );
                thet_p2=thet_p(sel);
                rho_p2=rho_p(sel);
                lv_p2=lv_p(sel);
                [lv_p2,i]=sort(lv_p2,'descend');
                rho_p2=rho_p2(i);
                thet_p2=thet_p2(i);
                
                % make sure enough peaks present by appending zeros
                rho_p1=  cat(1,rho_p1  ,   zeros(ctrl.nnod(1),1));
                rho_p2=  cat(1,rho_p2  ,   zeros(ctrl.nnod(2),1));
                thet_p1= cat(1,thet_p1 ,   zeros(ctrl.nnod(1),1));
                thet_p2= cat(1,thet_p2 , 90*ones(ctrl.nnod(2),1));
                lv_p1= cat(1,lv_p1 , max(rad(:))*ones(ctrl.nnod(1),1));
                lv_p2= cat(1,lv_p2 , max(rad(:))*ones(ctrl.nnod(2),1));
                
                % select peaks
                rho_p=[rho_p1(1:ctrl.nnod(1))
                    rho_p2(1:ctrl.nnod(2))];
                thet_p=[thet_p1(1:ctrl.nnod(1))
                    thet_p2(1:ctrl.nnod(2))];
                lv_p=[lv_p1(1:ctrl.nnod(1))
                    lv_p2(1:ctrl.nnod(2))];
                
                % plot
                hold on
                plot3(rho_p,thet_p,max(rad(:))+lv_p,'rsq')%
                hold off
                
                % select right one in sort list
                answ='No';
                while strcmp(answ,'No')
                    
                    % question
                    answ=questdlg(['Are ',num2str(sum(ctrl.nnod)),' peaks marked correct? Else relocate spurious by clicking'],'Check peaks','Yes','No','Yes');
                    
                    % plot peaks
                    if strcmp(answ,'No')
                        pause
                        
                        title('\textbf{Click incorrect peak, press spacebar and relocate by clicking}','interpreter','latex')
                        
                        sel=getCursorInfo(datacursormode(gcf));
                        
                        if isfield(sel,'DataIndex')
                            ind=sel.DataIndex;
                            if ind<=sum(ctrl.nnod)
                                delete(gco)
                                
                                pause
                                
                                sel=getCursorInfo(datacursormode(gcf));
                                
                                if isfield(sel,'Position')
                                    rho_p(ind)=sel.Position(1);
                                    thet_p(ind)=sel.Position(2);
                                    lv_p(ind)=sel.Position(3);
                                    
                                    % re-plot
                                    hold on
                                    plot3(rho_p,thet_p,max(rad(:))+lv_p,'rsq')
                                    hold off
                                end
                                
                            end
                        end
                    end
                    
                end
                
                % sort peaks
                [thet_p,i]=sort(thet_p,'descend'); % sort thet
                rho_p=rho_p(i);
                lv_p=lv_p(i);
                k=ctrl.nnod(2);%find(min(diff(thet_p))==diff(thet_p)); % indicator, can be replace by ctrl.nnod..
                
                [rho_p(1:k),i]=sort(rho_p(1:k),'ascend');
                thet_p(1:k)=thet_p(i);
                lv_p(1:k)=lv_p(i);
                
                [rho_p(k+1:end),i]=sort(rho_p(k+1:end),'ascend');
                thet_p(k+1:end)=thet_p(k+i);
                lv_p(k+1:end)=lv_p(k+i);
                
                % plot sorting at radon
                hold on
                for i=1:length(thet_p)
                    text(rho_p(i),thet_p(i),max(rad(:))+lv_p(i),num2str(i),'Color',[1 0 0])
                end
                hold off
                pause(0.5)
                
                % write to hom lines, first nodes along x (vertical lines) then y (horizontal lines)
                ly=[cosd(thet_p(k+1:end)) sind(thet_p(k+1:end)) -(rho_p(k+1:end)+cen(1)*cosd(thet_p(k+1:end))+cen(2)*sind(thet_p(k+1:end)))]';
                lx=[cosd(thet_p(1:k)) sind(thet_p(1:k)) -(rho_p(1:k)+cen(1)*cosd(thet_p(1:k))+cen(2)*sind(thet_p(1:k)))]';
                lcb=[lx ly]; % lines checkerboard
                
                % get outer perimeter board
                pcb=zeros(4,size(lcb,2));
                mval=(ctrl.nnod(1)+2)/ctrl.nnod(1);
                for i=1:size(lx,2)
                    xmin=homc2inhc(cross(lcb(:,i),ly(:,1))); % xmin
                    xmax=homc2inhc(cross(lcb(:,i),ly(:,end))); %xmax
                    xcen=(xmin+xmax)/2;
                    xdel=xmax-xmin;
                    pcb(1:2,i)=xcen-mval*xdel/2;
                    pcb(3:4,i)=xcen+mval*xdel/2;
                end
                mval=(ctrl.nnod(2)+2)/ctrl.nnod(2);
                for i=size(lx,2)+1:size(lx,2)+size(ly,2)
                    xmin=homc2inhc(cross(lcb(:,i),lx(:,1))); % xmin
                    xmax=homc2inhc(cross(lcb(:,i),lx(:,end))); %xmax
                    xcen=(xmin+xmax)/2;
                    xdel=xmax-xmin;
                    pcb(1:2,i)=xcen-mval*xdel/2;
                    pcb(3:4,i)=xcen+mval*xdel/2;
                end
                
                % transform back to camera image
                pcb=[homc2inhc(H\inhc2homc(pcb(1:2,:)))
                    homc2inhc(H\inhc2homc(pcb(3:4,:)))];
                lcb=H'*lcb; % lines transform along hom. conj.
                
                % get image coordinates and intensity values
                cbdat=zeros(5,nnz(msk));
                cbdat(1:3,:)=[xc(msk) yc(msk) ecb(msk)]';
                
                %figure; scatter(cbdat(1,:),cbdat(2,:),[],cbdat(3,:),'.')
                
                % divide data associated with each line
                [ D ] = abs(pntlindist(cbdat(1:2,:) , lcb )); % [px]
                [D,i]=min(D,[],2); % min index
                cbdat(4:5,:)=[D  i]'; % distance pixel location and indexing
                
                %figure; plot(D,'.')
                %figure; scatter(cbdat(1,:),cbdat(2,:),[],D,'.')
                
                % plot edge data checkerboard
                cla
                plt=ecb;
                surf(xc',yc',zeros(size(plt')),plt')
                view(2)
                shading interp
                axis tight
                axis square
                colormap gray
                xlabel('x [px]','interpreter','latex')
                ylabel('y [px]','interpreter','latex')
                title('Checkerboard Image','interpreter','latex')
                drawnow
                pause(0.5)
                
                % plot houghlines
                hold on
                for i=1:size(pcb,2)
                    plot([pcb(1,i) pcb(3,i)],[pcb(2,i) pcb(4,i)],'r','LineWidth',1)
                    text(pcb(1,i),pcb(2,i),num2str(i),'Color','r')
                end
                hold off
                drawnow
                
                % curve fit checkerboard and plot result
                t=linspace(0,1,100); % min(pixels linewidth,dev)
                dev=max([sqrt(range(cbdat(1,:))^2+range(cbdat(2,:))^2)...
                    /sqrt(prod(ctrl.nnod))*ctrl.pdis,2]);
                cfit=zeros(2*(ctrl.cord+1),sum(ctrl.nnod));
                ox=0:ctrl.cord; % order coef
                hold on
                for i=1:sum(ctrl.nnod)%number lines
                    % get line data, within confidence interval
                    ldat=cbdat(1:4,cbdat(5,:)==i); % cut of the ends by this value
                    
                    % rotate to base coordinates
                    nl=[lcb(1,i)
                        lcb(2,i)]/sqrt(lcb(1,i)^2+lcb(2,i)^2);
                    nx=[1
                        0];
                    ny=[0
                        1];
                    ang=-sign(dot(nl,nx))*acosd(dot(nl,ny));
                    R=rotz(ang);
                    R=R(1:2,1:2);
                    ldat(1:2,:)=R\ldat(1:2,:); % w.r.t. base line frame
                    
                    % crop domain line
                    xdom=sort([dot(nx,R\pcb(1:2,i)) dot(nx,R\pcb(3:4,i))]);
                    ldat(:,(ldat(1,:)<xdom(1)) | (ldat(1,:)>xdom(2)))=[];
                    %                 yran=range(xdom)*0.05;
                    ldat(:,ldat(4,:)>dev)=[];
                    %cla; plot(ldat(1,:),ldat(2,:),'.'); drawnow; pause(0.5)
                    
                    % define base
                    xdat=ldat(1,:)'; % x coordinates
                    base=[min(xdom) range(xdom)]; % define base with t in 0 1
                    
                    % fit/regress n-nd order curve w.r.t. base coordinates
                    ydat=ldat(2,:)'; % y coordinates
                    wdat=sqrt(abs(   ldat(3,:) )'); % weight for intensity
                    curv=regress(wdat.*ydat,wdat.*(xdat.^ox)); % curve
                    cdist=abs(ydat-(xdat.^ox)*curv);
                    %                 cdist=ldat(4,:)';
                    
                    % iterate 10 times by residual, not necc. conv..
                    for j=1:10
                        coutl=cdist>sqrt(mean(cdist.^2));
                        curv=regress(wdat(~coutl).*ydat(~coutl),wdat(~coutl).*(xdat(~coutl).^ox)); % curve
                        cdist=abs(ydat-(xdat.^ox)*curv);
                        %                     figure; plot(xdat,cdist,'.');hold on; plot(xdat,coutl,'.'); hold off
                    end
                    curv=curv';
                    
                    % define curve fit
                    fdat=zeros(2,ctrl.cord+1);
                    
                    % substitute base curve
                    fdat(1,1:2)=base;
                    
                    % substitute by binomial expansion
                    for j=0:ctrl.cord
                        for k=0:j
                            fdat(2,k+1)=fdat(2,k+1)+curv(j+1)*nchoosek(j,k)*base(1)^(j-k)*base(2)^k;
                        end
                    end
                    %cy=[cy(1)+cy(2)*cx(1)+cy(3)*cx(1)^2 cy(2)*cx(2)+cy(3)*cx(1)*cx(2)*2 cy(3)*cx(2)^2];
                    
                    % rotate (back) to image frame
                    fdat=R(1:2,1:2)*fdat; %i
                    
                    % plot
                    X=fdat*(t'.^ox)';
                    plot(X(1,:),X(2,:),'g','LineWidth',1)
                    
                    % save
                    cfit(:,i)=reshape(fdat',[],1);
                end
                hold off
                pause(0.5)
                
                % find intersections
                xcb=cell(ctrl.nnod(1),ctrl.nnod(2));
                ax=cfit(:,1:ctrl.nnod(2)); % lines along x
                ay=cfit(:,ctrl.nnod(2)+1:ctrl.nnod(2)+ctrl.nnod(1)); % lines along y
                for i=1:size(ax,2) % select line along y
                    for j=1:size(ay,2) % loop over lines along x
                        % coef data
                        adat=[ax(:,i);ay(:,j)];
                        
                        % numering assembly
                        asn=repmat((1:length(nu))',1,size(adat,2)); % assembly numbering
                        
                        % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
                        dat=[zeros(size(nk)) nk adat
                            ones(size(nu)) nu asn];
                        
                        % initiate solution vector
                        t0=[i/(ctrl.nnod(2)+2)
                            j/(ctrl.nnod(1)+2)]; % or do 1 ord coef a
                        
                        % solve global optimization (newton is here exact solution, but do generalize)
                        [tsol,~] = optobjf( [] , objf , gobjf , hobjf , t0 , dat , ctrl.optl ,'newton','no');
                        
                        inp=num2cell([ax(:,i);tsol(1)]',1);
                        
                        xcb{j,i}=xv(inp{:}); % eval vert line
                    end
                end
                xcb=cat(2,xcb{:,:}); % some numbering
                
                % plot checkerboard
                cla
                plt=log(abs(d0)+1);
                surf(xc'-1/2,yc'-1/2,zeros(size(plt')),plt') % shift correct centered color plot
                view(2)
                shading flat
                axis tight
                axis square
                colormap gray
                xlabel('x [px]','interpreter','latex')
                ylabel('y [px]','interpreter','latex')
                title('Identified checkerboard','interpreter','latex')
                drawnow
                
                % plot houghlines
                hold on
                for i=1:size(pcb,2)
                    plot([pcb(1,i) pcb(3,i)],[pcb(2,i) pcb(4,i)],'r','LineWidth',1)
                    text(pcb(1,i),pcb(2,i),num2str(i),'Color','r')
                end
                hold off
                drawnow
                
                % curve fit plot result
                hold on
                for i=1:sum(ctrl.nnod)%number lines
                    % plot
                    X=reshape(cfit(:,i),[],2)'*(t'.^ox)';
                    plot(X(1,:),X(2,:),'g','Linewidth',1)
                end
                hold off
                
                % plot coordinates with numbering on deformed image
                hold on
                for i=1:size(xcb,2)
                    plot(xcb(1,i),xcb(2,i),'msq')
                    text(xcb(1,i),xcb(2,i),num2str(i),'Color','magenta');
                end
                hold off
                pause(0.5)
                
                % question
                answ=questdlg('Is the identified checkerboard correct?','Check peaks','Yes','No','Yes');
                
                if strcmp(answ,'Yes')
                    count=count+1;
                    
                    % n frame, c view, l base
                    IMblin{count}=[n*ones(1,sum(ctrl.nnod))
                        c*ones(1,sum(ctrl.nnod))
                        (1:sum(ctrl.nnod))
                        lcb];
                    
                    % n frame, c view, coef x, coef y
                    IMcfit{count}=[n*ones(1,sum(ctrl.nnod))
                        c*ones(1,sum(ctrl.nnod))
                        (1:sum(ctrl.nnod))
                        cfit];
                    
                    % write the nodes
                    IMcnod{count}=[n*ones(1,prod(ctrl.nnod))
                        c*ones(1,prod(ctrl.nnod))
                        (1:prod(ctrl.nnod))
                        xcb];
                    
                    % display message
                    disp('Identified checkerboard correctly')
                    
                end % if correct / satisfied
                
            end % if image gradient
            
        end % if processing
        
    end % n
    
end % c

%% save
IMblin=cell2mat(IMblin);
IMcfit=cell2mat(IMcfit);
IMcnod=cell2mat(IMcnod);
if ~isempty(IMblin) && ~isempty(IMcfit) && ~isempty(IMcnod)
    fsave([folder date cal rec '_IMblin.dat'],IMblin,'w'); % append?
    fsave([folder date cal rec '_IMcfit.dat'],IMcfit,'w'); % append?
    fsave([folder date cal rec '_IMcnod.dat'],IMcnod,'w'); % write
end

end