function [Tlink,frm_new] = ind_Tlink(Cpln,frm_rem,frm_new)
%ind_Tlink Index temporal links to camera data
%
%   Input,
%       Cpln Camera plane data with indexing
%
%   Output,
%       Tlink Indexed temporal links between ellipse identification

%% Get global variables
global folder date rec prop ctrl plotting

%% check for available data create tspan and indexing

if exist([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Tlink.dat'],'file')~=0 % load
    
    % load previous timelink data
    Tlink=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Tlink.dat']);
    
    % maximum index
    tmax=max(Tlink(1,:));
    
    % remove processed frame if already passed this step
    frm_new=frm_new(:,~ismember(frm_new,Tlink(5,:)));
%     Cpln(3,...
%         ismember(Cpln(1,:),Tlink([4 6],:))) ));
    
else % initiate
    
    % start time link data
    Tlink=zeros(7,0); % [track camera frame(n-1) i-camera frame(n) j-camera avg-disp-ellipse]
    
    % temporal link index
    tmax=0;
    
end

%% Get data from image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop cameras to index data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop cameras
for c=unique(Cpln(2,:)) %loop different cameras
    
    % loop frames
    for n=frm_new
        
        %%%%%%%%%%%%%%%%%%%%%
        %%% Begin message %%%
        %%%%%%%%%%%%%%%%%%%%%
        
        disp(['Temporal linking camera ',num2str(c),' frame pair (',num2str(n-1),',',num2str(n),')'])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Load data in frame pair %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % open data 
        cpln=Cpln(:, Cpln(2,:)==c & ( Cpln(3,:)==(n-1) | Cpln(3,:)==n ) );
        
        % previous time link indexing
        tlink=Tlink(:,ismember(Tlink(6,:), cpln(1, cpln(3,:)==n-1 ) )); % In in Inm1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% define tracking adjacency %%% % within stats @ memory
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % define data previous frame
        Cnm1=cpln(4:9, cpln(3,:)==n-1 );
        Xnm1=cvec2eshp(Cnm1);
        Vnm1=cpln(10:11, cpln(3,:)==n-1 );
        
        %Xe=cvec2pnts(-Cnm1); % avoid imagionary
        %figure; plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],[squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-g','LineWidth',1)
        
        % shift n forward to n+1/2
        Xf=Xnm1+Vnm1/2;
        Cf=contrans(Cnm1,Vnm1/2,'displace'); % shift and then ..
        
        %Xe=cvec2pnts(-Cf); % avoid imagionary
        %hold on; plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],[squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-m','LineWidth',1); hold off;
        
        % expand
        Cf=contrans(Cf,sqrt(sum(Vnm1.^2,1))/2,'expand'); % expand to capture current and next positions
        
        %Xe=cvec2pnts(-Cf); % avoid imagionary
        %hold on; plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],[squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-r','LineWidth',1); hold off;
        
        % define data current frame
        Cn=cpln(4:9, cpln(3,:)==n );
        Xn=cvec2eshp(Cn);
        Vn=cpln(10:11, cpln(3,:)==n );
        
        %Xe=cvec2pnts(-Cn); % avoid imagionary
        %figure; plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],[squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-g','LineWidth',1)
        
        % shift n+1 backward to n+1/2
        Xb=Xn-Vn/2;
        Cb=contrans(Cn,-Vn/2,'displace');
        
        %Xe=cvec2pnts(-Cb); % avoid imagionary
        %hold on; plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],[squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-m','LineWidth',1); hold off;
        
        % expand
        Cb=contrans(Cb,sqrt(sum(Vn.^2,1))/2,'expand');
        
        %Xe=cvec2pnts(-Cb); % avoid imagionary
        %hold on; plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],[squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-r','LineWidth',1); hold off;
        
        % disparity in displacement at n-1/2
        Ef=ellipsedist(Cf,Xb,'matrix');
        Eb=ellipsedist(Cb,Xf,'matrix');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% solve the tracking adjacency problem %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % adjacency matrix - with front back consistency @ multiple scales?
        T_adj=sparse( Ef <= ctrl.dspl & Eb' <= ctrl.dspl ) ; % | is the complementary definition of being sure its uniquely tracked.
        
        %figure; spy(T_adj)
        
        % camera indexing
        Inm1=cpln(1, cpln(3,:)==n-1 );
        In=cpln(1, cpln(3,:)==n );
        
        % new time link indexing
        [i,j]=find( T_adj ); 
        tnew=[zeros(1,length(i)) % initiate track index
            c*ones(1,length(i)) % camera
            (n-1)*ones(1,length(i)) % frame n-1
            reshape(Inm1(i),1,[]) % index image n-1
            n*ones(1,length(j)) % frame n
            reshape(In(j),1,[]) % index image n
            reshape(Ef(sub2ind(size(Ef),i,j))+Eb(sub2ind(size(Eb),j,i)),1,[])/2]; % displacement between frames
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% link time links to previous links %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % link previous time link indexing to new links
        uc=unique(Inm1);
        
        % track index adjacency
        [utp,~,itp]=unique(tlink(6,:)); % In
        iuc=find(ismember(uc,utp));
        itp=iuc(itp);
        T_prev=sparse(itp,1:length(itp),ones(size(itp)),length(uc),length(itp));
        
        %figure; spy(T_prev)
        
        % track index adjacency
        [utn,~,itn]=unique(tnew(4,:)); % Inm1
        iuc=find(ismember(uc,utn));
        itn=iuc(itn);
        T_new=sparse(itn,1:length(itn),ones(size(itn)),length(uc),length(itn));
        
        %figure; spy(T_new)
        
        % temporal link adjacency
        T_lnk=T_new'*T_prev;
        
        %figure; spy(T_lnk)
        
        % unique time correspondence
        Unm1= sum(T_lnk,2)==1 ; % unique from nm1 to n
        Un= sum(T_lnk,1)==1 ; % unique from n to nm1
        
        % map index from uniquely adjacent ellipses
        [i,j]=find( Unm1 & T_lnk & Un ); % find
        tnew(1,i)=tlink(1,j); % write by mapping
        
        %figure; spy( Unm1 & T_lnk & Un );
        
        % index new tracks that are non adjacent
        t_new=(1:nnz(tnew(1,:)==0));
        tnew(1,tnew(1,:)==0)=tmax+t_new; % write new
        tmax=tmax+length(t_new);% raise imax level
        
        %figure; spy( T_lnk - ( Unm1 & T_lnk & Un ) );
        
        %%%%%%%%%%%%%%%%%%%%
        %%% plot results %%%
        %%%%%%%%%%%%%%%%%%%%
        
        if strcmp(plotting,'on')
            
            % figure
            figure(c)
            
            % clear axis
            cla 
            hold on
            
            % coordinates all
            X=cvec2eshp(cpln(4:9,:));
            
            % ellipses nm1
            Xe=cvec2pnts(-Cnm1); % avoid imagionary
            plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-g','LineWidth',1)
            
            % position
            plot(Xnm1(1,:),Xnm1(2,:),'.g','MarkerSize',5)
            
            % velocity
            quiver(Xnm1(1,:),Xnm1(2,:),Vnm1(1,:),Vnm1(2,:),0,'g','LineWidth',1)
            
            % ellipses n
            Xe=cvec2pnts(-Cn); % avoid imagionary
            plot([squeeze(Xe(1,:,:));squeeze(Xe(1,1,:))'],...
                [squeeze(Xe(2,:,:));squeeze(Xe(2,1,:))'],'-r','LineWidth',1)
            
            % position
            plot(Xn(1,:),Xn(2,:),'.r','MarkerSize',5)
            
            % velocity
            quiver(Xn(1,:),Xn(2,:),Vn(1,:),Vn(2,:),0,'r','LineWidth',1)
            
            % unique identifications
            [uc,~,ic]=unique(cpln(1,:));
            C_pln=sparse(ic,1:length(ic),ones(1,length(ic)),length(uc),length(ic));
            
            %figure; spy(C_pln)
            
            % unique tracks
            [utp,~,itp]=unique(tnew(4,:));
            iuc=find(ismember(uc,utp));
            itp=iuc(itp);
            [utn,~,itn]=unique(tnew(6,:));
            iuc=find(ismember(uc,utn));
            itn=iuc(itn);
            T_adj=sparse(itp,itn,ones(size(itn)),length(uc),length(uc));
            
            %figure; spy(T_adj)
            
            % adjaceny
            [i,j]=find(C_pln'*T_adj*C_pln);
            
            % plot
            quiver(X(1,i),X(2,i),X(1,j)-X(1,i),X(2,j)-X(2,i),0,'b.','LineWidth',1)
            
            % axis
            hold off
            xlabel('x-ref [m]','interpreter','latex')
            ylabel('y-ref [m]','interpreter','latex')
            title(['tracking cam ',num2str(c),' frm ',num2str(n)],'interpreter','latex')
            axis tight
            axis equal
%             xlim([-0.125 -0.075]);ylim([-0.025 0.025])
%             xlim([-0.1 -0.0]);ylim([-0.05 0.05])
            
            drawnow
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        %%% write results %%%
        %%%%%%%%%%%%%%%%%%%%%
        
        % write to unstructured data set with rowspace known
        Tlink=cat(2,Tlink,tnew);
        
        %%%%%%%%%%%%%%%%%%%
        %%% End message %%%
        %%%%%%%%%%%%%%%%%%%
        
        % display message
        disp(['Indexed ',num2str(nnz(~ismember(tnew(1,:),tlink(1,:)))),' new out of ',num2str(numel(tnew(1,:))),' tracks'])
        
    end % n
    
end % c

%% Remove data

% trash
Tlink=Tlink(:,~ismember(Tlink(5,:),frm_rem));

%% Plot camera stats
if strcmp(plotting,'on')
    
    for c=1:prop.res(4)
        
        % make figure
        figure(c)
        
        % select camera
        C=Tlink(2,:)==c;
        
        % average number of tracks per frame
        subplot(1,2,1)
        N= histcounts(Tlink(5,C),[unique(Tlink(5,C))-1/2,max(Tlink(5,C))+1/2]) ; % reference frame n
        histogram(N,linspace(0,max(N),50))
        xlabel('No objects $[\rm{\#}]$','interpreter','latex')
        ylabel('No frame occ. $[\rm{\#}]$','interpreter','latex')
        title('Number of tracked objects','interpreter','latex')
        
        % tracked number of frames
        subplot(1,2,2)
        T= histcounts(Tlink(1,C),[unique(Tlink(1,C))-1/2,max(Tlink(1,C))+1/2]) ; % reference frame n
        histogram(T,linspace(0,max(T),50))
        xlabel('No frames $[\rm{\#}]$','interpreter','latex')
        ylabel('No index occ. $[\rm{\#}]$','interpreter','latex')
        title('Number of tracked frames','interpreter','latex')
        
        drawnow
        
    end
    
end

end

