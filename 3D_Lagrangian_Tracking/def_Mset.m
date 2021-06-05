function [Mset,frm_new] = def_Mset(Cpln,frm_rem,frm_new)
%def_Mset Define matched sets between camera indices
%
%   Input,
%       Cpln Camera plane data with indexing
%
%   Output,
%       Mset Indexed matched sets including triangulation data for
%            instateneous reconstruction
%
%   Note,
%       This function makes all possible stereo to n-focal camera matches
%       which lie within the DOF and ellipse identification and then
%       optimizes a greedy solution (cost sorting) to keep a tractable set.
%
%       This step can be seen as a instateneous frame tomographic
%       reconstruction for the position. However, these matched sets seed
%       object trajectory and can be modified, this is where time
%       information will later be use to improve from this instateneous
%       frame matched sets.
%   
%   TODO
%		At high [optical] density the matching can experience memory 
%		spikes. Vectorize the logical matrix & statements line 258-280.
%		Also, pre-selection a random camera coordinate to branch from
%		will resolve the memory spike at cost of needing to rebranch 
%		looping the coordinates (major) and a sub-optimal solution
%		space (minor).
%		

%% Get globals
global folder date rec prop ctrl plotting Pmat Kmat

%% Initiate data

if exist([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Mset.dat'],'file')~=0 % load
    
    % load previous timelink data
    Mset=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Mset.dat']);
    
    % maximum index
    mmax=max(Mset(1,:));
    
    % remove processed frame if already passed this step
    frm_new=frm_new(:,~ismember(frm_new,Mset(3,:)));
%     frm_new=frm_new(:,~ismember(frm_new,Cpln(3,...
%         ismember(Cpln(1,:),Mset(2,:))) ));
    
else % initiate
    
    % initiate camera links
    Mset=zeros(11,0); % [match-index camera-index n-index identif-index triangulation displacement reproj]
    
    % define maximum index
    mmax=0;
    
end

%% Match nfocal sets

%-- #) cluster focality

%-- loop frames

for n=frm_new
	
	%%%%
	%%%
	%%%%
	
	% Message
	disp(['Correspondence matching frame ',num2str(n)])
	
	%%%%%
	%%% cluster over focality
	%%%%%
	
	% loop focality
	for f=2:ctrl.lfoc(2) % prop.res(4) 
		
		%-- 0) Begin message
		
		% message
		disp(['Matching ',num2str(f),'-focal sets'])
		
		%-- 3) Loop camera sets in focality
		
		% nfocal camera sets
		c_focf=nchoosek(1:prop.res(4),f)'; % nfocal
		
		% loop number of camera links
		for i=1:size(c_focf,2)
			
			%-- 1) get data
			
			% Select time frame
			cpln=Cpln(:, Cpln(3,:)==n & ismember(Cpln(2,:),c_focf(:,i)) );
			
            % existing matched sets
            mset=Mset(:, ismember(Mset(1,:), ...
                Mset(1,ismember(Mset(4,:),cpln(1,:)) )) ); % includes sets outside camera focality
            
            % remove outside focality
            mset=mset(:,~ismember(mset(1,:),mset(1,~ismember(mset(4,:),cpln(1,:)))));
            
			%-- #) Include camera disparity in camera data
			
			% change body size by disparity allowance
			cpln(4:9, : )=contrans(cpln(4:9, : ),ctrl.dspr,'resize');
			
			% Regularize conics at camera precision for evaluation
			for c=c_focf(:,i)'
				
				g_cam=cpln(2,:)==c;
				
				Dreg=prop.ecal(c)/sqrt(abs(det(Kmat{c}/Kmat{c}(end))));
				
				cpln(4:9,g_cam)=contrans(cpln(4:9,g_cam),Dreg,'expand');
				
			end
            
			%-- 2) define adjacencies
			
            % unique indices feasible track indices
            uc=unique(cpln(1,:));%[ mset(2,:)  ]
            
			% adjacency camera data
			cpi=(1:size(cpln,2))';
			[ucp,~,cp]=unique(cpln(1,:));
			iucp=find(ismember(uc,ucp));
			cp=iucp(cp);
			C_pln=sparse(cp,cpi,ones(size(cpi)),length(uc),length(cpi));
			
            %figure; spy(C_pln)
            
            c_cam=sum(C_pln*spdiags(cpln(2,:)',0,size(cpln,2),size(cpln,2)),2)';
            
			%figure; plot(c_cam,'.')
            
			% Define existing matched set
			[um,~,ms]=unique(mset(1, : )); %ismember(mset(2,:),ftrk(1,:)) sorting in is unique
			[umc,~,mc]=unique(mset(4, : )); %ismember(mset(2,:),ftrk(1,:)) sorting in is unique
			iumc=find(ismember(uc,umc)); % uft
			mc=iumc(mc);
			M_set=sparse(mc,ms,ones(size(ms)),length(uc),length(um));
			
			%figure; spy(M_set)
			%figure; plot(sum(M_set,1),'.')
			
			%-- #) Match correspondences
			
			if f==2 % initiate stereofocal
                
				%-- #) stereo focal linking
				
				% fundamental/essential matrix F m2s
				Fmat=pmat2fmat(Pmat{c_focf(1,i)},Pmat{c_focf(2,i)});
				
                % indices data in views and frames
                l=cpln(1, cpln(2,:)==c_focf(1,i) );%find(ismember(uc,))
                h=cpln(1, cpln(2,:)==c_focf(2,i) );%find(ismember(uc,))
                
				% get conics
				Cm=cpln(4:9, cpln(2,:)==c_focf(1,i) );
				Cs=cpln(4:9, cpln(2,:)==c_focf(2,i) );
				
				% camera coordinates current frame
				xm = cvec2eshp( Cm );
				xs = cvec2eshp( Cs );
				
				% compute lines of sights current frame
				lm=homc2inhc(Fmat*inhc2homc(xm)); % line from m
				ls=homc2inhc(Fmat'*inhc2homc(xs)); % line from s
				
				% point line distance on ellipse
				Dms=ellipselindist(Cm,ls)/2;%;%
				Dsm=ellipselindist(Cs,lm)/2;%;%
				
				% segment disparity
				[di,dj]= find( Dms<=1 & Dsm'<=1 ); % (Dms+Dsm')/2<=1 % Dms<=1 | Dsm'<=1 
				
				% write at indexing
				l=l(di);
				h=h(dj);
				
                % write new stereo sets
                [ulh,~,lh]=unique([l,h]);
                iulh=find(ismember(uc,ulh));
                lh=iulh(lh);
                M_new=sparse(lh,[1:length(l) 1:length(h)],ones(size(lh)),length(uc),length(l));
                
                %figure; spy(M_new)
                %figure; plot(sum(M_new,1),'.')
                
			else % cluster trifocal quafocal ... nfocal
                
                % sets of selected feasible tracks in previous focality
                g_focfm1=sum(M_set,1)==(f-1) ; % (g_pln,:)(ismember(c_ftrk,c_focn(:,i)),:)
                
                % indexing
                gind=find( g_focfm1 ); % g_cset &
                
                %figure; plot(g_focfm1,'.'); nnz(g_focfm1)
                
				% previous focality within current focality
				c_focfm1=nchoosek( c_focf(:,i) ,f-1)';
				
				% combinations subfocality
				c_subfoc=nchoosek( 1:size(c_focfm1,2) , 2)'; % pick two
				
                % initiate and match adjacency
                A=sparse( nnz( g_focfm1),nnz( g_focfm1) ); %g_cset &
                
                %figure; spy(A)
                
                % loop combinations
                for j=1:size(c_subfoc,2)
                    
                    % sets
                    g1_focfm1=sum(M_set( ismember( c_cam, c_focfm1(:,c_subfoc(1,j)) ) , g_focfm1 ) ,1) == (f-1); % g_cset
                    g2_focfm1=sum(M_set( ismember( c_cam , c_focfm1(:,c_subfoc(2,j)) ) , g_focfm1 ) ,1) == (f-1); %
                    
                    % indices
                    g1=find(g1_focfm1);
                    g2=find(g2_focfm1);
                    
                    %figure; plot(g1_focfm1,'.')
                    %figure; spy(M_set & g1_focfm1)
                    %figure; spy(M_set & g2_focfm1)
                    
%                     B=double( g1_focfm1 & M_set(:, g_focfm1 ) )'*double( g2_focfm1 & M_set(:, g_focfm1 ) ) == (f-2) ; %
                    
                    % initiate and match adjacency
                    [bi,bj,bk]=find( ( M_set(:, gind(g1_focfm1) )'*M_set(:, gind(g2_focfm1) ) ) == (f-2) );
                    B=sparse(g1(bi),g2(bj),bk',size(A,1),size(A,2))>0.5; % mild speed up
                    
                    %figure; spy(B)
                    
                    % append subbranch that reduces the number of combinations
                    A =A | B | B' ;
                    
                end % j
                
                %figure; spy(A)
                %figure; plot(sum(A,1),'.') % (1:floor(size(A,1)/2),:)
                
                % pick three subfocality to propagate connection
                c_subfoc=nchoosek( 1:size(c_focfm1,2) , 3)'; % pick two
                
                % propagate and clean inconsistency
                p=zeros(2,0);
                for j=1:size(c_subfoc,2) % slow
                    
                    % sets
                    g1_focfm1=sum(M_set( ismember( c_cam , c_focfm1(:,c_subfoc(1,j)) ) , g_focfm1 ) ,1) == (f-1); %
                    g2_focfm1=sum(M_set( ismember( c_cam , c_focfm1(:,c_subfoc(2,j)) ) , g_focfm1 ) ,1) == (f-1); %
                    g3_focfm1=sum(M_set( ismember( c_cam , c_focfm1(:,c_subfoc(3,j)) ) , g_focfm1 ) ,1) == (f-1); %
                    
                    g1=find(g1_focfm1);
                    g2=find(g2_focfm1);
                    g3=find(g3_focfm1);
                    
                    %figure; plot(g1_focfm1,'.')
                    %figure; plot(g2_focfm1,'.')
                    %figure; plot(g3_focfm1,'.')
                    %figure; plot(sum(M_set( ismember( c_cam , c_prop(:,k)) , g_focfm1 ) ,1),'.')
                                        
                    % first
%                     P=A(g3_focfm1,g1_focfm1) & double(A(g3_focfm1,g2_focfm1))*double(A(g2_focfm1,g1_focfm1)); % propagate
%                     [pi,pj]=find(P);
%                     p=cat(2,p,[g3(pi) ; g1(pj)],[g1(pj) ; g3(pi)]);
                    [ai,aj]=find(A(g3_focfm1,g1_focfm1));
                    ai=reshape(ai,1,[]);
                    aj=reshape(aj,1,[]);
                    val=sum(A(g2_focfm1,g3(ai)) & A(g2_focfm1,g1(aj)),1)>0.5; % propagate
                    p=cat(2,p,[g3(ai(val)) ; g1(aj(val))],[g1(aj(val)) ; g3(ai(val))]); % extra speed up (mild)
%                     A(g3_focfm1,g1_focfm1) = P  ;
%                     A(g1_focfm1,g3_focfm1) = P' ;
                    
                    %figure; spy(A(g2_focfm1,g1_focfm1))
                    %figure; spy(A(g3_focfm1,g2_focfm1))
                    %figure; spy(A(g3_focfm1,g2_focfm1)*A(g2_focfm1,g1_focfm1))
                    %figure; spy(A(g3_focfm1,g1_focfm1) & A(g3_focfm1,g2_focfm1)*A(g2_focfm1,g1_focfm1))
                    %figure; spy(A(g3_focfm1,g1_focfm1))
                    
                    % second
%                     P=A(g2_focfm1,g1_focfm1) & double(A(g2_focfm1,g3_focfm1))*double(A(g3_focfm1,g1_focfm1)); % propagate
%                     [pi,pj]=find(P);
%                     p=cat(2,p,[g2(pi) ; g1(pj)],[g1(pj) ; g2(pi)]); % speed up
                    [ai,aj]=find(A(g2_focfm1,g1_focfm1));
                    ai=reshape(ai,1,[]);
                    aj=reshape(aj,1,[]);
                    val=sum(A(g3_focfm1,g2(ai)) & A(g3_focfm1,g1(aj)),1)>0.5; % propagate
                    p=cat(2,p,[g2(ai(val)) ; g1(aj(val))],[g1(aj(val)) ; g2(ai(val))]); % extra speed up (mild)
%                     A(g2_focfm1,g1_focfm1) = P  ;
%                     A(g1_focfm1,g2_focfm1) = P' ;
                    
                    % third
%                     P=A(g2_focfm1,g3_focfm1) & double(A(g2_focfm1,g1_focfm1))*double(A(g1_focfm1,g3_focfm1)); % propagate
%                     [pi,pj]=find(P);
%                     p=cat(2,p,[g2(pi) ; g3(pj)],[g3(pj) ; g2(pi)]); % speed up
                    [ai,aj]=find(A(g2_focfm1,g3_focfm1));
                    ai=reshape(ai,1,[]);
                    aj=reshape(aj,1,[]);
                    val=sum(A(g1_focfm1,g2(ai)) & A(g1_focfm1,g3(aj)),1)>0.5; % propagate
                    p=cat(2,p,[g2(ai(val)) ; g3(aj(val))],[g3(aj(val)) ; g2(ai(val))]); % extra speed up (mild)
%                     A(g2_focfm1,g3_focfm1) = P  ;
%                     A(g3_focfm1,g2_focfm1) = P' ;
                    
                end
                
                P=sparse(p(1,:),p(2,:),ones(1,size(p,2)),size(A,1),size(A,2)); % need to iter A=P apparently
                
                %figure; spy(P)
                %figure; plot(sum(P,1),'.')
                
                % define pairing in ones of the set combination, here last
                [pi,pj]=find(tril(P,-1)); %find( A & B ); % 
                
                M_new=M_set(:,gind(pi)) | M_set(:,gind(pj)); % direct
                
                %figure; spy(M_new)
                %figure; plot(sum(M_new,1),'.')
                
                [mi,~]=find(M_new);
                
                miu=unique(reshape(mi,f,[])','rows')';
                
                mju=repmat(1:size(miu,2),f,1);
                
                M_new=sparse(miu(:),mju(:),ones(size(miu(:))),length(uc),size(mju,2));
                
                %figure; spy(M_new)
                %figure; plot(sum(M_new,1),'.')
                
			end
			
            %-- #) Triangulate matched sets and check validity
            
            % Assemble triangulation data
            xdat=cell(1,length(c_focf(:,i)));
            for j=1:length(xdat)
                
                g_cam=cpln(2,:)==c_focf(j,i);
                
                % Get midpoints for triangulation
                x=cvec2eshp(cpln(4:9,g_cam));
                
                % adjacencie camera data and solution
                A=C_pln(:,g_cam)'*M_new;
                
                % find
                [ai,aj]=find(A);
                
                xdat{j}=nan*zeros(3,size(M_new,2)); %0
                
                xdat{j}(:,aj)=inhc2homc(x(:,ai));
                
            end
            
            % Triangulate,~,r
            X=objtriang(xdat,Pmat(c_focf(:,i))); % lintriang
            
            %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
            
            %-- #) Find displacement vector
            
            % Assemble triangulation data
            xdat=cell(1,length(c_focf(:,i)));
            for j=1:length(xdat)
                
                g_cam=cpln(2,:)==c_focf(j,i);
                
                % Get midpoints for triangulation
                x=cvec2eshp(cpln(4:9,g_cam))+cpln(10:11,g_cam);
                
                % adjacencie camera data and solution
                A=C_pln(:,g_cam)'*M_new;
                
                % find
                [ai,aj]=find(A);
                
                xdat{j}=nan*zeros(3,size(M_new,2)); %0
                
                xdat{j}(:,aj)=inhc2homc(x(:,ai));
                
            end
            
            % Triangulate,~,r
            Xacc=objtriang(xdat,Pmat(c_focf(:,i))); % lintriang
            
            %figure; plot3(Xacc(1,:),Xacc(2,:),Xacc(3,:),'.'); axis equal; view(2)
            
            % displacement
            V=Xacc-X;
            
 			%-- #) new matched sets
			
			% find
			[fi,fj]=find(M_new);
			
			% pairing new camera data in tracking
			mnew=full([fj' % new index set
                c_cam(fi) % camera index
				n*ones(1,length(fi)) % end % n*ones(1,length(fj)) % start
				uc(fi') % index track uft
                X(:,fj) % triang
                V(:,fj) % displacement
				zeros(1,length(fj))]); % make full
			
			% Define set mapping
			[umn,~,mn]=unique(mnew(1,:)); % sorting in is unique
			G_new=sparse(mn,1:length(mn),ones(size(mn)),length(umn),length(mn));
			
			%figure; spy(G_new)
			
            %-- #) Check reprojection error
            
			% initiate exclusion set
			g_excl=false(1,size(mnew,2));
			
            % loop cameras
            for c=c_focf(:,i)'
                
                % camera set
                g_cam=cpln(2,:)==c;
                
                % new sets
                g_new=ismember(mnew(4,:),cpln(1,g_cam));
                dum=find(g_new);
                
                % find adjacency
                [ai,aj]=find( C_pln(:,g_cam)'*M_new*G_new(:,g_new) ); % can speed up outside t
                
                % loop image displacement
                for t=0:ctrl.tres
                    
                    % conic data
                    con=contrans(cpln(4:9,g_cam),t*cpln(10:11,g_cam),'displace');
                    
                    % velocity
                    V=mnew(8:10,g_new);
                    
                    % triangulation data
                    X=mnew(5:7,g_new)+t*V;
                    
                    % camera reference
                    X=Pmat{c}*inhc2homc(X);
                    
                    %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
                    
                    % segment dof
                    g_excl(dum(aj))=g_excl(dum(aj)) | X(3,aj)<ctrl.ldof(1) | X(3,aj)>ctrl.ldof(2);
                    
                    %figure; plot3(X(1,!~g_excl),X(2,~g_excl),X(3,~g_excl),'.'); axis equal; view(2)
                    
                    % project coordinates
                    x=homc2inhc( X ); % empty cell can comprimise
                    
                    %figure(1); plot(xp(1,:),xp(2,:),'.'); drawnow
                    
                    % compute ellipse distance
                    d=ellipsedist(con(:,ai),x(:,aj),'vector');
                    
%                     figure(f)
%                     subplot(size(c_focf,2),2,2*(i-1)+1)
%                     hold on
%                     plot(d(d<5),'.')
%                     hold off
%                     sp=subplot(size(c_focf,2),2,2*(i-1)+2);
%                     hold on
%                     histogram(d(d<5),50)
%                     plot(median(d(d<5))*[1 1],sp.YLim,'--g','LineWidth',1)
%                     plot(median(d(d<5))*[1 1]+median(abs(d(d<5)-median(d(d<5))))*[3 3 ],sp.YLim,'--r','LineWidth',1)
%                     hold off; drawnow
                    
%                     % reprojection to optimize within ellipse span (more convex with respect to anisotopy & midpoint)
%                     xm=cvec2eshp(con);
%                     r=sqrt(sum((xm(:,ai)-x(:,aj)).^2,1));
                    
                    % update exclusion
                    g_excl(dum(aj))=g_excl(dum(aj)) | d > 1;
                    
                    %figure(1); plot(xp(1,d <=1 & Xl(3,:)>=ctrl.ldof(1) & Xl(3,:)<=ctrl.ldof(2  )),...
                    %    xp(2,d <=1 & Xl(3,:)>=ctrl.ldof(1) & Xl(3,:)<=ctrl.ldof(2)),'.'); drawnow
                    
%                     % Swept reprojection error in object space
%                     d=d.*sqrt(sum(V(:,aj).^2,1));
                    
                    % update total disparity matched set
                    mnew(11,dum(aj))=mnew(11,dum(aj))+d/(2*ctrl.tres+1); % d; average value
                    
                end % t
                
            end % c
            
%             Davg=mnew(11,:)*G_new';%./sum(G_new,2)';
%             g_excl=g_excl | (Davg>1)*G_new;
            
            % full exclusion
            g_excl=ismember(mnew(1,:),mnew(1,g_excl)); % 
            
%             % iterate outliers
%             iter_cnt=0; % iteration count
%             outl_cnt=0; % outlier count
%             while outl_cnt~=nnz(g_excl) && iter_cnt<1
%                 
%                 % count
%                 iter_cnt=iter_cnt+1; % iterations
%                 outl_cnt=nnz(g_excl); % outliers
%                 
%                 % loop cameras
%                 for j=1:length(c_focf(:,i))
%                     
%                     % camera
%                     c=c_focf(j,i)';
%                     
%                     % camera set
%                     g_cam=cpln(2,:)==c;
%                     
%                     % new sets
%                     g_new= ~g_excl & ismember(mnew(2,:),cpln(1,g_cam));
%                     
%                     % disparity
%                     d=mnew(11,g_new);
%                     
%                     % outlier filter camera
%                     [~,~,bnds]=stattest(d,'median',3,1); % iter one time
%                     g_excl(g_new)=g_excl(g_new) | d<bnds(1) | d>bnds(2);
%                     
%                     figure(f)
%                     sp=subplot(size(c_focf,2),size(c_focf,1),length(c_focf(:,i))*(i-1)+j);
%                     hold on
%                     histogram(d,50)
%                     plot(bnds(1)*[1 1],sp.YLim,'--r','LineWidth',1)
%                     plot(bnds(2)*[1 1],sp.YLim,'--r','LineWidth',1)
%                     hold off; drawnow
%                     
%                 end % c
%                 
%                 % remove those indices
%                 g_excl=ismember(mnew(1,:),mnew(1,g_excl));
%                 
%             end
            
            %figure; plot(mnew(8,:),'.')
            %figure; plot(g_excl,'.')
            
            % remove exclusion
            mnew=mnew(:,~g_excl);%ismember(mnew(1,:),mnew(1,g_excl))
            
 			%-- #) Write new feasible matched sets
			
            % write new feasible matches
            [dumm,~,dum]=unique(mnew(1,:)); % quick fix
            mnew(1,:)=mmax+dum; % assign index
            mmax=mmax+length(dumm);%size(mnew,2)/f;
            
			%-- #) Write matched sets
			
			% write matched links
			Mset=cat(2,Mset, mnew );
			
			%-- #) End message
			
            % keep sets complete in previous focality
            disp(['Matched ',num2str(size(mnew,2)/f),' sets'])
            
		end % i
		
	end % f
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Optimize Frame Instateneous Matched Sets %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%-- 1) get data
	
	% Select time frame
	cpln=Cpln(:, Cpln(3,:)==n );
    
	% existing matched sets
	mset=Mset(:, ismember(Mset(1,:), ...
		Mset(1,ismember(Mset(4,:),cpln(1,:)) )) ); % includes sets outside camera focality
	
	%-- #) Open for reassessment
	
	% open selected matched sets
	Mset=Mset(:, ~ismember(Mset(1,:),mset(1,:)) );
	
	%-- %) define adjacencies
    
	% unique indices feasible track indices
	uc=unique(mset(4,:));
	
    % Define set mapping
    [um,~,ms]=unique(mset(1,:)); % sorting in is unique
    G_set=sparse(ms,1:length(ms),ones(size(ms)),length(um),length(ms));
    
    %figure; spy(G_set)
    
	% Define existing matched set
	[um,~,ms]=unique(mset(1, : )); %ismember(mset(2,:),ftrk(1,:)) sorting in is unique
	[umc,~,mc]=unique(mset(4, : )); %ismember(mset(2,:),ftrk(1,:)) sorting in is unique
	iumc=find(ismember(uc,umc)); % uft
	mc=iumc(mc);
	M_set=sparse(mc,ms,ones(size(ms)),length(uc),length(um));
	
    %     X=mset(5:7,ia); % ia
    
	%figure; spy(M_set)
	%figure; plot(sum(M_set,1),'.')
    %figure; plot(sum(M_set,2),'.')
    %figure; plot3(X(1,:),X(2,:),X(3,:),'.'); axis equal; view(2)
    
	%-- #) Define reward and cost functions
    
    % weight matched cameras
    cf=sum( M_set > 0.5 ,1); % average
    
    % figure; plot(cf,'.'); drawnow
    
    % disparity triangulation - this option is not better:-exp(-4*mset(11,:)^2) (as a brightness model)
    df=mset(11,:)*G_set'./sum(G_set,2)'; %.^2 2-norm (notworking sat..*sqrt(sum(mset(8:10,:).^2,1)))
    %max( G_set*spdiags(mset(11,:)',0,size(mset,2),size(mset,2)) ,[],2)';
    
	% figure; plot(df,'.'); drawnow
	% figure; histogram(df,50); drawnow
	
	%-- 4) Solve greedy optimizization
	
    % find open set
    g_ope=cf>=ctrl.lfoc(1); % true(1,size(M_set,2));% & cf>2;
    
    % initiate solution
    g_sol=false(1,size(M_set,2));
    
%     % optimize camera cover, reward, disparity
%     while nnz(g_ope)~=0
%         
%         % reward function free coordinates
%         rf=sum(M_set(sum(M_set(:,g_sol),2)==0,:),1);
        
        % figure; plot(rf,'.'); drawnow
        
        % initiate indices open solution set
        sol=find(g_ope);
        
        % sort global cost function descending ( -mean(df(g_sol)) )rf(sol)'  -of(sol)' 
        [~,ind]=sortrows([cf(sol)' -df(sol)' fliplr(1:length(sol))'],'descend');% rf(sol)' 
        
        % sort solution to cost
        sol=sol(ind);
        
        %figure; plot(df(sol),'.')
        %figure; plot(cf(sol),'.')
        
        %figure; spy(M_set(:,sol))
        %figure; plot(sum(M_set(:,sol),1),'.')
        %figure; plot(sum(M_set(:,sol),2),'.')
        
        % solve optimal cost by maximize the sortation
        [val,subsol]=max( double( M_set(:,sol)>0.5 )*spdiags(fliplr(1:length(sol))',0,length(sol),length(sol)) ,[],2); % first unique el
        
        % sub solution
        subsol=subsol(val>0);
        
        % count subsolution to find [most] indepence
        unsubsol=unique(subsol);
        
%         % count occurence subsolution
%         cnt=histcounts(subsol',[unsubsol'-1/2 max(unsubsol)'+1/2]);
%         
%         %figure; plot(cnt,'.')
%         %figure; plot(sum(M_set(:,sol(unsubsol)),1),'.')
%         
%         % select unique non-overlapping solution
%         sel=cnt==sum(M_set(:,sol(unsubsol)),1);
        sol=sol(unsubsol);%(sel));
        
        % solution
        g_sol(sol)=1;%g_sol | ismember(um,um(sol)); 
        
        % figure; histogram(df(g_sol),50); drawnow
        
%         % find complete overlapping
%         g_ope=cf>=ctrl.lfoc(1) & full(sum(M_set(sum(M_set(:,g_sol),2)==0,:),1)>0.5);% & cf>2;
%         
%     end
    
    % write solutions
    mset=mset(:, ismember(mset(1,:),um(g_sol)) );
    
	% figure; plot(df(g_sol),'.'); drawnow
	% figure; histogram(df(g_sol),50); drawnow
	
    %--# plotting
    if strcmp(plotting,'on')
        
        % plot triangulation
        figure; plot3(mset(5,:),mset(6,:),mset(7,:),'.'); axis equal; view(2); drawnow
        
        % plot triangulation
        hold on
        quiver3(mset(5,:),mset(6,:),mset(7,:),mset(8,:),mset(9,:),mset(10,:),0,'.')
        plot3(mset(5,:)+mset(8,:),mset(6,:)+mset(9,:),mset(7,:)+mset(10,:),'.'); axis equal; view(2); drawnow
        hold off
%         xlim([-3 3]); ylim([5 25]); zlim([-5 1])
        
        % plot camera cover
        [dum,~,dum2]=unique(mset(1,:));
        cnt=histcounts(mset(1,:),[dum-1/2 max(dum)+1/2]);
        cnt=cnt(dum2);
        figure; scatter3(mset(5,:),mset(6,:),mset(7,:),[],cnt,'.'); colorbar; axis equal; view(2); drawnow
%         xlim([-3 3]); ylim([5 25]); zlim([-5 1])
        
%         % plot triangulation
%         figure; scatter3(X(1,g_sol),X(2,g_sol),X(3,g_sol),[],cf(g_sol),'.'); colorbar; axis equal; view(2); drawnow
%         figure; scatter3(X(1,g_sol),X(2,g_sol),X(3,g_sol),[],df(g_sol),'.'); colorbar; axis equal; view(2); drawnow
        
    end
    
	%-- #) Write matched sets
	
	% write matched links
	Mset=cat(2,Mset, mset );
	
	%-- #) End message
	
	% end message
	disp(['Optimized ',num2str(nnz(g_sol)),' out of ',num2str(length(g_sol)),' matched sets'])
	
end % n

%% Remove data 

% trash
Mset=Mset(:,~ismember(Mset(3,:),frm_rem));
%     Mset=Mset(:,ismember(Mset(2,:),Cpln(1,:))); % remove relations

%% plot results
if strcmp(plotting,'on')
    
    %statistics
    
end

end

