function postproc_infotransfer
%postproc_infotransfer process information transfer by pair corrolation and
%   while doing a extension to the Procrustes analysis using least 
%   squares.
%
%   steps:
%       1) select timeshifted frame pair
%       2) select track pairs over time shift
%       3) compute Procrusters analysis using Kabsch algorithm (proof
%       omitted, might not be exact optimum but approx) find:
%           a) translation vector
%           b) rotationmatrix
%           c) rescaling
%           d) Procrustes distance
%       4) minimize the residual for the timeshift
%       5) interpolate timeshift
%       
%   @ several variabilities we trade speed over generality in underlying
%   least squares formalation, we gain clear interpretation.

%!also save the pair corrolation data in different dissimilarities
! spatio temporal corrolation data
! split in computing the spatio temporal corrolation and optimizing it

%% get globals
global folder date rec post prop pair

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer'])
end

%% Initiate files to write
PairIndex=zeros(1+1+2,0); % [ n | m | i j ]
TimeShift=zeros(1,0); % [dt]
Translation=zeros(3,0); % [t1 t2 t3]
Rotation=zeros(4,0); % [q1 q2 q3 q4]
Scaling=zeros(9,0); % [s1 ... s9]
Dissimilarity=zeros(3,0); % [d]

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for n=post.tproc(1)+(pair.cwin+pair.mran):post.tproc(2)-(pair.cwin+pair.mran) % post.tproc(1)+(cwin):post.tproc(2)-(cwin+mran)
    
    % get data frame range n
    N= Index(2,:)>=n-pair.cwin & Index(2,:)<=n+pair.cwin;
    
    % get the indexing
    indI=Index(1,N);
    
    % get frame indexing
    indN=Index(2,N);
    
    % get the positions
    posI=Position(:,N);
    
    % make sure that data is complete
    [ui,~,ic]=unique(indI);
    val=histcounts(indI,[ui-1/2,max(ui)+1/2])==2*pair.cwin+1;
    val=val(ic);
    
    % keep complete part
    indI=indI(val);
    indN=indN(val);
    posI=posI(:,val);
    
    %figure; scatter3(posI(1,:),posI(2,:),posI(3,:),[],indN,'.'); axis equal;   xlim([-2 2]); ylim([10 16]); zlim([-4 0.5]); view(2)
    
    % initiate write
    pairindex=zeros(1+2,0); % [ m | i j ]
    timeshift=zeros(1,0); % [dt]
    translation=zeros(3,0); % [t1 t2 t3]
    rotation=zeros(4,0); % [q1 q2 q3 q4]
    scaling=zeros(9,0); % [s1 ... s9]
    dissimilarity=zeros(3,0); % [d]

    % loop shifts in -m to +m
    for m=max(post.tproc(1)+pair.cwin,n-pair.mran):min(post.tproc(2)-pair.cwin,n+pair.mran) % n:n+mran
        
        % get data frame range n
        M= Index(2,:)>=m-pair.cwin & Index(2,:)<=m+pair.cwin;
        
        % get the indexing
        indJ=Index(1,M);
        
        % get frame indexing
        indM=Index(2,M);
        
        % get the positions
        posJ=Position(:,M);
        
        % make sure that data is complete
        [ui,~,ic]=unique(indJ);
        val=histcounts(indJ,[ui-1/2,max(ui)+1/2])==2*pair.cwin+1;
        val=val(ic);
        
        % keep complete part
        indJ=indJ(val);
        indM=indM(val);
        posJ=posJ(:,val);
        
        %figure; scatter3(posJ(1,:),posJ(2,:),posJ(3,:),[],indM,'.'); axis equal; xlim([-2 2]); ylim([10 16]); zlim([-4 0.5]); view(2)
        
        % make segmentation by midpoint for speed and redefine tvecij
        selN=indN==n;
        selM=indM==m;
        
        % define unique pairs posI posJ
        [pairI,pairJ]=meshgrid(indI(selN),indJ(selM)); % unique(indI),unique(indJ));
        pairIJ=unique([pairI(:) pairJ(:)],'rows')';
        
        % compute metric distance midpoint tracks on pairs
        [~,~,epI]=unique(pairIJ(1,:)); % enumerate pairI (we known that unique(pairI)=unique(indI))
        [~,~,eI]=unique(indI(selN)); % enumerate indI (thus has same enumeration)
        AdjI=( sparse((1:length(epI))',epI,ones(length(epI),1),length(epI),max(epI)) ...
            *sparse(eI,(1:length(eI))',ones(length(eI),1),max(eI),length(eI)) ) > 0 ; % replaces: AdjI=bsxfun(@eq,sparse(pairIJ(1,:))',sparse(indI(selN)));
        mposI=(AdjI*posI(:,selN)')';
        [~,~,epJ]=unique(pairIJ(2,:)); % enumerate pairI (we known that unique(pairI)=unique(indI))
        [~,~,eJ]=unique(indJ(selM)); % enumerate indI (thus has same enumeration)
        AdjJ=( sparse((1:length(epJ))',epJ,ones(length(epJ),1),length(epJ),max(epJ)) ...
            *sparse(eJ,(1:length(eJ))',ones(length(eJ),1),max(eJ),length(eJ)) )>0; % replaces: AdjJ=bsxfun(@eq,sparse(pairIJ(2,:))',sparse(indJ(selM)));
        mposJ=(AdjJ*posJ(:,selM)')';
        
        % segment pairs in metric range @ time consuption
        segIJ=sqrt(sum((mposJ-mposI).^2,1))<=pair.mdis; % segmentation selection
        pairIJ=pairIJ(:,segIJ); % pairs in range
        
        % selected data by pairs
        selI=ismember(indI,pairIJ(1,:));
        selJ=ismember(indJ,pairIJ(2,:));
        
        % define corresponding frame shift
        frmMIJ=m*ones(1,size(pairIJ,2)); % frame m pairIJ
        
        % assemble data for input to the Kabsch algorithm
        posPI=zeros(3,2*pair.cwin+1,size(pairIJ,2)); % pos pair I
        posPJ=zeros(3,2*pair.cwin+1,size(pairIJ,2)); % pos pair J
        for i=1:2*pair.cwin+1
            % selection
            selN=indN==(n-(i-1-pair.cwin)) & selI; % indN
            selM=indM==(m-(i-1-pair.cwin)) & selJ; % indM
            
            % data pos pair J
            [~,~,epI]=unique(pairIJ(1,:)); % enumerate pairI (we known that unique(pairI)=unique(indI))
            [~,~,eI]=unique(indI(selN)); % enumerate indI (thus has same enumeration)
            AdjI=( sparse((1:length(epI))',epI,ones(length(epI),1),length(epI),max(epI)) ...
                *sparse(eI,(1:length(eI))',ones(length(eI),1),max(eI),length(eI)) ) > 0 ; % replaces: AdjI=bsxfun(@eq,sparse(pairIJ(1,:))',sparse(indI(selN))); %vectorize by inproduct
            posPI(:,i,:)=(AdjI*posI(:,selN)')';
            
            % data pos pair J
            [~,~,epJ]=unique(pairIJ(2,:)); % enumerate pairI (we known that unique(pairI)=unique(indI))
            [~,~,eJ]=unique(indJ(selM)); % enumerate indI (thus has same enumeration)
            AdjJ=( sparse((1:length(epJ))',epJ,ones(length(epJ),1),length(epJ),max(epJ)) ...
                *sparse(eJ,(1:length(eJ))',ones(length(eJ),1),max(eJ),length(eJ)) )>0; % replaces: AdjJ=bsxfun(@eq,sparse(pairIJ(2,:))',sparse(indJ(selM))); %vectorize by inproduct
            posPJ(:,i,:)=(AdjJ*posJ(:,selM)')';
            
        end
        
        % perform Procrustes analysis
        [ t,R,s,d ] =procrust(posPI,posPJ,'matrix'); % vectorizing procrust will improve
        
        % write variables
        delTIJ=(frmMIJ-n)/prop.fps;
        tvecIJ=reshape(t,3,[]);
        qvecIJ=reshape(rotm2qcom(R),4,[]);
        svecIJ=reshape(permute(s,[2 1 3]),9,[]);
        disIJ=reshape(d,3,[]);
        
        %figure; plot(sort(disIJ),'.')
        
        %figure; plot(sqrt(sum(tvecIJ.^2,1)),disIJ,'.')
        
        %figure; scatter3(tvecIJ(1,:),tvecIJ(2,:),tvecIJ(3,:),[],disIJ,'.'); colorbar; xlim([0 1])
        
        %%% write variables beneath
        pairindex=cat(2,pairindex,[frmMIJ;pairIJ]); % [ m-n | i j ]
        timeshift=cat(2,timeshift,delTIJ); % [dt]
        translation=cat(2,translation,tvecIJ); % [t1 t2 t3]
        rotation=cat(2,rotation,qvecIJ); % [q1 q2 q3 q4]
        scaling=cat(2,scaling,svecIJ); % [s1 ... s9]
        dissimilarity=cat(2,dissimilarity,disIJ); % [d]

        % processing message
        disp(['processed frame pair (n,m)=(',num2str(n),',',num2str(m),')'])
    end
    
    % loop unique pairIJ to max. the corrolation by min. dissimilarity
    [~,~,ic]=unique(pairindex([2,3],:)','rows'); % infotransfer 
    for i=unique(ic)'
        % select frames and dissimilaritie to mimize over
        sel=find(ic==i);
        t=timeshift(1,sel); 
        d=dissimilarity(3,sel); % optimal fit
        
        %figure; plot(t,d,'o');
        
        % mark relevant local minima
        locm=zeros(size(t))==1;
        locm(2:end-1)=(d(2:end-1)<d(3:end) & d(2:end-1)<d(1:end-2) )...
            & ~( d(2:end-1)<10^(-6) & t(2:end-1) ==0);
        
        %hold on; plot(t(locm),d(locm),'o');
        
%         %%% write multiple local minima
%         PairIndex=cat(2,PairIndex,         [n*ones(1,nnz(locm))
%                                               pairindex(:,sel(locm))] ); % [ n | m-n | i j ]
%         TimeShift=cat(2,TimeShift,         t(locm)                    ); % [dt]
%         Translation=cat(2,Translation,     translation(:,sel(locm))   ); % [t1 t2 t3]
%         Rotation=cat(2,Rotation,           rotation(:,sel(locm))      ); % [q1 q2 q3 q4]
%         Scaling=cat(2,Scaling,             scaling(:,sel(locm))       ); % [s1 ... s9]
%         Dissimilarity=cat(2,Dissimilarity, d(locm)                    ); % [d]
%         
%         % message
%         %disp(['processed pairIJ ',num2str(i),' / ',num2str(max(ic))]) % slow
        
        if nnz(locm)>0
            % minimize the global minima
            glbm=d==min(d(locm));
            
            %hold on; plot(t(glbm),d(glbm),'o');
            
            % interpolate time shift
            t=t(bwmorph(glbm,'dilate',1));
            d=d(bwmorph(glbm,'dilate',1));
            coef=((t'.^(0:2))'*(t'.^(0:2)))\((t'.^(0:2))'*d'); % coefficients
            % hold on; plot(linspace(t(1),t(3),100),coef(1)+coef(2)*linspace(t(1),t(3),100)+coef(3)*(linspace(t(1),t(3),100).^2),'.')
            delTIJ=-coef(2)/(2*coef(3));
            disIJ=d(2);%coef(1)+coef(2)*delTIJ+coef(3)*(delTIJ^2); % can interpolate negtive values: use d(2)
            % hold on; plot(delTIJ,disIJ,'+')
            
            %%% write single global minima
            PairIndex=cat(2,PairIndex,         [n;pairindex(:,sel(glbm))] ); % [ n | m-n | i j ]
            TimeShift=cat(2,TimeShift,         delTIJ                     ); % [dt]
            Translation=cat(2,Translation,     translation(:,sel(glbm))   ); % [t1 t2 t3]
            Rotation=cat(2,Rotation,           rotation(:,sel(glbm))      ); % [q1 q2 q3 q4]
            Scaling=cat(2,Scaling,             scaling(:,sel(glbm))       ); % [s1 ... s9]
            Dissimilarity=cat(2,Dissimilarity, disIJ                      ); % [d]
            
            % message
            %disp(['processed pairIJ ',num2str(i),' / ',num2str(max(ic))]) % slow
            
        end
        
    end
    
    %figure; plot3(TimeShift,sqrt(sum(Translation.^2,1)),Dissimilarity,'.')
    
    %sel=PairIndex(4,:)~=PairIndex(3,:) & Dissimilarity<inf+0.001;
    %map=histcounts2(TimeShift(sel),sqrt(sum(Translation(:,sel).^2,1)),25);
    %figure; nor=1+0*repmat((1:size(map,1)).^2,size(map,1),1); surf((map./nor)'); view(2); shading flat
    %figure; plot3(TimeShift(sel),sqrt(sum(Translation(:,sel).^2,1)),Dissimilarity(sel),'.')
    
    % message
    disp(['post proc. info. trans. @ frm. ',num2str(n)])
    
end

toc

%% impose causality sort out the consistency between tracks
% then we would like indexing m-n to match
% ccpair=sum(PairIndex([1 2 3 4],:)==PairIndex([2 1 4 3],:),1)==2; % causal connected pair N M I J == M N J I

%figure; plot3(TimeShift(ccpair),sqrt(sum(Translation(:,ccpair).^2,1)),Dissimilarity(ccpair),'.')

%sel=ccpair;
%map=histcounts2(TimeShift(sel),sqrt(sum(Translation(:,sel).^2,1)),25);
%figure; nor=1+0*repmat((1:size(map,1)).^2,size(map,1),1); surf((map./nor)'); view(2); shading flat

%sel=PairIndex(3,:)==PairIndex(4,:) & ccpair;
%map=histcounts2(TimeShift(sel),sqrt(sum(Translation(:,sel).^2,1)),25);
%figure; nor=1+0*repmat((1:size(map,1)).^2,size(map,1),1); surf((map./nor)'); view(2); shading flat

%sel=PairIndex(3,:)~=PairIndex(4,:) & ccpair;
%map=histcounts2(TimeShift(sel),sqrt(sum(Translation(:,sel).^2,1)),25);
%figure; nor=1+0*repmat((1:size(map,1)).^2,size(map,1),1); surf((map./nor)'); view(2); shading flat

%% save
%PairIndex [ n | m | i j ]
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer' vsl 'PairIndex.dat'],PairIndex,'w');

% TimeShift [dt]
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer' vsl 'TimeShift.dat'],TimeShift,'w');

% Translation [t1 t2 t3]
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer' vsl 'Translation.dat'],Translation,'w');

% Rotation [q1 q2 q3 q4]
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer' vsl 'Rotation.dat'],Rotation,'w');

% Scaling [s1 ... s9]
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer' vsl 'Scaling.dat'],Scaling,'w');

% Dissimilarity [d]
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'InformationTransfer' vsl 'Dissimilarity.dat'],Dissimilarity,'w');

end


%%% some parameters
% post.mran=20; % range in m
% post.cwin=10; % corrolation window size intersection
% post.mdis=0.5; % maximum metric corrolation distance
