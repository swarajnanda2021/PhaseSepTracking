function [ filesel ] = sel_chkboard
%sel_chkboard Select checkerboard for further processing after image
%   processing, several use defined criteria are/can be implemented here.
%
%   Input:
%       Processed image data
%       "user defined criteria in this script"
%
%   Output:
%       filesel.mat structure containing field
%           .name Image set name
%           .dis  Files selected for distortion correction
%           .cal  Files selected for calibration

%% get globals
global folder date cal ctrl prop

%% files
fil=dir([folder date cal vsl '*_IMcnod.dat']); % keep

%% random board selection
nfil=inf; % number of files to consider
sel=1:length(fil);
while length(sel)>nfil % random board removal
    sel(randi(length(sel)))=[];
end

%% go true files
filesel=struct('name',cell(0),'dis',cell(0),'cal',cell(0));
count=1;
check=0;
for k=sel
    
    %name
    name=fil(k).name(1:end-11);
    
    %data    
    IMcnod=fload([folder date cal vsl fil(k).name]);
    
    for n=unique(IMcnod(1,:))
        N=IMcnod(1,:)==n;
        
        for c=unique(IMcnod(2,N))
            % camera selection
            C=IMcnod(2,:)==c;
            
            % chkboard
            X=IMcnod(4:5,C&N);
            
            % shape param
            xmin=min(X(1,:));
            ymin=min(X(2,:));
            xmax=max(X(1,:));
            ymax=max(X(2,:));
            
%             % area projection
%             [~,Aproj]=boundary(X');
            
            % write selection stuff
            if xmin>=ctrl.roi(1) && ... 
                    ymin>=ctrl.roi(2) && ...
                    xmax<=ctrl.roi(3) && ...
                    ymax<=ctrl.roi(4) % note some points can lie outside the image

                % check
                check=1;
                
                % name file
                filesel(count).name=name;
                
                % calibration      
%                 if strcmp(name(23:39),'Calibration_Route') 
%                     % 'calibration_route' / name(23:39) / (18:34)
                    
                    filesel(count).cal=[ filesel(count).cal
                        n c];
                    
%                 end
                                
                % distortion     strcmp(name(1:10),'2018_01_10') && 
%                 if ( sqrt(Aproj)>=(sqrt(prod(prop.res(1:2)))/3) ) || ...
%                         strcmp(name(23:43),'Distortion_Correction')
%                     %'distortion_calibration' / name(23:44) / 2 / (18:38)
                    
                    filesel(count).dis=[ filesel(count).dis
                        n c];
                    
%                 end
                
            end
        end
        
    end
    
    if 0
        % select random best board set n
        N=bsxfun(@eq,filesel(count).dis(:,1),filesel(count).cal(:,1)');
        A=sum(N,2)*sum(N,1).*N;
        S=A==max(A(:));
        S=filesel(count).dis(:,1)*filesel(count).cal(:,1)'.*S;
        v=unique(S(S>0));
        S=S==v(randi(length(v)));
        sd=sum(S,2)>0;
        sc=sum(S,1)>0;
        filesel(count).dis(~sd,:)=[];
        filesel(count).cal(~sc,:)=[];
    end
    
    % count #boards use dis or cal
    if check==1
        count=count+1;
    end
    check=0;
    
end

%% save
save([folder date cal vsl 'filesel.mat'],'filesel') 

end