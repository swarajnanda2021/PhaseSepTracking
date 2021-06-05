function postproc_freqanalysis
%postproc_freqananalysis
%       
%   [under construction]
%   

%% get globals
global folder date rec post prop traj

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'FrequencyAnalysis'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'FrequencyAnalysis'])
end

%% Initiate files to write
!Freq=;

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for n=post.tproc(1)+(post.cwin+post.mran):post.tproc(2)-(post.cwin+post.mran) % post.tproc(1)+(cwin):post.tproc(2)-(cwin+mran)
    
    % get data frame range n
    N= Index(2,:)>=n-post.cwin & Index(2,:)<=n+post.cwin;
    
    % get the indexing
    ind=Index(1:2,N);
    
    % get the positions
    pos=Position(:,N);
    
	!...
	
    % message
    disp(['post proc. info. trans. @ frm. ',num2str(n)])
    
end

toc

%% save

%PairIndex 
fsave([folder date rec vsl 'Postproc3DTracking' vsl 'FrequencyAnalysis' vsl 'WindowedFrequencies.dat'],Freq,'w');

end
