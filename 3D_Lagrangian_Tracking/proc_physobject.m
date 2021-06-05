function proc_physobject(frm_proc)
%proc_physobject Process physical object data.
%
%   Input:
%       Cpln Camera plane ellipse idenitifications
%       Fsol Feasible solution.
%   
%   Output:
%       Plink Physical links object indexing to camera indexing, including
%             triangulation data stored from Fsol, i.e. using the 
%             raytracing triangulation.
%       Pobj  Reconstructed object data using a regularized quadric
%             reconstruction. This caputures an object position size and 
%             shape and the triangulation can be compared to Fsol when
%             needed. Also the quadric as a levelset function can be used
%             to build an intensity volume in object space, that is usefull
%             for doing 3D PIV by image corrolation.
%   

%% Get globals
global folder date rec

%% Load data
Fsol=fload([folder date rec vsl 'Preproc3DTracking' vsl 'opt_Fsol.dat']);

%% loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruct the physical object data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    
    [Plink,Pobj]=rec_Pobj(Fsol,frm_proc);
    
end

%% write output and append object reconstruction

% physical solution
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'Plink.dat'],Plink,'a'); % append / rare case could append double when restarted

% physical object
fsave([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat'],Pobj,'a'); % append / rare case could append double when restarted

end

