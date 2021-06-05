Koen Muller

This code includes all the files for time-resolved three-dimensional lagrangian tracking.

This code tracks multiple objects from the image-plane as ellisoids, it requires time-resolved data for optimal prestations.

run main_3D_LPT for your project, add processing settings in proc_set, or impose variables externally from batch file

Some notes
 - Image processing: This sets what can be reconstructed and therefore dominates the quality of the reconstruction.
 - Trying different settings: 
    + first perform small time span 100 frame
    + then optimize settings, play around with them,
	- different sets of image filters
	- small and large trajectory fits
    + process data with half window time steps, in case of 100 frame window say 50 frames step.
 - Check image and identication quality over different views
 - etc.. to be added

Miscelaneous
 - This code is design for working with N-views (N>=2) but has been developed for four views.

Functionality and updates:

Revisions in the 3D_Lagrangian_Tracking. 1) Improments to backtracking, backtracking is now the same algorithmic implementation as forward tracking, improving the track quality. 2) Extended the feasible solutions branching into a 6-step algorithm: seed, branch, split, evaluate, merge, locally optimize tajectories. 3) brn_Fsol can now also displace an ellipse into the next frame in case there is no identification to add to the track. 4) In addition to the trajectory correction using the Savintsky-Golay filter we now also correct and smooth the ellipse shape, stabilizing the performance when in occlusion. 5) spl_Fsol.m Added a trajectory splitting script when violation the reprojection error, not just giving up that camera view but combinatoric choosing the best n-1 focal subset to branch into.