function postproc_trackqualitymetric

%% Postprocessing code that calculates the following quantities:
% 1) Total object density;
% 2) Ghost particle density; (For all track sizes)
% 3) Occlusion probability; (Estimated only from the Silhouette Sketches)


% get globals
global folder date rec prop ctrl post traj

% general trajectories
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']); % second row is time, first row is index
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat']); % physical time of the position
Pos=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']); % just positions of particles
Quadric_axes = fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'EllipsoidAxis.dat']); % Particle quadrics
Accuracy = fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accuracy.dat']); % Particle quadrics
FitRes = fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'FitResidual.dat']); % Polyfit residuals
%% Uncertainty of the polynomial fit

if 1
%     
%     pd = fitdist(FitRes','halfnormal');
%     x_values = linspace(0.00001,0.3,100);
%     y = pdf(pd,x_values);
%     
    thresh=0.0001;
    figure(69)
    hold all
    histogram(log(FitRes(FitRes>thresh)),'edgecolor','none')
%     plot(x_values,y,'LineWidth',2)
    
    xlabel('$\log(\epsilon_{f})$ (mm)','interpreter','latex','fontsize',20)
    ylabel('[\#]','interpreter','latex','fontsize',20)
    
    vline(mean(log(FitRes(FitRes>thresh))),'k-','\mu_{\epsilon_f}')
    vline(mean(log(FitRes(FitRes>thresh)))-2*std(log(FitRes(FitRes>thresh))),'k--','(\mu-2\sigma)_{\epsilon_f}')
    vline(mean(log(FitRes(FitRes>thresh)))+2*std(log(FitRes(FitRes>thresh))),'k--','(\mu+2\sigma)_{\epsilon_f}')
    
    legend('fit residuals','interpreter','latex','fontsize',20)
    
    box on
    set(gca,'linewidth',1)
%     xlim([0.0001 0.2])
    
    
    
%     set(gca,'xscale','log')
    saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'Fit_residual.fig'])
    
    % save bias and 2 sigma bound of random errors
    unc_fit.bias = exp(mean(log(FitRes(FitRes>thresh))));
    unc_fit.random.low = exp(mean(log(FitRes(FitRes>thresh)))-2*std(log(FitRes(FitRes>thresh))));
    unc_fit.random.high = exp(mean(log(FitRes(FitRes>thresh)))+2*std(log(FitRes(FitRes>thresh))));
    
end

%% Ghost track/outlier density

if 1
    t_spurious = zeros(1,0);
    for n=post.tproc(1):post.tproc(2)
        n

        % segment from Pos the particles that are at tstep n
        N_t=Index(2,:)==n; % pick out those particle positions at the current tstep
        N_x=Pos(1,:)>post.dom(1,1) & Pos(1,:)<post.dom(1,2);
        N_y=Pos(2,:)>post.dom(2,1) & Pos(2,:)<post.dom(2,2);
        N_z=Pos(3,:)>post.dom(3,1) & Pos(3,:)<post.dom(3,2);

        % segment from N those outside the bounds of the roi
        Pos_tstep = Pos(:,N_t & N_x & N_y & N_z);
        Ind_tstep = Index(1,N_t & N_x & N_y & N_z);
        %% Part 1: Total objects

        C_n(n) = size(Pos_tstep,2);

        
        
        if 1
            %% Part 2: Ghost particle density

            load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(n) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox')
            Vox_int = permute(Vox_int,[2 1 3]);
            
            % Query all particle positions for the visual hull voxel intensity

            Vq = interp3(X_vox,Y_vox,Z_vox,Vox_int,Pos_tstep(2,:),Pos_tstep(1,:),Pos_tstep(3,:),'nearest'); % use nearest neighbor to prevent interpolation
            
            % Find all track indices that have a particle inside the cavity
            t_spurious = cat(2,t_spurious,Ind_tstep(find(Vq==1)));
            
            % Cavity volume
            V_cav(n) = sum( Vox_int(:).*(X_vox(1,2,1)-X_vox(1,1,1)).^3 );
            
            % spurious tracks per unit volume per unit timestep
            rho_ghost_tracks(n) = length(unique(Ind_tstep(find(Vq==1))))./V_cav(n); % ghost tracks per cubic mm
            
            % Look for all intensities that have a value of 1
            
            C_n_ghost(n) = length(find(Vq==1));

            Prob_ghost(n) = C_n_ghost/C_n;
            
            

        end
    end
    
    % ghost track density 
    rho_ghost = (length(unique(t_spurious))./length(unique(Index(1,:)))) ;
    
    
    
    
end



disp('Finished calculated real/ghost particle density')


%% Plot Ghost particle probability
if 1
    
    figure(4)
    plot(Prob_ghost.*100)
%     plot(rho_ghost_tracks)
    xlabel('Timestep','interpreter','latex','fontsize',20)
    ylabel('$\%$ Ghost Particles','interpreter','latex','fontsize',20)
    set(gca,'linewidth',1)
    saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'Ghost_density.fig'])
    
    figure(420)
%     plot(Prob_ghost.*100)
    histogram(rho_ghost_tracks./length(V_cav))
%     xlabel('Timestep','interpreter','latex','fontsize',20)
    xlabel('Ghost Tracks/$mm^3$/tstep','interpreter','latex','fontsize',20)
    ylabel('[\#]','interpreter','latex','fontsize',20)
    set(gca,'linewidth',1)
    saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'Ghost_track_density.fig'])
    
    
    outlier_stat.rhogparts = Prob_ghost.*100; %ghost particle density (in percentage)
    outlier_stat.rhogtracks = (rho_ghost_tracks./length(V_cav)); %ghost tracks /mm^3 /tstep (cavity region)
    
    
end



%% Calculation of the theoretical Probability Density Function 

% Current state of derivation: Prob(L>=Nf) = Prob(L>=Nf)_jerry * Prob(Not Ghost) * Prob(Not occluded)
if 0
    % Take average of all occlusion rates
    avg_occlusion = sum(Occlusion_ratio,2)/size(Occlusion_ratio,2);
    % Triangulation possibilities
    triang_perms = nchoosek(1:5,3);
    % probability of occlusion by 3-focal set
    trifocal_nonocclusion = mean(prod(1 - avg_occlusion(triang_perms),2)); % averaged over each trifocal set as they are independent


    Prob_not_occluded = trifocal_nonocclusion;



    Prob_not_ghost = 1 - mean(Prob_ghost);
end


Nf=1:50;
epsilon=[0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
for i=1:length(epsilon)
    Prob_trackOver_n = epsilon(i).^Nf.*(1 - exp(-mean(C_n)));

    Prob(i,:) = Prob_trackOver_n;% .* Prob_not_occluded .* Prob_not_ghost;

end


for i=1:size(Prob,1)
    Prob(i,:) = (Prob(i,:)./Prob(i,1));
end


disp('Finished estimating theoretical probability')

%% Calculate track-length histogram

% ( track index | indentification index | camera index | frame number | camera conic | camera velocity | object position | object velocity | res )
Plink = fload([folder date rec vsl 'Preproc3DTracking' vsl 'Plink.dat']);
% ( identification index | camera index | frame number | conic vector )
Cpln = fload([folder date rec vsl 'Preproc3DTracking' vsl 'Cpln.dat']);
indIdentPlink = ismember(Index,Plink(1,:)); % From track indices, find identification index
Pobj=fload([folder date rec vsl 'Preproc3DTracking' vsl 'Pobj.dat']); % second row is time, first row is index
%%
% Long tracks
l = unique(Index(1,:));
L=histcounts(Index(1,:),[l-1/2,max(l)+1/2]);
L2=histcounts(L,1:30);

for iii=1:length(L2)
    Prob_exp(iii) = sum(L2(iii:end))/sum(L2);
end



% Long tracks
[~,ia,~]=unique(Pobj([1 2],:)','rows');
l22=unique(Pobj(1,ia));
L22=histcounts(Pobj(1,ia),[l22-1/2,max(l22)+1/2]);
L222=histcounts(L22,1:30);

for iii=1:length(L222)
    Prob_exp_raw(iii) = sum(L222(iii:end))/sum(L222);
end


% Prob_exp = L2./L2(1);

%%
if 1
    figure(1)
    hold all
    
    plot(Prob_exp_raw,'k+')
    plot(Prob_exp,'bs')
    plot(Prob','k--')
    
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    
    xlabel('$N_f$','interpreter','latex','fontsize',20)
    ylabel('Prob($L\geq N_f$)','interpreter','latex','fontsize',20)
    legend('Raw tracks','Processed tracks','interpreter','latex','fontsize',20)
    box on
    xlim([1 30])
    saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'Track_quality.fig'])
%     legend('M37',['1-P(occlusion)=' num2str(Prob_not_occluded)],'interpreter','latex','fontsize',20)
end

disp('Finished estimating theoretical and experimental probability')


%% Size of the reconstructed particle

if 1
    
   figure(2)
   hold all
   histogram(Quadric_axes(1,:),'FaceColor','r','EdgeColor','none')
   histogram(Quadric_axes(2,:),'FaceColor','b','EdgeColor','none')
   histogram(Quadric_axes(3,:),'FaceColor','g','EdgeColor','none')
   ylabel('[\#]','interpreter','latex','fontsize',20)
   xlabel('Particle axis length ($mm$)','interpreter','latex','fontsize',20)
   xlim([ 0 0.25])
   vline(0.05,'k-')
   vline(0.04,'k--')
   vline(0.06,'k--')
   box on
   saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'Object_size.fig'])
   
   figure(21)
   hold all
   histogram(Quadric_axes(1,:)-0.05,'FaceColor','r','EdgeColor','none')
   histogram(Quadric_axes(2,:)-0.05,'FaceColor','b','EdgeColor','none')
   histogram(Quadric_axes(3,:)-0.05,'FaceColor','g','EdgeColor','none')
   ylabel('[\#]','interpreter','latex','fontsize',20)
   xlabel('Particle size error ($mm$)','interpreter','latex','fontsize',20)
   xlim([ -0.1 0.25])
   vline(0.0,'k-')
   vline(-0.01,'k--')
   vline(0.01,'k--')
   box on
   saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'Object_size_uncertainty.fig'])
   
   unc_size.bias.x = mean(Quadric_axes(1,:)-0.05)+0.01; % add 10 micron standard deviation uncertainty
   unc_size.bias.y = mean(Quadric_axes(2,:)-0.05)+0.01; % add 10 micron standard deviation uncertainty
   unc_size.bias.z = mean(Quadric_axes(3,:)-0.05)+0.01; % add 10 micron standard deviation uncertainty
   unc_size.random.x = std(Quadric_axes(1,:)-0.05)+0.01; % add 10 micron standard deviation uncertainty
   unc_size.random.y = std(Quadric_axes(2,:)-0.05)+0.01; % add 10 micron standard deviation uncertainty
   unc_size.random.z = std(Quadric_axes(3,:)-0.05)+0.01; % add 10 micron standard deviation uncertainty
   
   
%    
%    figure(3)   
%    histogram((4*pi/3) .* (0.5.*Quadric_axes(1,:)) .* (0.5.*Quadric_axes(2,:)) .* (0.5.*Quadric_axes(3,:)),'FaceColor','b','EdgeColor','none')
%    xlim( [0 0.002] )
%    vline((4/3)*pi*0.025^3)
%    ylabel('[\#]','interpreter','latex','fontsize',20)
%    xlabel('Particle Volume ($mm^3$)','interpreter','latex','fontsize',20)
%    box on
%    saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'object_volume.fig'])
   
   figure(5)
   hold all
   histogram(Pobj(17,:),'FaceColor','k','EdgeColor','none')
   histogram(Accuracy,'FaceColor','b','EdgeColor','none')
   ylabel('[\#]','interpreter','latex','fontsize',20)
   xlabel('$\epsilon^*_{reproj}$','interpreter','latex','fontsize',20)
   legend('$\epsilon^*_{reproj}$ (raw)','$\epsilon^*_{reproj}$ (processed)','interpreter','latex','fontsize',20)
%    set(gca,'xscale','log')
   vline(exp(mean(log(Accuracy))),'b-')
   vline(exp(mean(log(Accuracy))+2*std(log(Accuracy))),'b--')
   vline(exp(mean(log(Accuracy))-2*std(log(Accuracy))),'b--')
   
   vline(exp(mean(log(Pobj(17,:)))),'k-')
   vline(exp(mean(log(Pobj(17,:)))+2*std(log(Pobj(17,:)))),'k--')
   vline(exp(mean(log(Pobj(17,:)))-2*std(log(Pobj(17,:)))),'k--')
   
   box on
   set(gca,'linewidth',1)
   saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'reproj_error.fig'])
   
   % save bias error and 2-sigma bounds of random error
   unc_reproj.bias = exp(mean(log(Accuracy)));
   unc_reproj.random.low = exp(mean(log(Accuracy))-2*std(log(Accuracy)));
   unc_reproj.random.high = exp(mean(log(Accuracy))+2*std(log(Accuracy)));
   
end


%% Displacement/Velocity field uncertainty estimate

% The uncertainty analysis done here substitutes the practice of
% mass-conservation law compliance that one can do with dense velocimetry.
% As our measurements are sparse, the uncertainty of the displacement is
% then given by three components: 1) uncertainty due to enlarged
% reconstructed size of the cavity; 2) uncertainty due to reprojection
% errors; and, 3) uncertainty due to residuals of the track fitting
% procedure.

% The displacement error is then given by following:
% unc_displacement 
%      = FitRes + Accuracy.*Quadric_axes + ((Quadric_axes) - 0.05 + 0.01);

unc_displacement = repmat(FitRes,3,1)  + ((Quadric_axes) - 0.05 + 0.01) + Accuracy.*Quadric_axes;

figure(51)

subplot(2,1,1)
hold on
histogram((unc_displacement(1,:)./(1/5000))./1000,'FaceColor','r','EdgeColor','none')
histogram((unc_displacement(2,:)./(1/5000))./1000,'FaceColor','b','EdgeColor','none')
histogram((unc_displacement(3,:)./(1/5000))./1000,'FaceColor','g','EdgeColor','none')
ylabel('[\#]','interpreter','latex','fontsize',20)
xlabel('$\epsilon_{U}$ (m/sec)','interpreter','latex','fontsize',20)
legend('$\epsilon_{Ux}$','$\epsilon_{Uy}$','$\epsilon_{Uz}$','interpreter','latex','fontsize',20)
xlim([0 1.5])
box on
set(gca,'linewidth',1)
hold off

subplot(2,1,2)
histogram((vecnorm(unc_displacement)./(1/5000))./1000,'FaceColor','r','EdgeColor','none')
ylabel('[\#]','interpreter','latex','fontsize',20)
xlabel('$\epsilon_{U}$ (m/sec)','interpreter','latex','fontsize',20)
xlim([0.5 2])
box on
set(gca,'linewidth',1)

saveas(gcf,[folder date rec vsl 'Postproc3DTracking' vsl 'velocity_error.fig'])
   
unc_velocity = ((unc_displacement)./(1/5000))./1000; % velocity uncertainty in m/sec for all three components

%% save data
save([folder date rec vsl 'Postproc3DTracking' vsl 'Uncertainties.mat'],'unc_fit','unc_size','unc_reproj','unc_velocity','outlier_stat')
save([folder date rec vsl 'Postproc3DTracking' vsl 'ParticleMetrics.mat'],'C_n','C_n_ghost')
save([folder date rec vsl 'Postproc3DTracking' vsl 'TrackQualityMetric.mat'],'Prob','epsilon','Prob_exp','rho_ghost','V_cav','t_spurious')


end