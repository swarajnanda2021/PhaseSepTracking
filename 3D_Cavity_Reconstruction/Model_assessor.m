%% Script for a thorough assessment of Deep-Learning models

% Assesses all models in Deep_Learning_Model_Database

clear all
close all
clc

%% Script for a thorough assessment of Deep-Learning models

folder='C:\Work\PostGraduation\Scripts'; % Desktop Location
date='\ParticleTracking'; % Date Folder
models='\Deep_Learning_Model_Database';


reordering_array = [1 8 6 7 4 5 2 3];

model_list = dir([folder date models '\*.mat']);
nocav=0;

model_list = model_list(reordering_array);


% Make table for book-keeping
name = {'thick_normal' ;'thick_even_less'; ...
    'thin_even_less' ;'thick_less' ;'thin_less' ;'noempty_thick' ;'noempty_thin' ;'thin_normal'};

name = name(reordering_array);

trainingsetsize = [1000; 488; 167; 640; 319; 808; 488; 616];

trainingsetsize = trainingsetsize(reordering_array);




for i=2%1:length(model_list)
    i
    net = load([model_list(i).folder vsl model_list(i).name]);
    
    disp(['evaluating ' model_list(i).name])
    
    if ismember(i,[1 3 5 7]) % thick cavity cases
        
        testset_all = 1001:1600;
        testset_cav = [1017:1080 1097:1160 1177:1240 1257:1320 1337:1400 1417:1480 1497:1560 1577:1600];
        testset_nocav = testset_all(~ismember(testset_all,testset_cav));
        
    else
        
        testset_all = 617:880;
        testset_cav = [617:680 697:760 777:840 857:880];
        testset_nocav = testset_all(~ismember(testset_all,testset_cav));
                
    end
    
    % evaluate cav cases
    counter=1;
    for ii=1:length(testset_cav)
        % Load test data
        testim_normal = imread([net.testImagesDir '\image_' num2str(testset_cav(ii)) '.tiff']);
        testim_noise1 = imnoise(testim_normal,'gaussian',0.01);
        testim_noise2 = imnoise(testim_normal,'salt & pepper',0.01);

        % Infer normal
        C1 = semanticseg(abs(testim_normal), net.net);
        B1 = labeloverlay(testim_normal,C1);
        
        % Infer gaussian noise
        C2 = semanticseg(abs(testim_noise1), net.net);
        B2 = labeloverlay(testim_noise1,C2);

        % Infer salt_pepper noise
        C3 = semanticseg(abs(testim_noise2), net.net);
        B3 = labeloverlay(testim_noise2,C3);


        labelim = imread([net.testLabelsDir '\labelled_image_' num2str(testset_cav(ii)) '.tiff']);

        figure(69)
        subplot(3,1,1)
        imagesc(B1)
        subplot(3,1,2)
        imagesc(B2)
        subplot(3,1,3)
        imagesc(B3)
        drawnow 
        pause
        
        % Gather similarity scores
        similarity_bfscore(1,counter) = bfscore(double(C1)-1, double(labelim));
        similarity_jaccardindex(1,counter) = jaccard(double(C1)-1, double(labelim));
        
        similarity_bfscore(2,counter) = bfscore(double(C2)-1, double(labelim));
        similarity_jaccardindex(2,counter) = jaccard(double(C2)-1, double(labelim));
        
        similarity_bfscore(3,counter) = bfscore(double(C3)-1, double(labelim));
        similarity_jaccardindex(3,counter) = jaccard(double(C3)-1, double(labelim));
        
        
        counter = counter+1;
    end
    
    figure(1)
    subplot(4,2,i)
    hold all
    
    [values, edges] = histcounts(similarity_jaccardindex(1,:),linspace(0.8,1,20));
    centers = (edges(1:end-1)+edges(2:end))/2;
    plot(centers, values, 'r-')
    
    [values, edges] = histcounts(similarity_jaccardindex(2,:),linspace(0.8,1,20));
    centers = (edges(1:end-1)+edges(2:end))/2;
    plot(centers, values, 'b-')

    [values, edges] = histcounts(similarity_jaccardindex(3,:),linspace(0.8,1,20));
    centers = (edges(1:end-1)+edges(2:end))/2;
    plot(centers, values, 'k-')
    
%     histogram(similarity_jaccardindex(1,:),'FaceAlpha',0.2,'DisplayStyle','stairs','EdgeColor',[0 0.4470 0.7410])
%     histogram(similarity_jaccardindex(2,:),'FaceAlpha',0.2,'DisplayStyle','stairs','EdgeColor',[0.4940 0.1840 0.5560])
%     histogram(similarity_jaccardindex(3,:),'FaceAlpha',0.2,'DisplayStyle','stairs','EdgeColor',[0.6350 0.0780 0.1840])
    xlim([0.8 1])
    xlabel('mIOU','interpreter','latex','FontSize',12)
    ylabel('[#]','interpreter','latex','FontSize',12)
    box on
    drawnow
    
    
    figure(2)
    subplot(4,2,i)
    hold all
    
    [values, edges] = histcounts(similarity_bfscore(1,:),linspace(0.5,1,20));
    centers = (edges(1:end-1)+edges(2:end))/2;
    plot(centers, values, 'r-')
    
    [values, edges] = histcounts(similarity_bfscore(2,:),linspace(0.5,1,20));
    centers = (edges(1:end-1)+edges(2:end))/2;
    plot(centers, values, 'b-')

    [values, edges] = histcounts(similarity_bfscore(3,:),linspace(0.5,1,20));
    centers = (edges(1:end-1)+edges(2:end))/2;
    plot(centers, values, 'k-')
    
%     histogram(similarity_bfscore(1,:),'FaceAlpha',0.2,'DisplayStyle','stairs','EdgeColor',[0 0.4470 0.7410])
%     histogram(similarity_bfscore(2,:),'FaceAlpha',0.2,'DisplayStyle','stairs','EdgeColor',[0.4940 0.1840 0.5560])
%     histogram(similarity_bfscore(3,:),'FaceAlpha',0.2,'DisplayStyle','stairs','EdgeColor',[0.6350 0.0780 0.1840])
    xlim([0.5 1])
    xlabel('bf-score','interpreter','latex','FontSize',12)
    ylabel('[#]','interpreter','latex','FontSize',12)
    box on
    drawnow
    
    similarity_bfscore_cav_mean(i,:) = mean(similarity_bfscore');
    similarity_jaccard_cav_mean(i,:) = mean(similarity_jaccardindex');
    
    similarity_bfscore_cav_std(i,:) = std(similarity_bfscore');
    similarity_jaccard_cav_std(i,:) = std(similarity_jaccardindex');
    
    
    similarity_bfscore = [];
    similarity_jaccardindex = [];
    
    if nocav==1
        similarity_bfscore = [];
        similarity_jaccardindex = [];

        counter=1;
        % evaluate nocav cases
        for ii=1:length(testset_nocav)
            % Load test data
            testim_normal = imread([net.testImagesDir '\image_' num2str(testset_nocav(ii)) '.tiff']);
            testim_noise1 = imnoise(testim_normal,'gaussian',0.01);
            testim_noise2 = imnoise(testim_normal,'salt & pepper',0.01);

            % Infer normal
            C1 = semanticseg(abs(testim_normal), net.net);
            B1 = labeloverlay(testim_normal,C1);

            % Infer gaussian noise
            C2 = semanticseg(abs(testim_noise1), net.net);
            B2 = labeloverlay(testim_noise1,C2);

            % Infer salt_pepper noise
            C3 = semanticseg(abs(testim_noise2), net.net);
            B3 = labeloverlay(testim_noise2,C3);


            labelim = imread([net.testLabelsDir '\labelled_image_' num2str(testset_nocav(ii)) '.tiff']);


            % Gather similarity scores
            similarity_bfscore(1,counter) = bfscore(double(C1)-1, double(labelim));
            similarity_jaccardindex(1,counter) = jaccard(double(C1)-1, double(labelim));

            similarity_bfscore(2,counter) = bfscore(double(C2)-1, double(labelim));
            similarity_jaccardindex(2,counter) = jaccard(double(C2)-1, double(labelim));

            similarity_bfscore(3,counter) = bfscore(double(C3)-1, double(labelim));
            similarity_jaccardindex(3,counter) = jaccard(double(C3)-1, double(labelim));


            counter = counter+1;
        end

        similarity_bfscore_nocav_mean(i,:) = mean(similarity_bfscore');
        similarity_jaccard_nocav_mean(i,:) = mean(similarity_jaccardindex');

        similarity_bfscore_nocav_std(i,:) = std(similarity_bfscore');
        similarity_jaccard_nocav_std(i,:) = std(similarity_jaccardindex');
    end
    
    
    
    
    % Calculate Epoch averaged training loss (dice coefficient)
    
    arr = net.info.TrainingLoss;
    n = length(net.info.TrainingLoss)./50;
    
    trainingloss_epoch_avgd(i,:) = mean(reshape(arr(1:n * floor(numel(arr) / n)), [], n), 2);
    
    
    

end


%%

figure(69420)
subplot(1,2,1)
plot(trainingloss_epoch_avgd([1 3 5 7],:)')
xlabel('[\#] Epochs','interpreter','latex','FontSize',14)
ylabel('Dice loss (epoch avg.)','interpreter','latex','FontSize',14)
legend('1000','808','640','488','interpreter','latex','FontSize',10)
% grid on

subplot(1,2,2)
plot(trainingloss_epoch_avgd([2 4 6 8],:)')
xlabel('[\#] Epochs','interpreter','latex','FontSize',14)
ylabel('Dice loss (epoch avg.)','interpreter','latex','FontSize',14)
legend('616','488','319','167','interpreter','latex','FontSize',10)
% grid on


%%

if nocav==0

    % cav
    mIOU_normal_cav_mean = similarity_jaccard_cav_mean(:,1);
    mIOU_normal_cav_std = similarity_jaccard_cav_std(:,1);

    bfs_normal_cav_mean = similarity_bfscore_cav_mean(:,1);
    bfs_normal_cav_std = similarity_bfscore_cav_std(:,1);

    mIOU_gaussian_cav_mean = similarity_jaccard_cav_mean(:,2);
    mIOU_gaussian_cav_std = similarity_jaccard_cav_std(:,2);

    bfs_gaussian_cav_mean = similarity_bfscore_cav_mean(:,2);
    bfs_gaussian_cav_std = similarity_bfscore_cav_std(:,2);


    mIOU_saltpep_cav_mean = similarity_jaccard_cav_mean(:,3);
    mIOU_saltpep_cav_std = similarity_jaccard_cav_std(:,3);

    bfs_saltpep_cav_mean = similarity_bfscore_cav_mean(:,3);
    bfs_saltpep_cav_std = similarity_bfscore_cav_std(:,3);
    
     % Consolidate into tables
    % cav
    T1 = table(name,trainingsetsize,mIOU_normal_cav_mean,mIOU_normal_cav_std...
        ,bfs_normal_cav_mean,bfs_normal_cav_std,mIOU_gaussian_cav_mean...
        ,mIOU_gaussian_cav_std,bfs_gaussian_cav_mean,bfs_gaussian_cav_std...
        ,mIOU_saltpep_cav_mean,mIOU_saltpep_cav_std,bfs_saltpep_cav_mean...
        ,bfs_saltpep_cav_std);

    T1_thick = T1([1,3,5,7],:);
    T1_thin = T1([2,4,6,8],:);


    save('summary_Deep_model_assessment.mat','T1_thin','T1_thick','trainingloss_epoch_avgd')
end
% nocav
if nocav==1

    mIOU_normal_nocav_mean = similarity_jaccard_nocav_mean(:,1);
    mIOU_normal_nocav_std = similarity_jaccard_nocav_std(:,1);

    bfs_normal_nocav_mean = similarity_bfscore_nocav_mean(:,1);
    bfs_normal_nocav_std = similarity_bfscore_nocav_std(:,1);

    mIOU_gaussian_nocav_mean = similarity_jaccard_nocav_mean(:,2);
    mIOU_gaussian_nocav_std = similarity_jaccard_nocav_std(:,2);

    bfs_gaussian_nocav_mean = similarity_bfscore_nocav_mean(:,2);
    bfs_gaussian_nocav_std = similarity_bfscore_nocav_std(:,2);


    mIOU_saltpep_nocav_mean = similarity_jaccard_nocav_mean(:,3);
    mIOU_saltpep_nocav_std = similarity_jaccard_nocav_std(:,3);

    bfs_saltpep_nocav_mean = similarity_bfscore_cav_mean(:,3);
    bfs_saltpep_nocav_std = similarity_bfscore_cav_std(:,3);

    % Consolidate into tables
    % cav
    T1 = table(name,trainingsetsize,mIOU_normal_cav_mean,mIOU_normal_cav_std...
        ,bfs_normal_cav_mean,bfs_normal_cav_std,mIOU_gaussian_cav_mean...
        ,mIOU_gaussian_cav_std,bfs_gaussian_cav_mean,bfs_gaussian_cav_std...
        ,mIOU_saltpep_cav_mean,mIOU_saltpep_cav_std,bfs_saltpep_cav_mean...
        ,bfs_saltpep_cav_std);

    T1_thick = T1([1,3,5,7],:);
    T1_thin = T1([2,4,6,8],:);

    % nocav
    T2 = table(name,trainingsetsize,mIOU_normal_nocav_mean,mIOU_normal_nocav_std...
        ,bfs_normal_nocav_mean,bfs_normal_nocav_std,mIOU_gaussian_nocav_mean...
        ,mIOU_gaussian_nocav_std,bfs_gaussian_nocav_mean,bfs_gaussian_nocav_std...
        ,mIOU_saltpep_nocav_mean,mIOU_saltpep_nocav_std,bfs_saltpep_nocav_mean...
        ,bfs_saltpep_nocav_std);


    T2_thick = T2([1,3,5,7],:);
    T2_thin = T2([2,4,6,8],:);


    save('summary_Deep_model_assessment.mat','T1_thin','T2_thin','T1_thick','T2_thick','trainingloss_epoch_avgd')

end


