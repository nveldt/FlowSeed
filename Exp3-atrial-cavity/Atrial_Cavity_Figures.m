% Load the MRI scan and the ground truth segmengation
addpath('nrrd_reader/')
[X,meta] = nrrdread('lgemri.nrrd');
[truth,truth_meta] = nrrdread('laendo.nrrd');

%% Load and visualize the target region

load data/MRI_Nmap
TargetRegion = find(truth);

%% Seed set
addpath('../algorithms')
load data/Seed100plusNeighbs
figure(10)
show_set_iso(X,Nmap,R,TargetRegion,1)
axis off
axis square
%print(gcf,sprintf('seed_heart.png'),'-dpng');

%% FlowSeed figures
load FlowSeed_refine_PPR_alpha_0.6_epsilon_0.1.mat
for t = 1:3
    fs = find(sets2(:,t));

    figure(t)
    show_set_iso(X,Nmap,fs,TargetRegion,0)
    axis off
    axis square
    print(gcf,sprintf(strcat('Figures/FlowSeed_heart',num2str(t),'refine_alpha60.png')),'-dpng');
end

%% PPR Figures

load PPR_Output/alpha_60sets_stats.mat
for t = 1:7
    ppr = find(sets(:,t));
    figure(t)
    show_set_iso(X,Nmap,ppr,TargetRegion,0)
    axis off
    axis square
    print(gcf,sprintf(strcat('Figures/PPR_heart',num2str(t+7),'alpha60.png')),'-dpng');
end

%% Show just the ground truth cavity

figure(10)
show_set_iso(X,Nmap,TargetRegion,[],1)
axis square
print(gcf,sprintf('Figures/left_atrial_cavity.png'),'-dpng');
