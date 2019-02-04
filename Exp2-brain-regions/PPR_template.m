% This file gives an example for how to run the PPR push method

addpath('../algorithms/pprpush/')
addpath('data/')
addpath('Brain_Regions')

%%
load BrainGraph
G = B;
volA = sum(nonzeros(B));
d = sum(G,1)';

%% Load the brain graph dataset
load BrainRegions_Train
Regions = RegionsTrain;
RegionSize = RegionSizeTrain;
RegionCond = RegionCondTrain;
RegionLabels = RegionLabelTrain;


%% Select a region

regi = 1;

Target = find(Regions(:,regi));

R1 = find(Seed1Train(:,regi));
R2 = find(Seed2Train(:,regi));

Rneighb = neighborhood(G,R1,1);


R =  unique([R1;Rneighb]);

[~,~,~,condR] = set_stats(G,R);

fprintf('Seed set has %d nodes and conductance %f \n',numel(R),condR);


%% Run ppr
tol = 1e-8;
alpha = .99;
tic
[PR,cond,cut,vol,y] = pprpush_weighted_mex(G, d, R, tol, alpha,-1-1.e-4);
tPR = toc;

[pPR,rPR,fPR] = AdjustedPRF(Target,PR,R);

fprintf('PPR push: Time = %f, Size %d, Conductance %f, PR = %f, Re = %f, F1 = %f \n',tPR,numel(PR),cond,pPR,rPR,fPR)
