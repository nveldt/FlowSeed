% Running PPR push experiments with the two smallest seed set sizes

addpath('../algorithms')
addpath('../algorithms/pprpush/')
addpath('data/')
addpath('Brain_Regions')

%%
load BrainGraph
G = B;
volA = sum(nonzeros(B));
d = sum(G,1)';
n = size(G,1);
%% Load the brain graph dataset
addpath('Brain_Regions/')
load BrainRegions_Train
Regions = RegionsTrain;
RegionSize = RegionSizeTrain;
RegionCond = RegionCondTrain;
RegionLabels = RegionLabelTrain;

%%
alphas = [.5,.6,.7,.8,.9];
tols = [1e-7,1e-8,1e-9,1e-10,1e-11];

for seedtype = 1:2
    
    % Select the type of seed set
    if seedtype == 1
        Seeds = Seed1Train;
    else
        Seeds = Seed2Train;
    end
    
    
    for alpha = alphas
        outputfile = strcat('PPR_Train_Output/Seed_',num2str(seedtype),'_alpha_',num2str(alpha*100),'.txt');
        fid = fopen(outputfile,'w');
        
        % For each file data (fixed alpha, fixed set of seeds), we will store the output
        sets_tol_7 = sparse(n,17);
        sets_tol_8 = sparse(n,17);
        sets_tol_9 = sparse(n,17);
        sets_tol_10 = sparse(n,17);
        sets_tol_11 = sparse(n,17);
        
        stats_tol_7 = zeros(7,17);
        stats_tol_8 = zeros(7,17);
        stats_tol_9 = zeros(7,17);
        stats_tol_10 = zeros(7,17);
        stats_tol_11 = zeros(7,17);
        
        
        for regi = 1:size(Regions,2)
            
           Target = find(Regions(:,regi));
           Label = RegionLabels(regi);
           Cond = RegionCond(regi);
           Size = RegionSize(regi);
           fprintf('\nRegion %d has conductance %f and size %d \n',Label,Cond,Size)
           fprintf(fid,'\nRegion %d has conductance %f and size %d \n',Label,Cond,Size);
            
           seed = find(Seeds(:,regi));
           Rneighbs = G(seed,:);
           [~,R] = find(Rneighbs);
           R = unique([R; seed]);
           [cut,volR,edges,condR] = set_stats(G,R,volA);
           fprintf('\tSeed set:conductance = %f, size = %d \n',condR,numel(R))
           fprintf('%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');

          fprintf(fid,'\tSeed set:conductance = %f, size = %d \n',condR,numel(R));
          fprintf(fid,'%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');
                             
           for t = 1:numel(tols)
               tol = tols(t);
               tic
               [PR,cond,cut,vol,y] = pprpush_weighted_mex(G, d, R, tol, alpha,-1-1.e-4);
               tPR = toc;
               [pPR,rPR,fPR] = AdjustedPRF(Target,PR,[]);
               fprintf('%d \t %d \t %f \t %f \t %f \t %f \t %f \n',round(log10(tol)),numel(PR),tPR,cond,pPR,rPR,fPR)
               fprintf(fid,'%d \t %d \t %f \t %f \t %f \t %f \t %f \n',round(log10(tol)),numel(PR),tPR,cond,pPR,rPR,fPR);
               
               switch t
                   case 1
                       sets_tol_7(PR,regi) = 1;
                       stats_tol_7(:,regi) = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];
                   case 2
                       sets_tol_8(PR,regi) = 1;
                       stats_tol_8(:,regi) = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];
                   case 3
                       sets_tol_9(PR,regi) = 1;
                       stats_tol_9(:,regi) = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];
                   case 4
                       sets_tol_10(PR,regi) = 1;
                       stats_tol_10(:,regi) = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];  
                   otherwise
                       sets_tol_11(PR,regi) = 1;
                       stats_tol_11(:,regi) = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];  
               end
               
           end
          
        end
        fclose(fid);
        save(strcat('PPR_Train_Output/Seed_',num2str(seedtype),'_alpha_',num2str(alpha*100),'sets_stats.mat'),'sets_tol_7','stats_tol_7','sets_tol_8','stats_tol_8','sets_tol_9','stats_tol_9','sets_tol_10','stats_tol_10')
    end
end