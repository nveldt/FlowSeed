% Running PPR push on the evaluation ("Test") regions

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
addpath('New_Region_Data/')
load BrainRegions_Test23
Regions = RegionsTest;
RegionSize = RegionSizeTest;
RegionCond = RegionCondTest;
RegionLabels = RegionLabelTest;

%% Load the training size-to-tolerance sets

load Training17_Results

%%
alpha = .6;

for seedtype = 1:4
    
    % Select the type of seed set
    if seedtype == 1
        Seeds = Seed1Test;
    elseif seedtype == 2
        Seeds = Seed2Test;
    elseif seedtype == 3
        Seeds = Seed2perTest;
    else
        Seeds = Seed3perTest;
    end
    
    
        outputfile = strcat('PPR_Test_Output/Seed_',num2str(seedtype),'_alpha_',num2str(alpha*100),'.txt');
        fid = fopen(outputfile,'w');
        
        % For each file data (fixed alpha, fixed set of seeds), we will store the output
        sets= sparse(n,78);
        stats = zeros(7,78);

        
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
           [~,volR,edges,condR] = set_stats(G,R,volA);
           fprintf('\tSeed set:conductance = %f, size = %d \n',condR,numel(R))
           fprintf('%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');

          fprintf(fid,'\tSeed set:conductance = %f, size = %d \n',condR,numel(R));
          fprintf(fid,'%3s\t%6s\t%8s\t %8s \t %8s \t %8s \t %8s \n', 'Tol','Size', 'Time','Cond', 'precision','recall','F1-score');
                             
     
           % Determine the best tolerance to use
           [~,closest_region] = min(abs(Sizes-Size));
           best_tol = BestParamsPPR(seedtype,closest_region);
           
           fprintf('Test region is closest in size to a training region with %d nodes \n',Sizes(closest_region))
           fprintf(fid, 'Test region is closest in size to a training region with %d nodes \n',Sizes(closest_region));

           % Use whichever tolerence worked best for that training region
           tol = Tols(best_tol);
           
           fprintf('Using tolerance = %f \n',tol)
           fprintf(fid,'Using tolerance = %f \n',tol);

           tic
           [PR,cond,cut,vol,y] = pprpush_weighted_mex(G, d, R, tol, alpha,-1-1.e-4);
           tPR = toc;
           [pPR,rPR,fPR] = AdjustedPRF(Target,PR,[]);
           fprintf('%d \t %d \t %f \t %f \t %f \t %f \t %f \n',round(log10(tol)),numel(PR),tPR,cond,pPR,rPR,fPR)
           fprintf(fid,'%d \t %d \t %f \t %f \t %f \t %f \t %f \n',round(log10(tol)),numel(PR),tPR,cond,pPR,rPR,fPR);

           sets(PR,regi) = 1;
           stats(:,regi) = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];

          
        end
        fclose(fid);
        save(strcat('PPR_Test_Output/Seed_',num2str(seedtype),'_alpha_',num2str(alpha*100),'sets_stats'),'sets','stats')
        
end