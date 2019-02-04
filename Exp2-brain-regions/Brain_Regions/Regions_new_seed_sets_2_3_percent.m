% Extract regions to work with.
% This script also saves two new seed sets, with 2% and 3% of the target 
% sets. The names 'Seedset1' and 'Seedset2' are used already, so these 
% will be called Seedset2per and Seedset3per since they correspond to 2% 
% and 3% of thetarget region

%% Load the brain graph 
addpath('data')
addpath('../algorithms')
load local_flow_icml_data
load BrainGraph
n = size(Brain,1);

%% Save all regions

% max(max(max(L))) = 207, the highest labeled region

next = 0;
Regions = sparse(n,207);
SeedStart2per = sparse(n,207);
SeedStart3per = sparse(n,207);
RegionCond = zeros(207,1);
RegionSize = zeros(207,1);
RegionLabel = zeros(207,1);
rng(1); % reset random numbers for reproduceability

for region = 1:207
    
    FRind = L == region;
    FRset = find(FRind);
    if numel(FRset) >= 3000

        % Get largest connected component of the region
        C = B(FRset,FRset);
        [ci,comps] = graphconncomp(C);
        Rset = FRset(comps==mode(comps));

        % Get some stats for the region
        [cut,vol,edges,cond] = set_stats(B,Rset,volA);
        fprintf('Region %d has %d nodes, cond = %f \n',region,numel(Rset), cond)
 
        next = next+1;
        
        % Save stats for region
        RegionCond(next) = cond;
        RegionLabel(next) = region;
        RegionSize(next) = numel(Rset);
        
        % Get two types of seed sets: 100 nodes, and 1% of the region
        Target = Rset;
        p = randperm(numel(Target));
        nseeds = round(numel(Target)*.02);
        nseeds2 = round(numel(Target)*.03);
        S2per = Target(p(1:nseeds));
        S3per = Target(p(1:nseeds2));

        Regions(Rset,next) = 1;
        SeedStart2per(S2per,next) = 1;
        SeedStart3per(S3per,next) = 1;
    
    end
end

Regions = sparse(Regions(:,1:next));
SeedStart2per = sparse(SeedStart2per(:,1:next));
SeedStart3per = sparse(SeedStart3per(:,1:next));
RegionSize = RegionSize(1:next);
RegionCond = RegionCond(1:next);
RegionLabel = RegionLabel(1:next);

%% Sort them by size

[RegionSize,order] = sort(RegionSize);
RegionCond = RegionCond(order);
RegionLabel = RegionLabel(order);
Regions = Regions(:,order);
SeedStart2per = SeedStart2per(:,order);
SeedStart3per = SeedStart3per(:,order);

for i = 1:size(Regions,2)
    
    Target = find(Regions(:,i));
    
    S3 = find(SeedStart2per(:,i));
    S4 = find(SeedStart3per(:,i));
    
    assert(all(ismember(S3,Target)))
    assert(all(ismember(S4,Target)))
    
end

% Running the above produces the seed sets
%save('Brain_All_Regions_Seedsets_2_3_percent','Regions','RegionLabel','RegionCond','RegionSize','SeedStart2per','SeedStart3per')