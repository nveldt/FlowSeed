% Extract new regions to work with
addpath('data')
addpath('../algorithms')

%% Load the graph
load local_flow_icml_data
load BrainGraph
n = size(A,1);

%% Save all regions

% max(max(max(L))) = 207, the highest labeled region
next = 0;
Regions = sparse(n,207);
SeedStart1 = sparse(n,207);
SeedStart2 = sparse(n,207);
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
        nseeds = 100;
        nseeds2 = round(numel(Target)*.01);
        S1 = Target(p(1:nseeds));
        S2 = Target(p(1:nseeds2));

        Regions(Rset,next) = 1;
        SeedStart1(S1,next) = 1;
        SeedStart2(S2,next) = 1;
    
    end
end

Regions = sparse(Regions(:,1:next));
SeedStart1 = sparse(SeedStart1(:,1:next));
SeedStart2 = sparse(SeedStart2(:,1:next));
RegionSize = RegionSize(1:next);
RegionCond = RegionCond(1:next);
RegionLabel = RegionLabel(1:next);

%% Sort them by size

[RegionSize,order] = sort(RegionSize);
RegionCond = RegionCond(order);
RegionLabel = RegionLabel(order);
Regions = Regions(:,order);
SeedStart1 = SeedStart1(:,order);
SeedStart2 = SeedStart2(:,order);

% Running the above produces the seed sets:
% save('Brain_All_Regions_Sorted','Regions','RegionLabel','RegionCond','RegionSize','SeedStart1','SeedStart2')