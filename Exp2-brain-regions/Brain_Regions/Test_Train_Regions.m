%% Dividing the regions into training and test sets.

load Brain_All_Regions_Sorted


Under5000 = find(RegionSize < 5000);

Under7500 = find(RegionSize < 7500);
Under10000 = find(RegionSize < 10000);
Under15000 = find(RegionSize < 15000);
Under20000 = find(RegionSize < 20000);
Under30000 = find(RegionSize < 30000);

Set5000to7499 = setdiff(Under7500,Under5000);
Set7500to9999 = setdiff(Under10000,Under7500);
Set10000to14999 = setdiff(Under15000,Under10000);
Set15000to19999 = setdiff(Under20000,Under15000);
Set20000to30000 = setdiff(Under30000,Under20000);

G1 = Under5000';
G2 = Set5000to7499';
G3 = Set7500to9999';
G4 = Set10000to14999';
G5 = Set15000to19999';
G6 = Set20000to30000';

rng(1);
p1 = randperm(numel(G1));
p2 = randperm(numel(G2));
p3 = randperm(numel(G3));
p4 = randperm(numel(G4));
p5 = randperm(numel(G5));
p6 = randperm(numel(G6));

% Take a random few from each category
TrainInds = sort([G1(p1(1:3)), G2(p2(1:3)), G3(p3(1:3)), G4(p4(1:2)), G5(p5(1:2)), G6(p6(1:2)), 92, 94])';
TestInds = setdiff(1:95,TrainInds);

RegionsTrain = Regions(:,TrainInds);
RegionsTest = Regions(:,TestInds);
Seed1Train = SeedStart1(:,TrainInds);
Seed1Test = SeedStart1(:,TestInds);
Seed2Train = SeedStart2(:,TrainInds);
Seed2Test = SeedStart2(:,TestInds);
RegionCondTest = RegionCond(TestInds);
RegionCondTrain = RegionCond(TrainInds);
RegionSizeTest = RegionSize(TestInds);
RegionSizeTrain = RegionSize(TrainInds);
RegionLabelTest = RegionLabel(TestInds);
RegionLabelTrain = RegionLabel(TrainInds);

%save('BrainRegions_Train','RegionsTrain','Seed1Train','Seed2Train','RegionCondTrain','RegionLabelTrain','RegionSizeTrain')
%asave('BrainRegions_Test','RegionsTest','Seed1Test','Seed2Test','RegionCondTest','RegionLabelTest','RegionSizeTest')
