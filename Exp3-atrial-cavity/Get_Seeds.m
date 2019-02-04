load('data/MRI_Graph.mat')

% Randomly select 100 nodes from the atrial cavity
nseeds = 100;
p = randperm(numel(cavity));
seed_start = cavity(p(1:nseeds));

% The starting set of 100 nodes used in experiments from 
% (Veldt, Klymko, Gleich 2019) are included in the repository.

%save('data/cavity_100_nodes.mat','seed_start')