%% Load the data
% nrrd reader can be obtained at:
% https://www.mathworks.com/matlabcentral/fileexchange/66658-addisonelliott-matnrrd
addpath('nrrd_reader/')

% MRI scan data comes from the 2018 Atrial Segmentation Challenge:
[X,meta] = nrrdread('lgemri.nrrd');
[truth,truth_meta] = nrrdread('laendo.nrrd');
cavity = find(truth);

% Save the MRI 3D image and the ground truth labeled cavity.
save('data/MRI_X_cavity','X','cavity','truth')

%% Build the graph, may take a long time
% ei, ej stores indices for an edge (i,j)
% evi is the edge intensity score
% evd is the physical distance between nodes (i.e voxels) in the image

addpath('../algorithms/')

[ei,ej,evi,evd,N,Nmap] = graph_from_mri(X,1); 

% Save this for later. This will generate around 2.5 GB of data
save('data/MRI_Graph_Data.mat','ei','ej','evi','evd', 'N', 'Nmap', 'X', 'truth', '-v7.3');

% Nmap is a 3 x n array that maps (i,j,k) spatial indices in the MRI to
% a single node index in the graph. We save it to use for visualizations
% and figures later
save('data/MRI_Nmap','Nmap')

%% Form the adjacency matrix -- takes a long time

ev = exp(-evi.^2/(0.05)^2);
% This threshold show most of the differences with adjacent pixels
% and keeps their weights close.
minweight = 0.1;
ec = ev.*(ev >= minweight).*1/minweight;

% The weighted adjacency matrix and volume of the graph
mri_A = sparse(ei,ej,ec, numel(X),numel(X));
mri_vol = sum(nonzeros(mri_A));

% Save the graph, volume, and taret cavity. 
% This will be around 2.2 GB of data.
save('data/MRI_Graph','mri_A','mri_vol', 'cavity','-v7.3');


%% Get stats for the set of nodes associated with the target atrial cavity
[cut,vol,edges,conductance] = set_stats(mri_A,cavity,mri_vol)

%% Show some visualizations
addpath('../algorithms')
addpath('../algorithms/vol3d')

figure(1)
show_set(X,Nmap,cavity,cavity,2)
axis off

figure(2)
show_set_iso(X,Nmap,cavity,cavity,1)
axis off