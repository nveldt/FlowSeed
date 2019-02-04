%% Load the data
% nii readers can be obtained at: 
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image/content/load_nii.m
addpath('nii_reader/')

% MRI scan data comes from MICCAI 2012 challenge:
% https://my.vanderbilt.edu/masi/workshops/
nii = load_nii('data/1000_3.nii');
niiL = load_nii('data/1000_3_glm.nii');

%% Save the MRI scan and the labels as X and L

X = nii.img;
L = niiL.img;

save('data/MRI_Labels','X','L')

%% Build the graph, may take a long time
addpath('../algorithms/')
[ei,ej,evi,evd,N,Nmap] = graph_from_mri(X,1);

save('data/local_flow_icml_data.mat',...
    'ei','ej','evi','evd', 'N', 'Nmap', 'X', 'L', '-v7.3');

%% Form the adjacency matrix -- takes a long time
ev = exp(-evi.^2/(0.05)^2);
% This threshold show most of the differences with adjacent pixels
% and keeps their weights close.
minweight = 0.1;
% 
ec = ev.*(ev >= minweight).*1/minweight;
B = sparse(ei,ej,ec, numel(X),numel(X));
volA = sum(nonzeros(B));

%% Get an example region
% Region 52 is the one used in (Veldt, Gleich, Mahoney, ICML 2016)
addpath('../algorithms/')

region = 52;
FRind = L == region;
[cut,vol,edges,cond] = set_stats(B,find(FRind),volA)
FRset = find(FRind);

%% Get largest connected component.
% There is a tiny secondary component that is significantly separated from 
% the rest of the ventricle. 

Cset = B(FRset,FRset);

[~,C] = graphconncomp(Cset,'DIRECTED',false,'Weak',true);
edges = mode(C);

Rset = FRset(find(C == edges));

save('data/BrainGraph','B','-v7.3')