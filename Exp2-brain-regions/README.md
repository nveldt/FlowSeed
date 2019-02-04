# README for Brain MRI experiments

This readme provides details for obtaining the data and reproducing the brain MRI region detection experiments from the paper

"Flow-Based Local Graph Clustering with Better Seed Set Inclusion"  
Nate Veldt, Christine Klymko, and David Gleich  
*Proceedings of the 2019 SIAM International Conference on Data Mining*

ArXiv preprint: [https://arxiv.org/abs/1811.12280](https://arxiv.org/abs/1811.12280)

Note that we use the same graph previously used in experiments in the paper

"A Simple and Strongly-Local Flow-Based Method for Cut Improvement" 
Nate Veldt, David Gleich, and Michael Mahoney
Proceedings of the 2016 International Conference on Machine Learning

## Original Data

The brain graph is constructed from MRI data we originally obtained as a part of the MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling

[http://www.miccai.org/news/miccai-2012-workshops-and-challenges-call-papers](http://www.miccai.org/news/miccai-2012-workshops-and-challenges-call-papers)
The original data came from the OASIS project: [https://www.oasis-brains.org/](https://www.oasis-brains.org/). Specifically, we used the OASIS-1 dataset:

* Marcus, Daniel S., et al. "Open Access Series of Imaging Studies (OASIS): cross-sectional MRI data in young, middle aged, nondemented, and demented older adults." Journal of cognitive neuroscience 19.9 (2007): 1498-1507.

As of 1-22-19, the form for requesting the MICCAI 2012 data can still be found here: [https://my.vanderbilt.edu/masi/workshops/](https://my.vanderbilt.edu/masi/workshops/)

The data from MICCAI 2012 challenge contains many training and testing images. The MRI we used to constuct the brain graph can from file `1000_3.nii`, and the corresonding region labels came from file `1000_3_glm.nii.gz`.

## Turning MRI into a Graph

In order to extract the MRI and region labels from the .nii files, we used the nii reader available at:

[https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image/content/load_nii.m](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image/content/load_nii.m)

We include a copy of this software in this GitHub repo.

Make sure that `1000_3.nii` and `1000_3_glm.nii.gz` have been copied to the "data" folder, or update the path to where they are saved.

Following the steps in `build_brain_graph.m` will allow you to obtain the adjacency matrix for the brain graph, and the set of indices for the target ventricle used in experiments in (Veldt, Gleich, Mahoney, 2016).

## Brain Regions and Seed Sets

The folder `Brain_Regions` stores all 95 regions of the brain, split into two sets: 17 example (or training) regions, and 68 evaluation (or testing) regions. We extracted seed sets and saved all .mat files using files
`Save_Regions.m`, `Regions_new_seed_sets_2_3_percent.m`, and `Test_Train_Regions.m`.

## Running Experiments

Experiments on example (training) and evaluation (testing) regions are performed in MATLAB for pprpush, and in Julia 0.6 for FlowSeed.

Experiments on training/example regions can be reproduced using `PPR_Train1.m` (ppr on first to types of seed sets), `PPR_Train2.m` (ppr with the largest two types of seed sets), and `FlowSeed_Train.jl`.

Experiments on the testing/evaluation regions are in `FlowSeed_Test.jl` and `PPR_Test.m`.

The file `Learn_Best_Parameters.jl` extracts the best locality parameter for FlowSeed, and the best PageRank tolerance parameter for PPR push, based on the 17 example regions. These are using in the experiments on the evaluation regions.

## Plots

The Plots folder contains several files for plotting results in Julia.
