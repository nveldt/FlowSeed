# FlowSeed

Code for algorithms and experiments included in 

"Flow-Based Local Graph Clustering with Better Seed Set Inclusion"  
Nate Veldt, Christine Klymko, and David Gleich  
*Proceedings of the 2019 SIAM International Conference on Data Mining*

ArXiv preprint: [https://arxiv.org/abs/1811.12280]()

## Datasets

We ran experiments on graphs with known communities from the SNAP repository. We also performed experiments on graphs constructed from MRI scans. The first is an MRI of the brain, provided as a part of the MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling ([https://my.vanderbilt.edu/masi/workshops/]()). The second graph was constructed from a full-body MRI scan provided as a part of the 2018 Atrial Segmentation Challenge ([http://atriaseg2018.cardiacatlas.org/]()).

This repository contains a subfolder for each experiment (snap experiments, brain experiment, atrial cavity experiments). Further details for each dataset and each experiment are contained in the appropxiate folder.


## Dependencies and Outside Code

For SimpleLocal (Veldt, Gleich, and Mahoney 2016) and Capacity Releasing Diffusion, we use the LocalGraphClustering package: 

[https://github.com/kfoynt/LocalGraphClustering]()

For PPR Push (Andersen, Chung, Lang 2006), we used a C++ implementation written by David Gleich, with a MATLAB front end. The implementation is included in the algorithms folder.

For HK-relax (Kloster and Gleich, 2014), we used the C++ implementation (with MATLAB interface) provided at

[https://github.com/kkloste/hkgrow]()

## Note on Julia Versions

The FlowSeed algorithm is implemented in Julia 0.6 and Julia 1.0.

The experiments with FlowSeed, SimpleLocal, and Capacity Releasing Diffusion were originally run using Julia 0.6 as a front end. Therefore, all code for these experiments is written for this version of Julia. Slight edits can be made to run a similar set of experiments in Julia 1.0 if desired.


