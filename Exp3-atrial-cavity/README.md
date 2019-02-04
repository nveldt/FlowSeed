## README for Atrial Cavity Experiments

This readme provides details for obtaining the data and reproducing the atrial cavity detection experiments from the paper

"Flow-Based Local Graph Clustering with Better Seed Set Inclusion"  
Nate Veldt, Christine Klymko, and David Gleich  
*Proceedings of the 2019 SIAM International Conference on Data Mining*

ArXiv preprint: [https://arxiv.org/abs/1811.12280]()

## Original Data

The data was downloaded as a part of the 2018 Atrial Segmentation Challenge:

[http://atriaseg2018.cardiacatlas.org/]()

The training data consists of 100 full body MRI scans and accompanying ground truth segmentations of the left atrial cavity. The MRI scan used in our experiments was the training scan labeled `0RZDK210BSMWAA6467LU`. The file `lgemri.nrrd` stores the full body MRI, and the file `laendo.nrrd` stores the ground truth segmentation of the atrial cavity. 

## Turning MRI into a Graph

To read the files, you need an nrrd reader. You can download a reader from here. We include a copy of the following reader:

[https://www.mathworks.com/matlabcentral/fileexchange/66658-addisonelliott-matnrrd]()

(available also at [https://github.com/addisonElliott/matnrrd]() )

## Visualizing sets in the MRI

One can visualize the atrial cavity using the vol3d Matlab package:

[https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/6110/versions/3/previews/toolbox_fast_marching/vol3d.m/index.html]()

along with the `show_set.m` and `show_set_iso.m` functions (written by David Gleich). An example is included at the end of `build_fullbodymri_graph.m`.

## Seed sets

See `Get_seeds.m`. The 100 "starter" seed nodes from the atrial cavity that were used in the original paper are contained in the data folder. We grow these by their two-hop neighborhood to form the seed set we used for PPR-Push. 

The exact 100 random target nodes we used and the Push algorithm seed set are stored in the `data` folder.

## Experiments

See code in:

`PPR_atrial_cavity.jl`  
`FlowSeed_Refine_Step.jl`  

## More Visualizations and Figures

`Heart_Figures.m` shows how to generate the figures from the paper. These again require the data from the 2018 Atrial Segmentation Challenge. 

We include a couple extra files:

`MRI_Nmap.mat`  
`Seed100plusNeighbs.mat`

needed for visualizing some sets of nodes.