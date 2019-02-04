## README for snap experiments

This readme provides details for obtaining the data and reproducing the community detection experiments from the paper

"Flow-Based Local Graph Clustering with Better Seed Set Inclusion"  
Nate Veldt, Christine Klymko, and David Gleich  
*Proceedings of the 2019 SIAM International Conference on Data Mining*

ArXiv preprint: [https://arxiv.org/abs/1811.12280]()

## Datasets

See datasets folder. 

### Original Data
The datasets used for our experiments were downloaded from the SNAP group on the SuiteSparse matrix collection ([https://sparse.tamu.edu/]()). The are also available directly on the SNAP website ([https://snap.stanford.edu/data/index.html]()).

SuiteSparse provides the data in .mat files.

com-Amazon.mat  
com-LiveJournal.mat  
com-DBLP.mat  
com-Orkut.mat

Each .mat files comes with a set of top 5000 communities as identified by Yang and Leskovec.

J. Yang and J. Leskovec. Defining and Evaluating Network Communities based 
on Ground-truth. ICDM, 2012.  [http://arxiv.org/abs/1205.6233 ]() 

### Extracting Largest 10 communities

For our experiments, we download .mat files for Amazon, LiveJournal, DBLP, and Orkut. Once these are downloaded and stored, we extract the 10 largest communities in MATLAB with the files

Top10communities.m
Amazon_top10.m

We have a separate file for Amazon, because some of the largest communties are in fact identical to each other. We do a little extra work to ensure we have the top 10 unique and distinct largest communities for this dataset.

The file community_stats_average_table was used to compute statistics for the communities from each network, to put them in Table 1 in the full version of the paper on arXiv: [https://arxiv.org/abs/1811.12280]()

### Creating edgelists

The .mat files can be used in Julia and Matlab.

The LocalGraphClustering package (https://github.com/kfoynt/LocalGraphClustering) which we used for some of our experiments requires the graphs be store in edgelist format. The MATLAB file save_edgelists.m takes the graphs stored in .mat form and outputs a .edgelist file. This is slow and isn't meant to be a very efficient way to accomplish this. However, this only needs to be run once.


### Seed sets

For reproduceability, after generating random seeds for each of the communities, we store them. This was done by running theh MATLAB file Extract_SeedSets.m. The seed sets are stored in 

{graph name}-seed-starter.mat


These are not the seed sets themselves, but a random sample (2%, 3%, and 5%), or each target community which we then grow by their neighborhood in order to produce seed sets. Results in our paper are reported for using 5% of the target set.

## Experiments

Experiments reported in the paper can be reproduced by running the code in the following files:

`HK_PR_standard.m`  
`CRD_Experiments.m `  
`SimpleLocal_FlowSeed_Experiments.jl`

Output from our experiments is stored in the folder Output5.

Before running these experiments, you will need to ensure the datasets have been downloaded to the datasets/ folder in the correct format, as outline above.

## Output

In order to summarize the experiment data and output results to a latex table, run

`print_output_to_table.m`

This gives code for tables that are used in both the full version of the paper on arXiv, and the 9-page SDM paper.