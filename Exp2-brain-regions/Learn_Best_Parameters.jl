# Determine the best parameters to set for FlowSeed and PPR Push

BestParams = zeros(Int64,4,17)
include("Plots/Plot_Helper.jl")

epsis = [0.05; .1; .25; .5]
Sizes = []
# There are 4 seed types:
# 1: 100 nodes from target region
# 2: 1% of target region
# 3: 2% of target region
# 4: 3% of target region
for seedtype = 1:4

    ## Get the F1 scores for different runs of FlowSeed, with both soft and strict penalties
    p1, r1, f1, Sizes, seednodes =  LoadFlowSeedTrain(.05,seedtype)
    p2, r2, f2, Sizes, seednodes =  LoadFlowSeedTrain(.1,seedtype)
    p3, r3, f3, Sizes, seednodes =  LoadFlowSeedTrain(.25,seedtype)
    p4, r4, f4, Sizes, seednodes =  LoadFlowSeedTrain(.5,seedtype)

    FS_f1 = [f1[4,:]'; f2[4,:]'; f3[4,:]'; f4[4,:]']

    # For each region, get the best F1 score, and which tolerance lead to it
    for i = 1:17
        f1b, loc = findmax(FS_f1[:,i],1)
        BestParams[seedtype,i] = loc[1]
    end

end
## Get the best of PPR, with various tolerances
BestParamsPPR = zeros(Int64,4,17)

alpha = .6
for seedtype = 1:4

    prPPR, rePPR, f1PPR, Sizes = LoadPPRTrain(alpha,seedtype)
    PR = f1PPR
    for i = 1:17
        f1b, loc = findmax(PR[:,i],1)
        BestParamsPPR[seedtype,i] = loc[1]
    end

end

Tols = [1e-7; 1e-8; 1e-9; 1e-10; 1e-11]
epsis = [0.05; .1; .25; .5]

using MAT

#matwrite("Training17_Results.mat", Dict( "Tols" => Tols, "epsis" => epsis,
#"BestParamsPPR" => BestParamsPPR, "BestParamsFlowSeed" => BestParams, "Sizes" => Sizes))
