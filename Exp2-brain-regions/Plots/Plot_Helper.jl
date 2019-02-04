# Load the training data for the specified alpha and seed type

    using MAT

function LoadFlowSeedTrain(epsilon,seedtype)
    if seedtype == 1
        seedtag = "1"
        seednodes = " 100 Target Nodes Given"
    elseif seedtype == 2
        seedtag = "2"
        seednodes = " 1% of Target Given"
    elseif seedtype == 3
        seedtag = "2per"
        seednodes = " 2% of Target Given"
    elseif seedtype == 4
        seedtag = "3per"
        seednodes = " 3% of Target Given"
    else
        println("Give a valid seedtype, 1 through 4")
        return
    end

    loadstring = "../FlowSeed_Train_Output/Seed_"*seedtag*"_epsilon_"*string(epsilon)*".mat"

    mat = matread(loadstring)

    stats1 = mat["stats1"]
    stats2 = mat["stats2"]
    stats3 = mat["stats3"]
    stats4 = mat["stats4"]


    # Load the Region sizes
    mat = matread("../Brain_Regions/BrainRegions_Train23.mat")
    Sizes = vec(round.(Int64,mat["RegionSizeTrain"]))


    # Explanation of what's stored in the "stats"
    # stats1[:,regi] = [Ssize;timer;Sstats[4];Rcond;pr;re;f1]

    # Save the Precision, Recall, and F1 scores
    numregions = size(stats1,2)
    prs = zeros(4,numregions)
    res = zeros(4,numregions)
    f1s = zeros(4,numregions)

    # Precision for 4 different versions of PageRank
    prs[1,:] = stats1[5,:]
    prs[2,:] = stats2[5,:]
    prs[3,:] = stats3[5,:]
    prs[4,:] = stats4[5,:]

    res[1,:] = stats1[6,:]
    res[2,:] = stats2[6,:]
    res[3,:] = stats3[6,:]
    res[4,:] = stats4[6,:]

    f1s[1,:] = stats1[7,:]
    f1s[2,:] = stats2[7,:]
    f1s[3,:] = stats3[7,:]
    f1s[4,:] = stats4[7,:]

    return prs, res, f1s, Sizes, seednodes

end

function LoadPPRTrain(alpha,seedtype)

    if seedtype == 1
        seedtag = "1"
    elseif seedtype == 2
        seedtag = "2"
    elseif seedtype == 3
        seedtag = "2per"
    elseif seedtype == 4
        seedtag = "3per"
    else
        println("Give a valid seedtype, 1 through 4")
        return
    end

    ppr_string = "../PPR_Train_Output/Seed_"*seedtag*"_alpha_"*string(round(Int64,alpha*100))*"sets_stats.mat"
    mat = matread(ppr_string)

    stats7 = mat["stats_tol_7"]
    sets7 = mat["sets_tol_7"]
    stats8 = mat["stats_tol_8"]
    sets8 = mat["sets_tol_8"]
    stats9 = mat["stats_tol_9"]
    sets9 = mat["sets_tol_9"]
    stats10 = mat["stats_tol_10"]
    sets10 = mat["sets_tol_10"]

    # Load the Region sizes
    mat = matread("../Brain_Regions/BrainRegions_Train23.mat")

    Sizes = vec(round.(Int64,mat["RegionSizeTrain"]))


    ## Explanation of what's stored in the "stats"
    # % [tol;numel(PR);tPR;cond;pPR;rPR;fPR];
    # % Row 1: epsilon tolerance (ppr locality parameter)
    # % Row 2: size of output
    # % Row 3: runtime
    # % Row 4: conductance of output
    # % Row 5: precision
    # % Row 6: recall
    # % Row 7: f1 score

    ## Save the Precision, Recall, and F1 scores
    numregions = size(stats7,2)
    prs = zeros(4,numregions)
    res = zeros(4,numregions)
    f1s = zeros(4,numregions)

    # Precision for 4 different versions of PageRank
    prs[1,:] = stats7[5,:]
    prs[2,:] = stats8[5,:]
    prs[3,:] = stats9[5,:]
    prs[4,:] = stats10[5,:]

    res[1,:] = stats7[6,:]
    res[2,:] = stats8[6,:]
    res[3,:] = stats9[6,:]
    res[4,:] = stats10[6,:]

    f1s[1,:] = stats7[7,:]
    f1s[2,:] = stats8[7,:]
    f1s[3,:] = stats9[7,:]
    f1s[4,:] = stats10[7,:]

    return prs, res, f1s, Sizes

end
