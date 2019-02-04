# Run CRD using default parameters on the snap datasets

using MAT

include("FlowSeed.jl")
using PyCall
@pyimport localgraphclustering as lgc

function snapexp(percentofset::Int64,IndividualSeeds::Bool = false)

datasets = ["Amazon","DBLP","LiveJournal","Orkut","Youtube"]

for graph = datasets

    if IndividualSeeds
        outputstring = "Output"*string(percentofset)*"/"*graph*"_crd_individual.txt"
        outputmatrix = "Output"*string(percentofset)*"/"*graph*"_crd_individual.mat"
    else
        outputstring = "Output"*string(percentofset)*"/"*graph*"_crd.txt"
        outputmatrix = "Output"*string(percentofset)*"/"*graph*"_crd.mat"
    end

    open(outputstring,"w") do f
         write(f,"CRD. We assume we know $percentofset percent of the target cluster.\n")
    end

    # Load the edgelist of the graph
    tic()
    g = lgc.GraphLocal(homedir()*"/data/snap-top10/"*graph*".edgelist","edgelist"," ")
    toc()

    # Load the communities
    mat = matread(homedir()*"/data/snap-top10/"*graph*"-top10.mat")
    A = mat["A"]

    # Load the seed sets
    matS = matread(homedir()*"/data/snap-top10/"*graph*"-seed-starter.mat")
    if percentofset == 2
         S = matS["S2"]
    elseif percentofset == 3
         S = matS["S3"]
    elseif percentofset == 5
         S = matS["S5"]
    else
         println("Choose 2,3, or 5 percent of the target set. Defaulting to 5.")
         S = mat["S5"]
    end

    comm = mat["C"]
    Conds = mat["Conds"]
    d = sum(A,1)'
    volA = sum(nonzeros(A))

    n = size(A,1)
    d = sum(A,1)'
    println("nodes = $n")
    numcom = size(comm,2)

    crd_sets = spzeros(n,numcom)
    crd_b_sets = spzeros(n,numcom)
    crd_stats = zeros(6,numcom)
    crd_b_stats = zeros(6,numcom)

    for commID = collect(1:numcom)

        ## Choose one community to be the target cluster
        Target = find(comm[:,commID])
        TarSize = length(Target)
        TargetStats = set_stats(A,vec(Target),volA)
        Tcond = round(TargetStats[4],3)
        Tsize = length(Target)
        println("Community $commID has $Tsize nodes and a conductance of $Tcond")

        ## Get the starter for the seed set for this community
        Rstart = find(S[:,commID])
        R = neighborhood(A,Rstart,1)      # Grow by a one-hope neighborhood
        Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
        Rn = setdiff(Rn,R)          # ...but we exclude R itself
        inRc = ones(n)
        inRc[R] = 0
        Rc = find(inRc)             # complement of R
        numR = length(R)
        volR = sum(d[R])
        STATS = set_stats(A,R,volA)
        condR = round(STATS[4],4)
        pr,re,f1 = round.(PRF(Target,R),4)
        println("Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1")

        open(outputstring,"a") do f
             write(f,"\nCommunity $commID has $Tsize nodes and a conductance of $Tcond \n")
             write(f, "Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1 \n")
        end

        if IndividualSeeds
            # Run CRD trying each known target node as a seed
            # (Takes a long time and returns poor results)
            ## Next do CRD with default parameters and the full seed set
            condCRD = 1
            crd_S = []
            tic()
            for rr = Rstart
                 crd = lgc.flow_clustering(g,[rr-1],method="crd")
                 crd_temp = convert(Array{Int64},crd[1])+1
                 STATS = set_stats(A,crd_temp,volA)
                 cond_temp = round(STATS[4],3)

                 if cond_temp < condCRD
                      condCRD = cond_temp
                      crd_S = crd_temp
                 end
            end
            tCRD = toc()

            # Output some stats
            sizeCRD = length(crd_S)
            prf = PRF(Target,crd_S)
            prCRD = round(prf[1],3)
            reCRD = round(prf[2],3)
            f1CRD = round(prf[3],3)
            crd_sets[crd_S,commID] = 1
            if isnan(f1CRD)
                 f1CRD = 0
            end
            crd_stats[:,commID] = [tCRD; sizeCRD; condCRD; prCRD; reCRD; f1CRD]
            println("CRD Indiv: \tPR = $prCRD \t RE = $reCRD \t F1 = $f1CRD, Size = $sizeCRD, Cond = $condCRD \t Time = $tCRD")

            open(outputstring,"a") do f
                 write(f, "CRD Indiv: \tPR = $prCRD \t RE = $reCRD \t F1 = $f1CRD \t Size = $sizeCRD \t Cond = $condCRD \t Time = $tCRD\n")
            end

        else

            ## Run CRD with default parameters and the full seed set (best results)
            tic()
            crd = lgc.flow_clustering(g,R-1,method="crd",U=3,h=10,w=2,iterations=20)
            tCRD = toc()
            crd_S = convert(Array{Int64},crd[1])+1
            STATS = set_stats(A,crd_S,volA)
            condCRD = round(STATS[4],4)
            sizeCRD = length(crd_S)

            # Output some stats
            prf = PRF(Target,crd_S)
            prCRD = round(prf[1],3)
            reCRD = round(prf[2],3)
            f1CRD = round(prf[3],3)
            crd_sets[crd_S,commID] = 1
            if isnan(f1CRD)
                 f1CRD = 0
            end
            crd_stats[:,commID] = [tCRD; sizeCRD; condCRD; prCRD; reCRD; f1CRD]
            println("CRD Stand: \tPR = $prCRD \t RE = $reCRD \t F1 = $f1CRD, Size = $sizeCRD, Cond = $condCRD \t Time = $tCRD")

            open(outputstring,"a") do f
                 write(f, "CRD Stand: \tPR = $prCRD \t RE = $reCRD \t F1 = $f1CRD \t Size = $sizeCRD \t Cond = $condCRD \t Time = $tCRD\n")
            end
        end

    end
    matwrite(outputmatrix, Dict("crd_sets" => crd_sets,"crd_stats" => crd_stats))

end

end
