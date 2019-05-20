
using MAT

include("../algorithms/FlowSeed-0.6.jl")

# Set resolution parameters
epsi = 0.1

# Choose 2,3, or 5 percent of the target community to use.
# Relative results of algorithms are comparable for all choices of seed set size.
# Choosing 5% leads to larger seed sets and better results for all methods.
percentofset = 5

datasets = ["DBLP","Amazon","LiveJournal","Orkut"]

for graph = datasets

     # "nosoft" refers to the fact that we include no soft penalties on
     # excluding seed nodes from the output set
     outputstring = "Output"*string(percentofset)*"/"*graph*"slfs_"*string(epsi)*"epsnosoft_truecond.txt"
     outputmatrix = "Output"*string(percentofset)*"/"*graph*"slfs_"*string(epsi)*"epsnosotf_truecond.mat"

    open(outputstring,"w") do f
         write(f,"Flowseed. Strict penalties on excluding starter nodes
         (which are known to be in the target cluster) and soft penalties of 1 on excluding immediate neighbors.
         We assume we know $percentofset percent of the target cluster,
         and we set epsilon to $epsi.\n")
    end

    mat = matread("datasets/"*graph*"-top10.mat")
    A = mat["A"]

    matS = matread("datasets/"*graph*"-seed-starter.mat")
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

    fs_sets = spzeros(n,numcom)
    fs_stats = zeros(6,numcom)

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
        fR = volR/(volA - volR)
        epsilon = max(epsi,2*fR)
        epsilon = epsi
        STATS = set_stats(A,R,volA)
        condR = round(STATS[4],4)
        pr,re,f1 = round.(PRF(Target,R),4)
        println("Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1")

        open(outputstring,"a") do f
             write(f,"\nCommunity $commID has $Tsize nodes and a conductance of $Tcond \n")
             write(f, "Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1 \n")
        end

        ## Run some community detection experiments
        RinS = zeros(numR,1)
        pR = zeros(numR,1)
        # For standard flow seed procedure, include
        # Hard constraint on the starter seed nodes.
        for r = 1:numR
           if in(R[r],Rstart)
                RinS[r] = 1
           end
        end

        tic()
        S1, cond1 = FlowSeed(A,R,Rn,Rc,epsilon,pR,RinS,d,volA,volR,true)
        t1 = toc()
        c1 = round(cond1,4)
        size1 = length(S1)
        prf = PRF(Target,S1)
        pr1 = round(prf[1],3)
        re1 = round(prf[2],3)
        f11 = round(prf[3],3)
        Stats = set_stats(A,S1,volA)
        truecond = round(Stats[4],4)
        println("FlowSeed: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t Size = $size1 \t Cond = $c1 \t TrueCond = $truecond\t Time = $t1")

        fs_sets[S1,commID] = 1
        fs_stats[:,commID] = [t1; length(S1); c1; pr1; re1; f11]


        open(outputstring,"a") do f

          write(f, "Flow-Seed: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t Size = $size1 \t Cond = $c1 \t TrueCond = $truecond\t Time = $t1 \n")
        end

    end

    matwrite(outputmatrix, Dict( "fs_stats" => fs_stats,"fs_sets" => fs_sets))

end
