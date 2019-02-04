
using MAT

include("../algorithms/FlowSeed-0.6.jl")
using PyCall
@pyimport localgraphclustering as lgc

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
     outputstring = "Output"*string(percentofset)*"/"*graph*"slfs_"*string(epsi)*"epsnosoft.txt"
     outputmatrix = "Output"*string(percentofset)*"/"*graph*"slfs_"*string(epsi)*"epsnosotf.mat"

    open(outputstring,"w") do f
         write(f,"Flowseed and SimpleLocal. Strict penalties on excluding starter nodes
         (which are known to be in the target cluster) and soft penalties of 1 on excluding immediate neighbors.
         We assume we know $percentofset percent of the target cluster,
         and we set epsilon to $epsi.\n")
    end

    println(graph)
    tic()
    g = lgc.GraphLocal(homedir()*"/data/snap-top10/"*graph*".edgelist","edgelist"," ")
    toc()

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

    # fs_sl = FlowSeed-SimpleLocal, i.e. does the same thing as simplelocal,
    # and we just want to run it to check runtime.
    fs_sets = spzeros(n,numcom)
    fs_sl_sets = spzeros(n,numcom)
    sl_sets = spzeros(n,numcom)
    fs_stats = zeros(6,numcom)
    sl_stats = zeros(6,numcom)
    fs_sl_stats = zeros(6,numcom)

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
        println("FlowSeed: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t Size = $size1 \t Cond = $c1 \t Time = $t1")

        fs_sets[S1,commID] = 1
        fs_stats[:,commID] = [t1; length(S1); c1; pr1; re1; f11]


        # Now for runtime comparison against SimpleLocal, use no penalties
        RinS = zeros(numR,1)
        pR = zeros(numR,1)
        tic()
        S2, cond2 = FlowSeed(A,R,Rn,Rc,epsilon,pR,RinS,d,volA,volR,false)
        t2 = toc()
        c2 = round(cond2,4)
        size2 = length(S2)
        prf = PRF(Target,S2)
        pr2 = round(prf[1],3)
        re2 = round(prf[2],3)
        f12 = round(prf[3],3)
        println("FS SimLoc: \tPR = $pr2 \t RE = $re2 \t F1 = $f12 \t Size = $size2 \t Cond = $c2 \t Time = $t2")

        fs_sl_sets[S2,commID] = 1
        fs_sl_stats[:,commID] = [t1; length(S2); c2; pr2; re2; f12]

        ## Now run the LocalGraphClustering implementation of SimpleLocal (SL)
        del = epsilon-fR
        @show maximum(R)
        tic()
        simplelocal = lgc.flow_clustering(g,R-1,method="sl",delta=del)
        tSL = toc()
        SL = convert(Array{Int64},simplelocal[1])+1
        STATS = set_stats(A,SL,volA)
        condSL = round(STATS[4],4)
        si = length(SL)
        prf = PRF(Target,SL)
        pr = round(prf[1],3)
        re = round(prf[2],3)
        f1 = round(prf[3],3)
        sl_sets[SL,commID] = 1
        sl_stats[:,commID] = [tSL; si; condSL; pr; re; f1]
        println("SimpleLocal \tPR = $pr \t RE = $re \t F1 = $f1, Size = $si, Cond = $condSL, Time = $tSL")

        println("\nCommunity $commID has $Tsize nodes and a conductance of $Tcond ")
        # println("Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1")
        println("SimpleLocal \tPR = $pr \t RE = $re \t F1 = $f1, Size = $si, Cond = $condSL, Time = $tSL")
        println("FS SimLoc: \tPR = $pr2 \t RE = $re2 \t F1 = $f12 \t Size = $size2 \t Cond = $c2 \t Time = $t2")
        println("Flow-Seed: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t Size = $size1 \t Cond = $c1 \t Time = $t1")

        open(outputstring,"a") do f
          write(f, "SimpleLoc: \tPR = $pr  \t RE = $re  \t F1 = $f1  \t Size = $si    \t Cond = $condSL \t Time = $tSL \n")
          write(f, "FS SimLoc: \tPR = $pr2 \t RE = $re2 \t F1 = $f12 \t Size = $size2 \t Cond = $c2 \t Time = $t2 \n")
          write(f, "Flow-Seed: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t Size = $size1 \t Cond = $c1 \t Time = $t1 \n")
        end

    end

    matwrite(outputmatrix, Dict( "fs_stats" => fs_stats, "fs_sl_stats" => fs_sl_stats, "sl_stats" => sl_stats,
        "fs_sets" => fs_sets, "fs_sl_sets" => fs_sl_sets, "sl_sets" => sl_sets))

end
