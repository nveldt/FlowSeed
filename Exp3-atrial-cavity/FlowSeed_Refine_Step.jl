# Run FlowSeed on the full body MRI graph, refining PPR outputs

include("../algorithms/FlowSeed-0.6.jl")
using MAT

mat = matread("data/MRI_Graph.mat")
A = mat["mri_A"]
Target = mat["cavity"]
Target = round.(Int64,Target)
volA = mat["mri_vol"]
d = sum(A,1)'
n = size(A,1)

# Output stats about the target atrial cavity
TarSize = length(Target)
TargetStats = set_stats(A,vec(Target),volA)
Tcond = round(TargetStats[4],3)
Tsize = length(Target)
println("Target Cavity has $Tsize nodes and a conductance of $Tcond ")

mat = matread("data/cavity_100_nodes.mat")
seed_start = mat["seed_start"] # Seed100

alpha = .6
epsilon = .1

# These 100 nodes we know to be in the atrial cavity
R0 = vec(round.(Int64,seed_start))

# One-hop neighborhood
R1 = neighborhood(A,R0,1)

# Load the output from PPR
mat = matread("PPR_Output/alpha_"*string(round(Int64,100*alpha))*"sets_stats.mat")

sets = mat["sets"]
stats = mat["stats"]

outputstring = "FlowSeed_refine_PPR_alpha_"*string(alpha)*"_epsilon_"*string(epsilon)*".txt"
outputmat = "FlowSeed_refine_PPR_alpha_"*string(alpha)*"_epsilon_"*string(epsilon)*".mat"

open(outputstring,"w") do f
    write(f,"Running FlowSeed on the PPR push algorithm output.\n")
end

sets1 = spzeros(n,4)
stats1 = zeros(7,4)
sets2 = spzeros(n,4)
stats2 = zeros(7,4)

# Refine the output of PPR for a subset of the PPR tolerances used
for t = 2:size(sets,2)-3

    tol = t-1

    # Output from PPR
    Rppr = find(sets[:,t])
    prf = PRF(Target,Rppr)
    pr = round(prf[1],3)
    re = round(prf[2],3)
    f1 = round(prf[3],3)

    # The seed set for FlowSeed will be the output from PPR, plus
    # the original 100 nodes from the atrial cavity, and their neighbors
    R = unique([R0;R1;Rppr])
    Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
    Rn = setdiff(Rn,R)          # ...but we exclude R itself
    inRc = ones(n)
    inRc[R] = 0
    Rc = find(inRc)             # complement of R

    pprSize = length(Rppr)
    open(outputstring,"a") do f
        write(f,"\nPPR tolerance Number $tol in [1e-9, 1e-10, 1e-11]\n")
        write(f,"PPR output: $pprSize nodes, PR = $pr, RE = $re, F1 = $f1 \n")
    end

    numR = length(R)
    prf = PRF(Target,R)
    pr = round(prf[1],3)
    re = round(prf[2],3)
    f1 = round(prf[3],3)
    println("Seed set: \tPR = $pr \t RE = $re \t F1 = $f1, Size = $numR")
    open(outputstring,"a") do f
        write(f,"Seed set: $numR nodes, PR = $pr, RE = $re, F1 = $f1 \n")
    end

    println("Penalty \t _Size_ \t _Time_ \t _Cond_ \t R-cond \t precision \t recall \t F1-score")
    open(outputstring,"a") do f
         write(f,"Penalty \t _Size_ \t _Time_ \t _Cond_ \t R-cond \t precision \t recall \t F1-score  \n")
    end

    # Type 1 seed exclusion penalties
    RinS = zeros(numR,1)
    pR = zeros(numR,1)
    # Hard constraint on the starter seed nodes.
    # Soft penalty on the 1st neighborhood
    for r = 1:numR
        if in(R[r],R0)
            RinS[r] = 1
        end
        if in(R[r],R1)
            pR[r] = 1
        end
    end

    # FlowSeed
    volR = sum(d[R])

    tic()
    S, Rcond = FlowSeed(A,R,Rn,Rc,epsilon,pR,RinS,d,volA,volR,true)
    timerfull = toq()
    timer = round(timerfull,3)
    Sstats = set_stats(A,S,volA)
    condS = round(Sstats[4],4)
    RcondS = round(Rcond,4)
    Ssize = length(S)
    prf = PRF(Target,S)
    pr = round(prf[1],4)
    re = round(prf[2],4)
    f1 = round(prf[3],4)
    println("Type 1 \t $Ssize\t$timer\t $condS \t $RcondS \t $pr \t $re \t $f1 ")
    open(outputstring,"a") do f
         write(f,"Type 1 \t $Ssize\t$timer\t $condS \t $RcondS \t $pr \t $re \t $f1 \n")
    end

    sets1[S,tol] = 1
    stats1[:,tol] = [Ssize;timerfull;Sstats[4];Rcond;pr;re;f1]

    # Type 2 seed exclusion penalties, more strict
    RinS = zeros(numR,1)
    pR = zeros(numR,1)
    # Hard constraint on the starter seed nodes.
    # Soft penalty on the 1st neighborhood
    # Softer penalty on ppr output
    for r = 1:numR
        if in(R[r],R0)
            RinS[r] = 1
        end
        if in(R[r],R1)
            pR[r] = 1
        end
        if in(R[r],Rppr)
            pR[r] = .5
        end
    end

    tic()
    S, Rcond = FlowSeed(A,R,Rn,Rc,epsilon,pR,RinS,d,volA,volR,true)
    timerfull = toq()
    timer = round(timerfull,3)
    Sstats = set_stats(A,S,volA)
    condS = round(Sstats[4],4)
    RcondS = round(Rcond,4)
    Ssize = length(S)
    prf = PRF(Target,S)
    pr = round(prf[1],4)
    re = round(prf[2],4)
    f1 = round(prf[3],4)
    println("Type 2 \t $Ssize\t$timer\t $condS \t $RcondS \t $pr \t $re \t $f1 ")
    open(outputstring,"a") do f
         write(f,"Type 2 \t $Ssize\t$timer\t $condS \t $RcondS \t $pr \t $re \t $f1 \n")
    end
    sets2[S,tol] = 1
    stats2[:,tol] = [Ssize;timerfull;Sstats[4];Rcond;pr;re;f1]

end

matwrite(outputmat, Dict( "sets1" => sets1, "stats1" => stats1,
"sets2" => sets2, "stats2" => stats2))
