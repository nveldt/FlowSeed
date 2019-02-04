# Test the FlowSeed local graph clustering implementation

using MAT

include("algorithms/FlowSeed-1.0.jl")

mat = matread("Exp1-snap-communities/datasets/DBLP-top10.mat")
A = mat["A"]

matS = matread("Exp1-snap-communities/datasets/DBLP-seed-starter.mat")
S = matS["S5"]

comm = mat["C"]
Conds = mat["Conds"]
d = sum(A,dims = 1)'
volA = sum(A.nzval)

n = size(A,1)
println("nodes = $n")
numcom = size(comm,2)

for commID = 1:10

# Extract one of the target communities
Target = findall(x->x!=0,comm[:,commID])
TarSize = length(Target)
TargetStats = set_stats(A,vec(Target),volA)
Tcond = round(TargetStats[4],digits = 3)
Tsize = length(Target)
println("Community $commID has $Tsize nodes and a conductance of $Tcond")

## Get a seed set for the community and report some stats about it
Rstart = findall(x->x!=0,S[:,commID])
R = neighborhood(A,Rstart,1)      # Grow by a one-hope neighborhood
epsilon = 0.1
STATS = set_stats(A,R,volA)
condR = round(STATS[4],digits = 4)
pr,re,f1 = round.(PRF(Target,R),digits = 4)
numR = length(R)
println("Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1")

# Set up the penalties on excluding seed nodes
RinS = zeros(numR,1)
pR = zeros(numR,1)
for r = 1:numR
  # Strictly penalize the exclusion of nodes that are known
  # to be in the target set
  if in(R[r],Rstart)
    RinS[r] = 1
  else
    # Place a soft penalty on excluding other seed nodes
    pR[r] = 1
  end
end

# Run FlowSeed and check how well it is able to detect the target community
start = time()
S1, cond1 = FlowSeed(A,R,epsilon,pR,RinS)
t1 = time() - start
c1 = round(cond1,digits = 4)
size1 = length(S1)
prf = PRF(Target,S1)
pr1 = round(prf[1],digits = 3)
re1 = round(prf[2],digits = 3)
f11 = round(prf[3],digits = 3)
print("FlowSeed: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t Size = $size1 \t ")
println("Cond = $c1 \t Time = $t1")

end
