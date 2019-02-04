# Test the FlowSeed local graph clustering implementation

using MAT

include("algorithms/FlowSeed-0.6.jl")

mat = matread("Exp1-snap-communities/datasets/DBLP-top10.mat")
A = mat["A"]

matS = matread("Exp1-snap-communities/datasets/DBLP-seed-starter.mat")
S = matS["S5"]

comm = mat["C"]
Conds = mat["Conds"]
d = sum(A,1)'
volA = sum(nonzeros(A))

n = size(A,1)
d = sum(A,1)'
println("nodes = $n")
numcom = size(comm,2)

for commID = 1:10

# Extract one of the target communities
Target = find(comm[:,commID])
TarSize = length(Target)
TargetStats = set_stats(A,vec(Target),volA)
Tcond = round(TargetStats[4],3)
Tsize = length(Target)
println("Community $commID has $Tsize nodes and a conductance of $Tcond")

## Get a seed set for the community and report some stats about it
Rstart = find(S[:,commID])
R = neighborhood(A,Rstart,1)      # Grow by a one-hope neighborhood
epsilon = 0.1
STATS = set_stats(A,R,volA)
condR = round(STATS[4],4)
pr,re,f1 = round.(PRF(Target,R),4)
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
tic()
S1, cond1 = FlowSeed(A,R,epsilon,pR,RinS)
t1 = toc()
c1 = round(cond1,4)
size1 = length(S1)
prf = PRF(Target,S1)
pr1 = round(prf[1],3)
re1 = round(prf[2],3)
f11 = round(prf[3],3)
print("FlowSeed: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t Size = $size1 \t ")
println("Cond = $c1 \t Time = $t1")
end
