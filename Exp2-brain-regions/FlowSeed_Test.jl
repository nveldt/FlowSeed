include("algoriths/FlowSeed-0.6.jl")

using MAT

mat = matread("data/BrainGraph")
A = mat["B"]

n = size(A,1)
d = sum(A,1)'
volA = sum(nonzeros(A))

## Load Regions, these are the target clusters we want to recover.
# There are 17 in the training set and 78 in the testing set

mat = matread("Brain_Regions/BrainRegions_Test23.mat")

Regions = mat["RegionsTest"]
RegionSize = round.(Int64,mat["RegionSizeTest"])
RegionCond = mat["RegionCondTest"]
RegionLabels = round.(Int64,mat["RegionLabelTest"])
Seed2per = mat["Seed2perTest"]
Seed3per = mat["Seed3perTest"]
Seed1 = mat["Seed1Test"]
Seed2 = mat["Seed2Test"]

## Load the best results

mat = matread("Brain_Regions/Training17_Results.mat")

Epsilons = mat["epsis"]
TrainingSizes = mat["Sizes"]

# This is a mapping from cluster size, to best parameter to use for that cluster size
BestEps = mat["BestParamsFlowSeed"]


## Load a region, and its designated seed sets
# Seed sets were chosen and saved in advance for experimental reproduceability

for seedtype = 1:4

    # Load the appropriate seed set
    if seedtype == 1
        Seeds = Seed1
    elseif seedtype == 2
        Seeds = Seed2
    elseif seedtype == 3
        Seeds = Seed2per
        # epsilon is fixed, seed type is fixed (i.e. an integer from 1 to 4)
    elseif seedtype == 4
        Seeds = Seed3per
    else
        println("For seedtype: Input 1,2,3, or 4 (100 nodes, 1%, 2%, or 3% of target set)")
        return
    end

    outputstring = "FlowSeed_Test_Output/Seed_"*string(seedtype)*".txt"
    outputmatrix = "FlowSeed_Test_Output/Seed_"*string(seedtype)*".mat"

    open(outputstring,"w") do f
         write(f,"Brain Dataset testing region experiments for FlowSeed. Epsilon is chosen based on training experiments\n")

         if seedtype == 1
             write(f,"Seed sets are made up of a set of 100 initial starter nodes inside the target region, which we grow by their neighborhood.\n")
         elseif seedtype == 2
             write(f,"Seed sets are made up of a set of initial random starter nodes that make up 1% of the target region, which we grow by their neighborhood.\n")
         elseif seedtype == 3
             write(f,"Seed sets are made up of a set of initial random starter nodes that make up 2% of the target region, which we grow by their neighborhood.\n")
         else
             write(f,"Seed sets are made up of a set of initial random starter nodes that make up 3% of the target region, which we grow by their neighborhood.\n")
         end
         write(f,"We include strict constraints on excluding initial starter nodes from the seed set, and a soft penalty of p_r = 1 for excluding any other seed nodes.\n")
    end


        # Store output sets and stats for each set of experiments
        # There are 4 different types of penalties we will test
        sets = spzeros(n,78)
        stats = zeros(7,78)


        for regi = 1:size(Regions,2)

            # Load target set and output results
            Target = find(Regions[:,regi])
            Label = RegionLabels[regi]
            Size = length(Target)
            Cond = round(RegionCond[regi],4)
            println("\nRegion $Label has conductance $Cond and size $Size")
            open(outputstring,"a") do f
                 write(f,"\nRegion $Label has conductance $Cond and size $Size. \n")
            end
            T = Target

            # Determine how to set epsilon by comparing against training regions.
            # Find the training region (of 17 choices) with size closest to
            # current testing regions size
            m,closest_region = findmin(abs.(TrainingSizes-Size))

            tsize = TrainingSizes[closest_region]
            println("\n Test region is closest in size to a training region of size $tsize")
            open(outputstring,"a") do f
                 write(f,"It is closest in size to a training region with $tsize nodes.")
            end
            best_eps = BestEps[seedtype,closest_region]
            # Use whichever epsilon worked best for that TrainingSize
            epsilon = Epsilons[best_eps]
            println("Using epsilon = $epsilon")
            open(outputstring,"a") do f
                 write(f,"Using epsilon = $epsilon\n")
            end

            # Load the seed set, grow it by its neighborhood
            Rstart = find(Seeds[:,regi])
            R = neighborhood(A,Rstart,1)    # Grow "starter" seed by their neighborhood
            Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
            Rn = setdiff(Rn,R)          # ...but we exclude R itself
            inRc = ones(n)
            inRc[R] = 0
            Rc = find(inRc)             # complement of R

            cut,volR,edges,condR = set_stats(A,R,volA)
            volR = sum(d[R])
            numR = length(R)
            cR = round(condR,4)
            fR = volR/(volA - volR)
            println("\tSeed set: conductance = $cR, size = $numR")
            open(outputstring,"a") do f
                 write(f,"\tSeed set: conductance = $cR, size = $numR \n")
            end

            println("Penalty \t _Size_ \t _Time_ \t _Cond_ \t R-cond \t precision \t recall \t F1-score")
            open(outputstring,"a") do f
                 write(f,"Penalty \t _Size_ \t _Time_ \t _Cond_ \t R-cond \t precision \t recall \t F1-score  \n")
            end

            # Set up the soft and strict penalties
            RinS = zeros(numR,1)
            for i = 1:numR
                if in(R[i],Rstart)
                    RinS[i] = 1
                end
            end
            pR = ones(numR,1)

            # Run the algorithm
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
            println("FlowSeed \t\t $Ssize\t\t$timer\t\t $condS \t $RcondS \t $pr \t $re \t $f1 ")
            open(outputstring,"a") do f
                 write(f,"FlowSeed \t\t $Ssize\t\t$timer\t\t $condS \t $RcondS \t $pr \t $re \t $f1  \n")
            end
                # Save the output, both the set and the main output stats

            sets[S,regi] = 1
            stats[:,regi] = [Ssize;timer;Sstats[4];Rcond;pr;re;f1]

        matwrite(outputmatrix, Dict( "sets" => sets, "stats" => stats))

    end
end
