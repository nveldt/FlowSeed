include("algoriths/FlowSeed-0.6.jl")

using MAT

mat = matread("data/BrainGraph")
A = mat["B"]

n = size(A,1)
d = sum(A,1)'
volA = sum(nonzeros(A))

mat = matread("Brain_Regions/BrainRegions_Train23.mat")

Regions = mat["RegionsTrain"]
RegionSize = round.(Int64,mat["RegionSizeTrain"])
RegionCond = mat["RegionCondTrain"]
RegionLabels = round.(Int64,mat["RegionLabelTrain"])
Seed2per = mat["Seed2perTrain"]
Seed3per = mat["Seed3perTrain"]

mat = matread("Brain_Regions/BrainRegions_Train.mat")
Seed1 = mat["Seed1Train"]
Seed2 = mat["Seed2Train"]

## Load a region, and its designated seed sets
# Seed sets were chosen and saved in advance for experimental reproduceability

# Set a single value for epsilon
epsilon = .1
tag = ""
for seedtype = 1:4

    TTA = "Train"
    if seedtype == 1
        Seeds = Seed1
        outputstring = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype)*"_epsilon_"*string(epsilon)*tag*".txt"
        outputmatrix = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype)*"_epsilon_"*string(epsilon)*tag*".mat"

    elseif seedtype == 2
        Seeds = Seed2
        outputstring = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype)*"_epsilon_"*string(epsilon)*tag*".txt"
        outputmatrix = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype)*"_epsilon_"*string(epsilon)*tag*".mat"

    elseif seedtype == 3
        Seeds = Seed2per
        # epsilon is fixed, seed type is fixed (i.e. an integer from 1 to 4)
        outputstring = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype-1)*"per_epsilon_"*string(epsilon)*tag*".txt"
        outputmatrix = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype-1)*"per_epsilon_"*string(epsilon)*tag*".mat"
    elseif seedtype == 4
        Seeds = Seed3per
        outputstring = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype-1)*"per_epsilon_"*string(epsilon)*tag*".txt"
        outputmatrix = "FlowSeed_"*TTA*"_Output/Seed_"*string(seedtype-1)*"per_epsilon_"*string(epsilon)*tag*".mat"
    else
        println("For seedtype: Input 1,2,3, or 4 (100 nodes, 1%, 2%, or 3% of target set)")
        return
    end

    open(outputstring,"w") do f
         write(f,"Brain Dataset training region experiments for FlowSeed when epsilon = $epsilon\n")

         if seedtype == 1
             write(f,"Seed sets are made up of a set of 100 initial starter nodes inside the target region, which we grow by their neighborhood.\n")
         elseif seedtype == 2
             write(f,"Seed sets are made up of a set of initial random starter nodes that make up 1% of the target region, which we grow by their neighborhood.\n")
         elseif seedtype == 3
             write(f,"Seed sets are made up of a set of initial random starter nodes that make up 2% of the target region, which we grow by their neighborhood.\n")
         else
             write(f,"Seed sets are made up of a set of initial random starter nodes that make up 3% of the target region, which we grow by their neighborhood.\n")
         end
         write(f,"We use four different types of penalties on excluding seed nodes:\n")
         write(f,"Type 1: there is no penalty on excluding seed nodes, but no additional soft penalties.\n")
         write(f,"Type 2: there is a strict penalty on excluding any nodes from the initial starter set.\n")
         write(f,"Type 3: strict penalty for excluding starter set. Soft penalties for other nodes, of weight neighbs/2 \nwhere neighbs[r] = number of nodes in the starter set that are adjacent to node r in the seed set.\n")
         write(f,"Type 4: strict penalty for excluding starter set. Soft penalty of 1 for every other node.\n")
    end

        # Store output sets and stats for each set of experiments
        # There are 4 different types of penalties we will test
        sets1 = spzeros(n,17)
        sets2 = spzeros(n,17)
        sets3 = spzeros(n,17)
        sets4 = spzeros(n,17)

        stats1 = zeros(7,17)
        stats2 = zeros(7,17)
        stats3 = zeros(7,17)
        stats4 = zeros(7,17)

        for regi = 1:size(Regions,2)

            # Load target set and output results
            Target = find(Regions[:,regi])
            Label = RegionLabels[regi]
            Size = length(Target)
            Cond = round(RegionCond[regi],4)
            println("\nRegion $Label has conductance $Cond and size $Size")
            open(outputstring,"a") do f
                 write(f,"\nRegion $Label has conductance $Cond and size $Size \n")
            end
            T = Target

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

            for Type = 1:4

                # Set up four different versions of strict and soft penalties
                if Type == 1
                    pR = zeros(numR,1)
                    RinS = zeros(numR,1)
                elseif Type == 2
                    pR = zeros(numR,1)
                    RinS = zeros(numR,1)
                    for i = 1:numR
                        if in(R[i],Rstart)
                            RinS[i] = 1
                        end
                    end
                elseif Type == 3
                    pR = zeros(numR,1)
                    NeighbTime = zeros(n,1)    # count the number of nodes in Rstart are adjacent to each node
                    for i = 1:length(Rstart)
                        N = neighborhood(A,[Rstart[i]],1)
                        NeighbTime[N] += 1
                    end
                    ~, ~, times = findnz(NeighbTime)
                    pR = times/2    # penalize proportional to the percentage of nodes in Rstart that neighbor it

                    RinS = zeros(numR,1)
                    for i = 1:numR
                        if in(R[i],Rstart)
                            RinS[i] = 1
                        end
                    end
                else
                    RinS = zeros(numR,1)
                    for i = 1:numR
                        if in(R[i],Rstart)
                            RinS[i] = 1
                        end
                    end
                    pR = ones(numR,1)
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
                println("Type $Type \t\t $Ssize\t\t$timer\t\t $condS \t $RcondS \t $pr \t $re \t $f1 ")
                open(outputstring,"a") do f
                     write(f,"Type $Type \t\t $Ssize\t\t$timer\t\t $condS \t $RcondS \t $pr \t $re \t $f1  \n")
                end
                # Save the output, both the set and the main output stats

                if Type == 1
                    sets1[S,regi] = 1
                    stats1[:,regi] = [Ssize;timer;Sstats[4];Rcond;pr;re;f1]
                elseif Type ==2
                    sets2[S,regi] = 1
                    stats2[:,regi] = [Ssize;timer;Sstats[4];Rcond;pr;re;f1]
                elseif Type ==3
                    sets3[S,regi] = 1
                    stats3[:,regi] = [Ssize;timer;Sstats[4];Rcond;pr;re;f1]
                else
                    sets4[S,regi] = 1
                    stats4[:,regi] = [Ssize;timer;Sstats[4];Rcond;pr;re;f1]
                end
            end
        end

        matwrite(outputmatrix, Dict( "sets1" => sets1, "stats1" => stats1,
            "sets2" => sets2, "stats2" => stats2,"sets3" => sets3,
            "stats3" => stats3,"sets4" => sets4, "stats4" => stats4))

end
