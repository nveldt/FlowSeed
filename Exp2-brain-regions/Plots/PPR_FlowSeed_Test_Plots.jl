using Plots

# Specify which experiment you want
alpha = .6
seedtype = 1
# Load the data for that experiment
using MAT
ppr_string = "../PPR_Test_Output/Seed_"*string(seedtype)*"_alpha_"*string(round(Int64,alpha*100))*"sets_stats.mat"
matPPR = matread(ppr_string)

# stats = [tol;numel(PR);tPR;cond;pPR;rPR;fPR];
ppr_stats = matPPR["stats"]

fs_string = "../FlowSeed_Test_Output/Seed_"*string(seedtype)*".mat"
matFS = matread(fs_string)

# [Ssize;timer;Sstats[4];Rcond;pr;re;f1]
fs_stats = matFS["stats"]

flowseed_f1 = fs_stats[7,1:2:end]
ppr_f1 = ppr_stats[7,1:2:end]

mat = matread("../Brain_Regions/BrainRegions_Test.mat")
Sizes = vec(round.(Int64,mat["RegionSizeTest"]))
Sizes = Sizes[1:2:end]

x_label = "Number of Nodes in Region"
y_label = "F1 Score"
numregions = length(ppr_f1)
l_place = :bottomleft
lstyle = :solid
lw = 1.5
s1 = 400
s2 = 300
ms = 3
xt = 1:3:numregions
sz = Sizes[xt]

plot(flowseed_f1, labels = "FlowSeed",grid = false,size = (s1,s2),
    xlabel = x_label,ylabel = y_label,legend = l_place,linewidth = lw,
    xticks = (xt,sz), xrotation=50,markershape = :circle,markersize = ms,
    linestyle = lstyle,color = :blue,ylim = [0,1]
)

plot!(ppr_f1, labels = "Push",grid = false,size = (s1,s2),
    xlabel = x_label,ylabel = y_label,legend = l_place,linewidth = lw,
    xticks = (xt,sz), xrotation=50,markershape = :cross,markersize = ms,
    linestyle = lstyle, color = :green,ylim = [0,1]
)



savefig("TestPlots_Seed_"*string(seedtype)*".pdf")
