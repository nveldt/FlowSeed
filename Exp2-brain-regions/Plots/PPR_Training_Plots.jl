# Visualize output from running multiple experiments on regions of the brain
include("Plot_Helper.jl")
using Plots
using LaTeXStrings
# Specify which experiment you want
epsi = .1
for seedtype = 2:4
alpha = .6

prPPR, rePPR, f1PPR, Sizes = LoadPPRTrain(alpha,seedtype)

ppr_best,~ = findmax(f1PPR,1)

x_label = "Number of Nodes in Region"
y_label1 = "Precision"
y_label2 = "Recall"
y_label3 = "F1 Score"
numregions = size(f1s,2)
l_place = :false

# F1 Plots
lstyle = :solid
lw = 2
s1 = 400
s2 = 300

plot(f1PPR[1,:], labels = "\\varepsilon_{pr} = 10^{ -7}",grid = false,size = (s1,s2),
    xlabel = x_label,ylabel = y_label3,legend = l_place,linewidth = lw,
    xticks = (1:numregions,Sizes), xrotation=50,markershape = :circle,linestyle = lstyle,
)

plot!(f1PPR[2,:], labels = "\\varepsilon_{pr} = 10^{ -8}",
    xlabel = x_label,ylabel = y_label3,legend = l_place,linewidth = lw,size = (s1,s2),
    xticks = (1:numregions,Sizes), xrotation=50,markershape = :circle,linestyle = lstyle,
)

plot!(f1PPR[3,:], labels = "\\varepsilon_{pr} = 10^{ -9}",
    xlabel = x_label,ylabel = y_label3,legend = l_place,linewidth = lw,size = (s1,s2),
    xticks = (1:numregions,Sizes), xrotation=50,markershape = :circle,linestyle = lstyle,
)

plot!(f1PPR[4,:], labels = "\\varepsilon_{pr} = 10^{ -10}",
    xlabel = x_label,ylabel = y_label3,legend = l_place,linewidth = lw,size = (s1,s2),
    xticks = (1:numregions,Sizes), xrotation=50,markershape = :circle,linestyle = lstyle,
)

savefig("JustPush_SeedType_"*string(seedtype)*"alpha60.pdf")
end
