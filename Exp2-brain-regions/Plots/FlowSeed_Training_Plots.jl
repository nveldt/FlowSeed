# Visualize output from running multiple experiments on regions of the brain
include("Plot_Helper.jl")
using Plots
using LaTeXStrings

epsi = .1
numregions = 17

for seedtype = 1:4

    prs1, res1, f1s1, Sizes, seednodes =  LoadFlowSeedTrain(epsi,seedtype)

    x_label = "Number of Nodes in Region"
    y_label1 = "Precision"
    y_label2 = "Recall"
    y_label3 = "F1 Score"

    l_place = :bottomright
    s1 = 400
    s2 = 300
    ms = 4
    lw = 2

    # F1 Plots
    plot(f1s1[1,:], labels = "No Penalties",ylim = [0,1],grid = false,size = (s1,s2),
        xlabel = x_label,ylabel = y_label3,legend = l_place,
        xticks = (1:numregions,Sizes), xrotation=50,linewidth = lw,
        color = :blue,markershape = :circle, markersize = ms)

    plot!(f1s1[2,:], labels = "Some Strict",size = (s1,s2),
        color = :red,markershape = :diamond,linewidth = lw,markersize = ms)

    plot!(f1s1[4,:], labels = "Strict+Soft",size = (s1,s2),
        color = :green,markershape = :square,linewidth = lw,markersize = ms)


    savefig("FlowSeed_Train_Seed_"*string(seedtype)*"_eps_"*string(epsi)*"_F1.pdf")

    # Recall Plots
    plot(res1[1,:], labels = "No Penalties",ylim = [0,1],grid = false,size = (s1,s2),
        xlabel = x_label,ylabel = y_label3,legend = l_place,
        xticks = (1:numregions,Sizes), xrotation=50,linewidth = lw,
        color = :blue,markershape = :circle, markersize = ms)

    plot!(res1[2,:], labels = "Some Strict",size = (s1,s2),
        color = :red,markershape = :diamond,linewidth = lw,markersize = ms)

    plot!(res1[4,:], labels = "Strict+Soft",size = (s1,s2),
        color = :green,markershape = :square,linewidth = lw,markersize = ms)


    savefig("FlowSeed_Train_Seed_"*string(seedtype)*"_eps_"*string(epsi)*"_Recall.pdf")

end
