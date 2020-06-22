using Plots
using DelimitedFiles

function plot_obj(a,b,c)
    p = plot()
    plot!(a, linewidth = 2, label = "Centralized")
    plot!(b, linewidth = 2, ls = :dash, label = "OA")
    plot!(c, linewidth = 2, ls = :dash, label = "Benders")    
    ylabel!("Valor Objetivo")
    xlabel!("Iterações MPC")

    fn = plot(p)
    savefig(fn, "objs.pdf")
end
#  plot_obj(Objetivo_cent,Objetivo_OA,Objetivo_benders)
Obj_cent = readdlm("battery_charging\\teste pc casa T=5 seed=5\\centobj.csv")
Obj_bend = readdlm("battery_charging\\teste pc casa T=5 seed=5\\bendersobj.csv")
Obj_OA = readdlm("battery_charging\\teste pc casa T=5 seed=5\\OAobj.csv")

plot_obj(Obj_cent,Obj_OA,Obj_bend)