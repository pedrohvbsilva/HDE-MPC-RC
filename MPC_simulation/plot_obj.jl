using Plots

function plot_obj(a,b,c)
    p = plot()
    plot!(a, linewidth = 2, label = "Centralized")
    plot!(b, linewidth = 2, ls = :dash, label = "Bilevel")
    plot!(c, linewidth = 2, ls = :dash, label = "Benders")    
    ylabel!("Valor Objetivo")
    xlabel!("Iterações MPC")

    fn = plot(p)
    savefig(fn, "objs.pdf")
end
#  plot_obj(Objetivo_cent,Objetivo_bilevel,Objetivo_benders)