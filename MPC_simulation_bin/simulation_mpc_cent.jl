using Plots, LinearAlgebra, Random
using JuMP, Gurobi
using DelimitedFiles
# using MAT

Random.seed!(5)

include("problem.jl")
include("create_m_models.jl")
include("sdp_qp_cent.jl")

global const GRB_ENV = Gurobi.Env()

T = 25;
G = Array{NamedTuple}(undef, nM);
SS = Array{NamedTuple}(undef, nM, T);
Objetivo_cent = Float64[];

u0 = zeros(nM, nU);
xm = zeros(nM, nX, T);
xm[:,:,1] = rand(nM, nX);
ym = zeros(nM, nY, T-1);

G = create_m_models(A, B, C, D, nM, nX, nY, nU, nR, nT, xm[:,:,1], u0);

for m in 1:nM
    SS[m, 1] = (du = 0., u = u0[m,1], y = 0., x = xm[m,:,1]);
end

tempo = @elapsed let SS = SS, G = G, ym = ym, xm = xm
    for i in 2:T
 
        obj, Sm = sdp_qp_cent(G);
         # global S;
        @show i
        push!(Objetivo_cent, obj)

        for m in 1:G.nM
            SS[m,i] = (du = Sm.Du[m,1,1], u = Sm.u[m,1,2], y = Sm.y[m,1,2], x = Sm.x[m,:,2]);

            xm[m, :, i] = A .* xm[m ,:, i-1] .+ B .* SS[m, i-1].u;
            ym[m, :, i-1] = C .* xm[m, :, i-1] .+ D .* SS[m, i-1].u;

            u0[m,1] = SS[m,i].u;
        end
        # if i == 10
        #     nM = 5
        # end
        G = create_m_models(A, B, C, D, nM, nX, nY, nU, nR, nT, xm[:, :, i], u0);
    end
end
;
println("Elapsed T: ", tempo, "(s)")

upl = zeros(nM, T);
ypl = zeros(nM, T);
dupl = zeros(nM, T);
wmpl = 1*ones(1,T);

for m in 1:G.nM
    for i in 1:T-1
        upl[m,i] = SS[m,i].u;
        ypl[m,i] = ym[m, 1, i];
        dupl[m,i] = SS[m,i].du;
    end
end

# matwrite("var_plot_cent.mat", Dict("u" => upl, "y"=>ypl, "obj_cent"=>Objetivo_cent))

# pyplot()
# PyPlot.pygui(true)
# plotlyjs()
gr()
p1 = plot();
for i in 1:nM
    plot!((1:T-1).*Tsh, upl[i,2:end], linewidth = 2);
end
plot!(xlabel="tempo (h)", ylabel="Sinal de Controle");

# p2 = plot();
# for i in 1:nM
#     plot!((1:T).*(Ts/60),dupl[i,:],label="du$i", linewidth = 2);
# end
# plot!(xlabel="tempo (min)",ylabel="Sinal da Variação Controle");

p3 = plot();
for i in 1:nM
    plot!((1:T-1).*Tsh, ypl[i,1:end-1], linewidth = 2);
end
plot!((1:T-1).*Tsh, wmpl[1:end-1]);
plot!(xlabel="tempo (h)", ylabel="Saida do Sistema");

l = @layout [a; b]
fn = plot(p1, p3, layout = l)
savefig(fn, "m4cent.pdf")

writedlm("centy.csv", ypl,',')
writedlm("centu.csv", upl,',')
writedlm("centobj.csv", Objetivo_cent,',')