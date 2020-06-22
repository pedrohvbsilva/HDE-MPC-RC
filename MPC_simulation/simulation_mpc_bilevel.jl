using Plots, LinearAlgebra, Random
using JuMP, Ipopt, Gurobi
using DelimitedFiles
# using ProgressMeter
# using MAT

Random.seed!(2)

include("problem.jl")
include("create_m_models.jl")
include("sdp_bilevel_master.jl")
include("sdp_bilevel_sm.jl")
include("sdp_bilevel_s.jl")

global const GRB_ENV = Gurobi.Env()

T = 40;
G = Array{NamedTuple}(undef, nM);
SS = Array{NamedTuple}(undef, nM, T);
Objetivo_bilevel = Float64[];

u0 = zeros(nM, nU);
xm = zeros(nM, nX, T);
xm[:,:,1] = rand(nM, nX);
ym = zeros(nM, nY, T-1);

G = create_m_models(A, B, C, D, nM, nX, nY, nU, nR, nT, xm[:,:,1], u0);

for m in 1:nM
    SS[m, 1] = (du = 0., u = u0[m,1], y = 0., x = xm[m,:,1]);
end

tempo = @elapsed let SS = SS, G = G, ym = ym, xm = xm
#    @showprogress 1 "Computing..." for i in 2:T
    for i in 2:T

        obj, S, Sm = sdp_bilevel_master(G);
         # global S;

        push!(Objetivo_bilevel, obj)
        @show i
        for m in 1:G.nM
            SS[m,i] = (du = Sm[m].Du[1,1], u = Sm[m].u[1,2], y = Sm[m].y[1,2], x = Sm[m].x[:,2]);

            xm[m, :, i] = A .* xm[m ,:, i-1] .+ B .* SS[m, i-1].u;
            ym[m, :, i-1] = C .* xm[m, :, i-1] .+ D .* SS[m, i-1].u;

            u0[m,1] = SS[m,i].u;
        end
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

# matwrite("var_plot_bilevel.mat", Dict("u" => upl, "y"=>ypl, "obj_bilevel"=>Objetivo_bilevel))

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
savefig(fn, "m4bilevel.pdf")

writedlm("bilevely.csv", ypl,',')
writedlm("bilevelu.csv", upl,',')
writedlm("bilevelobj.csv", Objetivo_bilevel,',')