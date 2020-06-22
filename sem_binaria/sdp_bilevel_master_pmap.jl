# Julia v1.2.0 & JuMP v0.20
using Distributed, DelimitedFiles
@everywhere using JuMP, Ipopt, Gurobi

# @everywhere include("sdp_setpb_dissertation.jl")
# @everywhere include("sdp_setpb.jl")
# @everywhere include("sdp_setprob_sbai.jl")
# @everywhere include("sdp_setpbcamps.jl")
@everywhere include("sdp_bilevel_sm_pmap.jl")
@everywhere include("sdp_bilevel_s_pmap.jl")
include("read_objective_bilevel.jl")

@everywhere global const GRB_ENV = Gurobi.Env()

if isdir(pwd()*"/bilevel_data") == false
    mkdir("bilevel_data")
end

master_model = Model(optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-4, "max_cpu_time" => 600., "output_file" => "bilevel_data/M$(nM)T$(nT).txt"))
tempo_sub = 0;

function subProblem(x...)
    global grad, Sm, tempo_sub
    s, grad, sm, Sm, tempo = sdp_bilevel_s(x)
    tempo_sub += tempo
    return s
end

custo(x...) = subProblem(x...)

function ∇custo(g, x...)
    for i in 1:length(grad)
        g[i] = grad[i]
    end
    return g
end

# @variable(master_model, Smaster[1:nM[1],1:nT[1] - 1] >= 0);
@variable(master_model, Smaster[1:nM[1],1:nT[1] - 1] >= 0, start = r_max[1,1]/nM[1]);

for r in 1:nR[1]
    for t in 1:nT[1] - 1
        @constraint(master_model, sum(Smaster[m,t] for m in 1:nM[1]) <= r_max[r,t])
    end
end

d = nM[1] * (nT[1] - 1) * nR[1];

JuMP.register(master_model, :custo, d, custo, ∇custo, autodiff = false);
@NLobjective(master_model, Min, custo(Smaster...));

tempo_total = @elapsed JuMP.optimize!.(master_model);
# JuMP.optimize!.(master_model);

objetivo = JuMP.objective_value.(master_model); 

println("The sub time value is: ", tempo_sub)
println("The master time value is: ", tempo_total - tempo_sub)
println("The total time value is: ", tempo_total)

println("\nThe Objective Value is: ", objetivo)

UB = read_objective_bilevel(nM, nT)

arIt = collect(0:length(UB)-1)
artime = range(0,stop=tempo_total,length=length(UB))

open("bilevel_data/M$(nM)T$(nT).csv"; write=true) do f
    write(f, "it,CPU_TIME,UB\n")
    writedlm(f,[arIt artime UB],',')
end;

open("bilevel_data/bilevel.csv", "a") do f
    println(f, "$tempo_total, $objetivo")
end;

# Base.rm("bilevel_data/M$(nM)T$(nT).txt")