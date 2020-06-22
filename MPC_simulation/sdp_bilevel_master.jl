# Julia v1.1.0 & JuMP v0.19
# using Distributed
# @everywhere using JuMP, Ipopt, Gurobi
#
# @everywhere include("sdp_setprob_sbai.jl")
# @everywhere include("sdp_bilevel_sm_pmap.jl")
# @everywhere include("sdp_bilevel_s_pmap.jl")

function sdp_bilevel_master(G)

    master_model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0,"tol" => 1e-6, "max_iter" => 200))

    tempo_sub = 0;

    function subProblem(x...)
        global grad, Sm, tempo_sub
        tempo = @elapsed s, grad, sm, Sm = sdp_bilevel_s(x, G)
        # tempo_sub += tempo
        return s
    end

    custo(x...) = subProblem(x...)

    function ∇custo(g, x...)
        for i in 1:length(grad)
            g[i] = grad[i]
        end
        return g
    end

    # @variable(master_model, Smaster[1:G.nM, 1:G.nT - 1] >= 0, start = G.r_max[1,1]/G.nM);
    @variable(master_model, Smaster[1:G.nM, 1:G.nT - 1] >= 0);

    for r in 1:G.nR
        for t in 1:G.nT - 1
            @constraint(master_model, sum(Smaster[m,t] for m in 1:G.nM) <= G.r_max[r,t])
        end
    end

    d = G.nM * (G.nT - 1);

    JuMP.register(master_model, :custo, d, custo, ∇custo, autodiff = false)
    @NLobjective(master_model, Min, custo(Smaster...))

    tempo_total = @elapsed JuMP.optimize!.(master_model);

    objetivo = JuMP.objective_value.(master_model)

    # println("The sub time value is: ", tempo_sub)
    # println("The master time value is: ", tempo_total - tempo_sub)
    # println("The total time value is: ", tempo_total)

    # println("\nThe Objective Value is: ", objetivo)
    
    return objetivo, Smaster, Sm
end
