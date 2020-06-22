using Distributed
@everywhere using DelimitedFiles
@everywhere using JuMP, Gurobi, LinearAlgebra
# @everywhere include("sdp_setprob_sbai.jl")
# @everywhere include("sdp_setpb_dissertation.jl")
# @everywhere include("sdp_setpb.jl")
# @everywhere include("sdp_benders_sub_pmap.jl")
@everywhere include("sdp_benders_sub.jl")
@everywhere include("sdp_benders_submf.jl")
@everywhere include("sdp_benders_submi.jl")
@everywhere const GRB_ENV = Gurobi.Env()

# Benders Master Problem

max_it = 10000;

global I = Array{NamedTuple}(undef, max_it + 1, 1);   # Indices of infeasibility cuts
global F = Array{NamedTuple}(undef, max_it + 1, 1);   # Indices of feasibility cuts
lF = 0;
lI = 0;

UB = 1e8;
LB = -1e8;      # Cost is nonnegative, quadratic convex
eps = 1e-1;

arLB = [];
arUB = [];
arIt = [];

it = 0;
tempo1 = time()
tempo = @elapsed let it = it, UB = UB, LB = LB, eps = eps, arLB = arLB, arUB = arUB, lF = lF, lI = lI, arIt = arIt, I = I, F = F, tempo1 = tempo1
    while ((UB - LB) / (UB + 1e-20) * 100 >= eps) & (time() - tempo1 < 600) & (it <= max_it)  
        println("------- BEGIN ITERATION ------");
        println("it = ", it);
        it = it + 1;

        # master_model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "mumps", max_iter = 50000))
        master_model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))

        ss_var  = @variable(master_model, ss_var[1:nR[1], 1:nT[1] - 1, 1:nM[1]] >= 0);
        α_B = @variable(master_model, α_B >= LB);  # Cost is nonnegative, quadratic convex

        @objective(master_model, Min, α_B);
     
        for r in 1:nR[1]
            for t in 1:nT[1] - 1
                @constraint(master_model, sum(ss_var[r,t,m] for m in 1:nM[1]) <= r_max[r,t]);
            end
        end

        # Introduce Optimality cuts
        @constraint(master_model, [i in 1:lF], α_B >= F[i].cut_cte + dot(F[i].cut_coef, ss_var));

        # Introduce Feasibility cuts
        @constraint(master_model, [i in 1:lI], I[i].cut_cte + dot(I[i].cut_coef, ss_var) <= 0);
  
        println("-----------  SOLVE MASTER -------------");

        JuMP.optimize!.(master_model);

        println("LB = ", LB)
        println("UB = ", UB)

        status = JuMP.termination_status(master_model)
        
        if status == MOI.OPTIMAL            
            println("----------- IMPROVE LOWER BOUND -------------");
            LB = max(LB, JuMP.objective_value(master_model));
        end

        # ss = max.(value.(ss_var), 0)
        ss = value.(ss_var)

        println("-----------  SOLVE SUBPROBLEMS -------------");
        (sm, Sm, flg_all, flg, cut_benders) = sdp_benders_sub(ss);

        if flg_all == 0 # All feasible problems       
            println("----------- IMPROVE UPPER BOUND -------------");
            s = sum(sm);
            UB = min(UB, s);    
                       
            println("---- Introduce Optimality Cut ----");

            lF = lF + 1;
            F[lF] = (cut_cte = cut_benders[1].cut_cte, cut_coef = cut_benders[1].cut_coef);
        else # At least one problem is infeasible
            println("---- Introduce Feasibility Cut ----");

            lI = lI + 1;
            I[lI] = (cut_cte = cut_benders[1].cut_cte, cut_coef = cut_benders[1].cut_coef);
        end

        push!(arLB, LB);
        push!(arUB, UB);
        push!(arIt, it)
        println("LB = ", LB)
        println("UB = ", UB)
        println("Gap % = ", round((UB - LB) / (UB + 1e-12) * 100, digits = 4))
        global Sm
    end
end
;
println("Solved in: ", tempo, "(s)")


if arUB[end] - arLB[end] < eps
    sol = arLB[end];
else
    sol = arUB[end];
end

gap = (arUB[end] - arLB[end]) / (arUB[end] + 1e-12) * 100;


if isdir(pwd() * "/benders_data") == false
    mkdir("benders_data")
end

open("benders_data/benders.csv", "a") do f
    println(f, "$tempo, $sol, $gap")
end;

open("benders_data/M$(nM)T$(nT).csv"; write = true) do f
    write(f, "it,LB,UB\n")
    writedlm(f, [arIt arLB arUB], ',')
end;