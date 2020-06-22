using Distributed
using DelimitedFiles
@everywhere using JuMP, Gurobi, LinearAlgebra
# @everywhere include("sdp_setprob_sbai.jl")
# @everywhere include("sdp_setpb_dissertation.jl")
# @everywhere include("sdp_setpb.jl")
# # @everywhere include("sdp_benders_sub_pmap.jl")
@everywhere include("sdp_benders_sub_bin.jl")
@everywhere include("sdp_benders_submf_bin.jl")
@everywhere include("sdp_benders_submi_bin.jl")
@everywhere const GRB_ENV = Gurobi.Env()

# Benders Master Problem

max_it = 10000;

I = Array{NamedTuple}(undef, max_it + 1);   # Indices of infeasibility cuts
F = Array{NamedTuple}(undef, max_it + 1);   # Indices of feasibility cuts
lF = 0;
lI = 0;

UB = 1e12;
LB = -1e12;     
eps = 1e-2;

arLB = [];
arUB = [];
arIt = [];
global arTempo = [0.0];
global tempo_total = 0.0;
it = 0;
tempos = time()
tempo = @elapsed let it = it, UB = UB, LB = LB, eps = eps, arLB = arLB, arUB = arUB, lF = lF, lI = lI, arIt = arIt, I = I, F = F, tempos = tempos
    while ((UB - LB) / (UB + 1e-12) * 100 > eps) & (tempo_total < 3600.) & (it <= max_it)  
        println("------- BEGIN ITERATION ------");
        println("it = ", it);
        it += 1;

        master_model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0, "Threads" => 1))
        # master_model = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.MSG_OFF))

        @variable(master_model, δ[1:nM[1], 1:nT[1] - 1], Bin);
        @variable(master_model, ss_var[1:nR[1], 1:nT[1] - 1, 1:nM[1]] >= 0);
        @variable(master_model, α_B >= LB);

        @objective(master_model, Min, α_B);

        @constraint(master_model, [r in 1:nR[1], t in 1:nT[1] - 1], sum(ss_var[r,t,m] for m in 1:nM[1]) <= r_max[r,t]);

        # Introduce Optimality cuts
        @constraint(master_model, [i in 1:lF], α_B >= F[i].cut_cte + F[i].cut_cte2 + dot(F[i].cut_coef, ss_var) + dot(F[i].cut_coef2', δ));

        # Introduce Feasibility cuts
        @constraint(master_model, [i in 1:lI], I[i].cut_cte + I[i].cut_cte2 + dot(I[i].cut_coef, ss_var) + dot(I[i].cut_coef2', δ) <= 0);
     
        println("-----------  SOLVE MASTER -------------");

        tempo_master = @elapsed optimize!(master_model);

        # println("LB = ", LB)
        # println("UB = ", UB)

        status = termination_status(master_model)            
        
        if status == MOI.OPTIMAL            
            println("----------- IMPROVE LOWER BOUND -------------");
            LB = max(LB, objective_value(master_model));
            # LB = objective_value(master_model)
        end

        println("-----------  SOLVE SUBPROBLEMS -------------");
        (sm, Sm, flg_all, flg, cut_benders, tempo_sub) = sdp_benders_sub_bin(value.(ss_var), value.(δ));
        # println("FLAG ALL: ", flg_all)

        if flg_all == 0 # All feasible problems       
            println("----------- IMPROVE UPPER BOUND -------------");
            s = sum(sm);
            UB = min(UB, s);    
                                 
            println("---- Introduce Optimality Cut ----");

            lF += 1;
            F[lF] = (cut_cte = cut_benders[1].cut_cte, cut_coef = cut_benders[1].cut_coef, cut_cte2 = cut_benders[1].cut_cte2, cut_coef2 = cut_benders[1].cut_coef2);
        else # At least one problem is infeasible
            println("---- Introduce Feasibility Cut ----");

            lI += 1;
            I[lI] = (cut_cte = cut_benders[1].cut_cte, cut_coef = cut_benders[1].cut_coef, cut_cte2 = cut_benders[1].cut_cte2, cut_coef2 = cut_benders[1].cut_coef2);
        end

        global tempo_total += tempo_master + tempo_sub

        push!(arLB, LB);
        push!(arUB, UB);
        push!(arIt, it)
        push!(arTempo,  round(tempo_total, digits = 4))
        # println("δ: ", value.(δ))
        println("LB = ", LB)
        println("UB = ", UB)
        println("Gap % = ", round((UB - LB) / (UB + 1e-12) * 100, digits = 4))
        # global Sm
    end
end
;
println("Solved in: ", tempo_total, "(s)")

gap = (arUB[end] - arLB[end]) / (arUB[end] + 1e-12) * 100;

if gap < eps
    sol = arLB[end];
else
    sol = arUB[end];
end

gap = round(gap, digits = 4)

if isdir(pwd() * "/benders_data") == false
    mkdir("benders_data")
end

open("benders_data/benders.csv", "a") do f
    println(f, "$tempo_total, $sol, $gap")
end;

open("benders_data/M$(nM)T$(nT).csv"; write = true) do f
    write(f, "it,tempo,LB,UB\n")
    writedlm(f, [arIt arTempo[2:end] arLB arUB], ',')
end;