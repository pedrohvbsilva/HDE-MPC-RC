using Distributed
using DelimitedFiles
@everywhere using JuMP, Gurobi, LinearAlgebra
# @everywhere include("sdp_setprob_sbai.jl")
# @everywhere include("sdp_setpb_dissertation.jl")
# @everywhere include("sdp_setpb.jl")
# # @everywhere include("sdp_benders_sub_pmap.jl")
@everywhere include("sdp_OA_sub.jl")
@everywhere include("sdp_OA_submf.jl")
@everywhere include("sdp_OA_submi.jl")
@everywhere const GRB_ENV = Gurobi.Env()

# Outer Approximation Master Problem
max_it = 10000;

I = Array{NamedTuple}(undef, max_it + 1);   # Indices of infeasibility cuts
F = Array{NamedTuple}(undef, max_it + 1);

global lF = 0;
global lI = 0;

UB = 1e12;
LB = -1e12;    
eps = 1e-2;

arLB = [];
arUB = [];
arIt = [];
global arTempo = [0.0];
global tempo_total = 0.0;

it = 0;
tempo1 = time()
tempo = @elapsed let it = it, UB = UB, LB = LB, eps = eps, arLB = arLB, arUB = arUB,  arIt = arIt, I = I, F = F, tempo1 = tempo1
    while ((UB - LB) / (UB + 1e-12) * 100 > eps) & (tempo_total < 3600.) & (it <= max_it)  
        println("------- BEGIN ITERATION ------");
        println("it = ", it);
        it += 1;

        if it == 1
            ss = zeros(nR[1], nT[1] - 1, nM[1]);
            delta = zeros(nM[1], nT[1] - 1);

            for r in 1:nR[1]
                for t in 1:nT[1] - 1
                    for m in 1:nM[1]
                        ss[r,t,m] = r_max[r,t] / nM[1];
                    end
                end
            end
        else
            ss = ss_save;
            delta = d_save
        end

        println("-----------  SOLVE SUBPROBLEMS -------------");
        (sm, Sm, flg_all, flgc1, flgc2, cut_OA, tempo_sub) = sdp_OA_sub(ss, delta);

        if flg_all == 0 # All feasible problems       
            println("----------- IMPROVE UPPER BOUND -------------");
            s = sum(sm);
            UB = min(UB, s);    
                                 
            println("---- Introduce Optimality Cut ----");

            global lF += 1;
            F[lF] = (smh = cut_OA[1].smh, Smh = cut_OA[1].Smh);
        else # At least one problem is infeasible
            println("---- Introduce Feasibility Cut ----");

            global lI += 1;
            I[lI] = (smh = cut_OA[1].smh, Smh = cut_OA[1].Smh);
        end

        master_model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0, "Threads" => 1))

        @variable(master_model, y[1:nM[1],  1:nY[1], 1:nT[1]]);
        @variable(master_model, x[1:nM[1], 1:nX[1], 1:nT[1]]);
        @variable(master_model, u[1:nM[1], 1:nU[1], 1:nT[1] - 1]);
        @variable(master_model, Du[1:nM[1], 1:nU[1], 1:nT[1] - 1]);

        @variable(master_model, δ[1:nM[1], 1:nT[1] - 1], Bin);
        @variable(master_model, ss_var[1:nR[1], 1:nT[1] - 1, 1:nM[1]] >= 0);
        @variable(master_model, α_OA >= LB);  

        @objective(master_model, Min, α_OA);

        for r in 1:nR[1]
            for t in 1:nT[1] - 1
                @constraint(master_model, sum(ss_var[r,t,m] for m in 1:nM[1]) <= r_max[r,t]);
            end
        end
        
        for k in 1:lF
            @constraint(master_model, α_OA >= F[k].smh + 
                                              sum(rm[m] * 2 * F[k].Smh[m].Du[i,t] * (Du[m,i,t] - F[k].Smh[m].Du[i,t]) 
                                                                        for m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1) + 
                                              sum(qm[m] * 2 * (F[k].Smh[m].y[i,t] - wm[m,i,t]) * (y[m,i,t] - F[k].Smh[m].y[i,t])
                                                                        for m in 1:nM[1], t in N1m[m]:N2m[m], i in 1:nY[1]));
        end
        
        # Initial states
        @constraint(master_model, x0_cons[m in 1:nM[1], i in 1:nX[1]], x[m,i,1] == x0[m,i]);

        # Control variation
        @constraint(master_model, Δu_cons[m in 1:nM[1], i in 1:nU[1], 1], Du[m,i,1] == u[m,i,1] - u0[m,i]);

        @constraint(master_model, Δu_cons2[m in 1:nM[1], i in 1:nU[1], t in 2:nT[1] - 1], Du[m,i,t] == u[m,i,t] - u[m,i,t - 1]);

        # State Dynamic Equations
        @constraint(master_model, x_eq_cons[m in 1:nM[1], i in 1:nX[1], t in 2:nT[1]], x[m,i,t] == dot(Am[m,i,:], x[m,:,t - 1]) + dot(Bm[m,i,:], u[m,:,t - 1]));

        # Output Equations
        @constraint(master_model, y_eq_cons[m in 1:nM[1], i in 1:nY[1], t in 2:nT[1]], y[m,i,t] == dot(Cm[m,i,:], x[m,:,t]) + dot(Dm[m,i,:], u[m,:,t - 1]));

        # Output bounds 
        @constraint(master_model, y_min_cons[m in 1:nM[1], i in 1:nY[1], t in N1m[m]:N2m[m]], y_min[m,i] - y[m,i,t] <= 0);
        @constraint(master_model, y_max_cons[m in 1:nM[1], i in 1:nY[1], t in N1m[m]:N2m[m]], y[m,i,t] - y_max[m,i] <= 0);

        # bounds on control variation
        @constraint(master_model, Δu_min_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], Du_min[m,i] - Du[m,i,t] <= 0);
        @constraint(master_model, Δu_max_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], Du[m,i,t] - Du_max[m,i] <= 0);
        
        # bounds on control signal
        
        for k in 1:lI
            for m in 1:nM[1]
                for t in 1:Num[m] - 1
                    for i in 1:nU[1]
                        if flgc1[m][t] == 0
                            @constraint(master_model, (δ[m,t] * u_min[m,i]) - I[k].Smh[m].u[i,t] - (u[m,i,t] - I[k].Smh[m].u[i,t]) <= 0);
                            @constraint(master_model, I[k].Smh[m].u[i,t] + (u[m,i,t] - I[k].Smh[m].u[i,t] ) - (u_max[m,i] * δ[m,t]) <= 0);
                        end
                    end
                end
            end
        end

        if flg_all == 0      
            @constraint(master_model, u_min_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], (δ[m,t] * u_min[m,i]) - u[m,i,t] <= 0);
            @constraint(master_model, u_max_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], u[m,i,t] - (δ[m,t] * u_max[m,i]) <= 0);
        end

        # Resource constraints
        @constraint(master_model, res_cons[m in 1:nM[1], r in 1:nR[1], t in 1:nT[1] - 1], dot(rmu[m,r,:], u[m,:,t]) <= ss_var[r,t,m]);
       
        println("-----------  SOLVE MASTER -------------");
        # println(master_model)

        tempo_master = @elapsed optimize!.(master_model);

        println("LB = ", LB)
        println("UB = ", UB)

        status = termination_status(master_model)

        if status == MOI.OPTIMAL
            println("----------- IMPROVE LOWER BOUND -------------");
            LB = max(LB, objective_value(master_model));
        end
        
        global ss_save = value.(ss_var);
        global d_save = value.(δ)

        global tempo_total += tempo_master + tempo_sub

        push!(arLB, LB);
        push!(arUB, UB);
        push!(arIt, it)
        push!(arTempo,  round(tempo_total, digits = 4))
        println("δ: ",  value.(δ))
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
tempo_total = round(tempo_total, digits = 4)

if isdir(pwd() * "/OA_data") == false
    mkdir("OA_data")
end

open("OA_data/OA.csv", "a") do f
    println(f, "$tempo_total, $sol, $gap")
end;

open("OA_data/M$(nM)T$(nT).csv"; write = true) do f
    write(f, "it,tempo,LB,UB\n")
    writedlm(f, [arIt arTempo[2:end] arLB arUB], ',')
end;
