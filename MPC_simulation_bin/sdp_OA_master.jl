function sdp_OA_master(G)
# Outer Approximation Master Problem
max_it = 10000;

I = Array{NamedTuple}(undef, max_it + 1);   # Indices of infeasibility cuts
F = Array{NamedTuple}(undef, max_it + 1);

global lF = 0;
global lI = 0;

UB = 1e12;
LB = -1e12;    
eps = 1e-6;

arLB = [];
arUB = [];
arIt = [];
arsm = [];
arSOL = []

it = 0;
tempos = time()
tempo = @elapsed let it = it, UB = UB, LB = LB, eps = eps, arLB = arLB, arUB = arUB, lF = lF, lI = lI, arIt = arIt, I = I, F = F, tempos = tempos, arsm = arsm, arSOL = arSOL
    while ((UB - LB) > eps) & (time()-tempos < 300.) & (it <= max_it)  
        # println("------- BEGIN ITERATION ------");
        # println("it = ", it);
        it += 1;

        if it == 1
            ss = zeros(G.nR, G.nT - 1, G.nM);
            delta = zeros(G.nM, G.nT - 1);

            for r in 1:G.nR
                for t in 1:G.nT - 1
                    for m in 1:G.nM
                        ss[r,t,m] =  G.r_max[r,t] / G.nM;
                    end
                end
            end
        else
            ss = ss_save;
            delta = d_save
        end

        # println("-----------  SOLVE SUBPROBLEMS -------------");
        (sm, Sm, flg_all, flgc1, flgc2, cut_OA) = sdp_OA_sub(ss, delta, G);
        push!(arsm, Sm)

        if flg_all == 0 # All feasible problems       
            # println("----------- IMPROVE UPPER BOUND -------------");
            s = sum(sm);
            UB = min(UB, s);    
                                 
            # println("---- Introduce Optimality Cut ----");

            global lF += 1;
            F[lF] = (smh = cut_OA[1].smh, Smh = cut_OA[1].Smh);
        else # At least one problem is infeasible
            # println("---- Introduce Feasibility Cut ----");

            global lI += 1;
            I[lI] = (smh = cut_OA[1].smh, Smh = cut_OA[1].Smh);
        end

        master_model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))

        @variable(master_model, y[1:G.nM,  1:G.nY, 1:G.nT]);
        @variable(master_model, x[1:G.nM, 1:G.nX, 1:G.nT]);
        @variable(master_model, u[1:G.nM, 1:G.nU, 1:G.nT - 1]);
        @variable(master_model, Du[1:G.nM, 1:G.nU, 1:G.nT - 1]);

        @variable(master_model, δ[1:G.nM, 1:G.nT - 1], Bin);
        @variable(master_model, ss_var[1:G.nR, 1:G.nT - 1, 1:G.nM] >= 0);
        @variable(master_model, α_OA >= LB);  

        @objective(master_model, Min, α_OA);

        for r in 1:G.nR
            for t in 1:G.nT - 1
                @constraint(master_model, sum(ss_var[r,t,m] for m in 1:G.nM) <= G.r_max[r,t]);
            end
        end
        
        for k in 1:lF
            @constraint(master_model, α_OA >= F[k].smh + 
                                              sum(G.rm[m] * 2 * F[k].Smh[m].Du[i,t] * (Du[m,i,t] - F[k].Smh[m].Du[i,t]) 
                                                                        for m in 1:G.nM, i in 1:G.nU, t in 1:G.Num[m] - 1) + 
                                              sum(G.qm[m] * 2 * (F[k].Smh[m].y[i,t] - G.wm[m,i,t]) * (y[m,i,t] - F[k].Smh[m].y[i,t])
                                                                        for m in 1:G.nM, t in G.N1m[m]:G.N2m[m], i in 1:G.nY));
        end
        
        # Initial states
        @constraint(master_model, x0_cons[m in 1:G.nM, i in 1:G.nX], x[m,i,1] == G.x0[m,i]);

        # Control variation
        @constraint(master_model, Δu_cons[m in 1:G.nM, i in 1:G.nU, 1], Du[m,i,1] == u[m,i,1] - G.u0[m,i]);

        @constraint(master_model, Δu_cons2[m in 1:G.nM, i in 1:G.nU, t in 2:G.nT - 1], Du[m,i,t] == u[m,i,t] - u[m,i,t - 1]);

        # State Dynamic Equations
        @constraint(master_model, x_eq_cons[m in 1:G.nM, i in 1:G.nX, t in 2:G.nT], x[m,i,t] == dot(G.Am[m,i,:], x[m,:,t - 1]) + dot(G.Bm[m,i,:], u[m,:,t - 1]));

        # Output Equations
        @constraint(master_model, y_eq_cons[m in 1:G.nM, i in 1:G.nY, t in 2:G.nT], y[m,i,t] == dot(G.Cm[m,i,:], x[m,:,t]) + dot(G.Dm[m,i,:], u[m,:,t - 1]));

        # Output bounds 
        @constraint(master_model, y_min_cons[m in 1:G.nM, i in 1:G.nY, t in G.N1m[m]:G.N2m[m]], G.y_min[m,i] - y[m,i,t] <= 0);
        @constraint(master_model, y_max_cons[m in 1:G.nM, i in 1:G.nY, t in G.N1m[m]:G.N2m[m]], y[m,i,t] - G.y_max[m,i] <= 0);

        # bounds on control variation
        @constraint(master_model, Δu_min_cons[m in 1:G.nM, i in 1:G.nU, t in 1:G.Num[m] - 1], G.Du_min[m,i] - Du[m,i,t] <= 0);
        @constraint(master_model, Δu_max_cons[m in 1:G.nM, i in 1:G.nU, t in 1:G.Num[m] - 1], Du[m,i,t] - G.Du_max[m,i] <= 0);
        
        # bounds on control signal
        
        for k in 1:lI
            for m in 1:G.nM
                for t in 1:G.Num[m] - 1
                    for i in 1:G.nU
                        if flgc1[m][t] == 0
                            @constraint(master_model, (δ[m,t] * G.u_min[m,i]) - I[k].Smh[m].u[i,t] - (u[m,i,t] - I[k].Smh[m].u[i,t]) <= 0);
                            @constraint(master_model, I[k].Smh[m].u[i,t] + (u[m,i,t] - I[k].Smh[m].u[i,t] ) - (G.u_max[m,i] * δ[m,t]) <= 0);
                        end
                    end
                end
            end
        end

        if flg_all == 0      
            @constraint(master_model, u_min_cons[m in 1:G.nM, i in 1:G.nU, t in 1:G.Num[m] - 1], (δ[m,t] * G.u_min[m,i]) - u[m,i,t] <= 0);
            @constraint(master_model, u_max_cons[m in 1:G.nM, i in 1:G.nU, t in 1:G.Num[m] - 1], u[m,i,t] - (δ[m,t] * G.u_max[m,i]) <= 0);
        end

        # Resource constraints
        @constraint(master_model, res_cons[m in 1:G.nM, r in 1:G.nR, t in 1:G.nT - 1], dot(G.rmu[m,r,:], u[m,:,t]) <= ss_var[r,t,m]);
       
        # println("-----------  SOLVE MASTER -------------");
        # println(master_model)

        tempo_master = @elapsed optimize!.(master_model);

        # println("LB = ", LB)
        # println("UB = ", UB)

        status = termination_status(master_model)

        if status == MOI.OPTIMAL
            # println("----------- IMPROVE LOWER BOUND -------------");
            LB = max(LB, objective_value(master_model));
        end

        push!(arSOL, LB)

        global ss_save = value.(ss_var);
        global d_save = value.(δ)

        # global tempo_total += tempo_master + tempo_sub

        push!(arLB, LB);
        push!(arUB, UB);
        push!(arIt, it)
        # push!(arTempo,  round(tempo_total, digits = 4))
        # println("δ: ",  value.(δ))
        # println("LB = ", LB)
        # println("UB = ", UB)
        # println("Gap % = ", round((UB - LB) / (UB + 1e-12) * 100, digits = 4))
        # global Sm
    end
end
;
    return arSOL[end], arsm[end]
end
# println("Solved in: ", tempo_total, "(s)")

# gap = (arUB[end] - arLB[end]) / (arUB[end] + 1e-12) * 100;

# if gap < eps
#     sol = arLB[end];
# else
#     sol = arUB[end];
# end

# gap = round(gap, digits = 4)
# tempo_total = round(tempo_total, digits = 4)

# if isdir(pwd() * "/OA_data") == false
#     mkdir("OA_data")
# end

# open("OA_data/OA.csv", "a") do f
#     println(f, "$tempo_total, $sol, $gap")
# end;

# open("OA_data/M$(nM)T$(nT).csv"; write = true) do f
#     write(f, "it,tempo,LB,UB\n")
#     writedlm(f, [arIt arTempo[2:end] arLB arUB], ',')
# end;
