function sdp_benders_master(G)
    
    # Benders Master Problem

    max_it = 3500;

    I = Array{NamedTuple}(undef, max_it+1, 1);   # Indices of infeasibility cuts
    F = Array{NamedTuple}(undef, max_it+1, 1);   # Indices of feasibility cuts
    lF = 0;
    lI = 0;

    UB = 1e8;
    LB = -1e8;      # Cost is nonnegative, quadratic convex
    eps = 1e-6;

    arLB = [];
    arUB = [];
    arSOL = [];
    arsm = [];
    time1 = time()
    it = 0;
    tempo = @elapsed let it = it, UB = UB, LB = LB, eps = eps, arLB = arLB, arUB = arUB, lF = lF, lI = lI, arSOL = arSOL, G = G, arsm = arsm, time1 = time1
            while (it <= max_it) & ((UB - LB) > eps) & (time() - time1 < 60)
                # println("------- BEGIN ITERATION ------");
                # println("it = ", it);
                it = it+1;

                master_model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))

                ss_var  = @variable(master_model, ss_var[1:G.nR, 1:G.nT - 1, 1:G.nM] >= 0);
                α_B = @variable(master_model, α_B >= 0);  # Cost is nonnegative, quadratic convex

                @objective(master_model, Min, α_B);
            
                for r in 1:G.nR
                    for t in 1:G.nT - 1
                        @constraint(master_model, sum(ss_var[r,t,m] for m in 1:G.nM) <= G.r_max[r,t]);
                    end
                end

                # Introduce Optimality cuts
                @constraint(master_model, [i in 1:lF], α_B >= F[i].cut_cte + dot(F[i].cut_coef, ss_var));

                # Introduce Feasibility cuts
                @constraint(master_model, [i in 1:lI], I[i].cut_cte + dot(I[i].cut_coef, ss_var) <= 0);
        
                # println("-----------  SOLVE MASTER -------------");

                JuMP.optimize!.(master_model);

                # println("LB = ", LB)
                # println("UB = ", UB)

                status = JuMP.termination_status(master_model)
                
                if status == MOI.OPTIMAL            
                    # println("----------- IMPROVE LOWER BOUND -------------");
                    LB = max(LB, JuMP.objective_value(master_model));
                end
                push!(arSOL, LB)
                ss = value.(ss_var)

                # ss = max.(value.(ss_var), 0)
                # if it == 1
                #     ss = (G.r_max[1,1]/G.nM) * ones(G.nR, G.nT - 1, G.nM)
                # else
                #     ss = value.(ss_var)
                # end
                # println("-----------  SOLVE SUBPROBLEMS -------------");
                
                (sm, Sm, flg_all, flg, cut_benders) = sdp_benders_sub(ss, G);
                
                push!(arsm, Sm)

                if flg_all == 0 # All feasible problems       
                    # println("----------- IMPROVE UPPER BOUND -------------");
                    s = sum(sm);
                    UB = min(UB, s);    
                            
                    # println("---- Introduce Optimality Cut ----");

                    lF = lF + 1;
                    F[lF] = (cut_cte = cut_benders[1].cut_cte, cut_coef = cut_benders[1].cut_coef);
                else # At least one problem is infeasible
                    # println("---- Introduce Feasibility Cut ----");

                    lI = lI + 1;
                    I[lI] = (cut_cte = cut_benders[1].cut_cte, cut_coef = cut_benders[1].cut_coef);
                end

                push!(arLB, LB);
                push!(arUB, UB);
                # println("LB = ", LB)
                # println("UB = ", UB)
                # println("Gap % = ", round((UB - LB) / (UB + 1e-12) * 100, digits = 4))
                
            end
        end
    ;
    # println("Solved in: ", tempo, "(s)")

    return arSOL[end], arsm[end];
end