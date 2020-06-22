#############
# Benders subproblem BF_m(ss_m)

function sdp_benders_submf(ss, G, m)
    
    sub_model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))
    # sub_model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "mumps", max_iter = 800))

    @variable(sub_model, y[1:G.nY, 1:G.nT]);
    @variable(sub_model, x[1:G.nX, 1:G.nT]);
    @variable(sub_model, u[1:G.nU, 1:G.nT - 1]);
    @variable(sub_model, Du[1:G.nU, 1:G.nT - 1]);

    # -----------------------
    # QP OPtimization
    @objective(sub_model, Min, sum(G.rm[m] * Du[i,t]^2 for i in 1:G.nU, t in 1:G.Num[m] - 1) + 
                                    sum(G.qm[m] * (y[i,t] - G.wm[m,i,t])^2 for i in 1:G.nY, t in G.N1m[m]:G.N2m[m]));
    # --------------------------
    # Constraints

    # Initial states
    @constraint(sub_model, x0_cons[i in 1:G.nX], x[i,1] == G.x0[m,i]);

    # Control variation
    @constraint(sub_model, Δu_cons[i in 1:G.nU, 1], Du[i,1] == u[i,1] - G.u0[m,i]);

    @constraint(sub_model, Δu_cons2[i in 1:G.nU, t in 2:G.nT - 1], Du[i,t] == u[i,t] - u[i,t - 1]);

    # State Dynamic Equations
    @constraint(sub_model, x_eq_cons[i in 1:G.nX, t in 2:G.nT], x[i,t] == dot(G.Am[m,i,:], x[:,t - 1]) + dot(G.Bm[m,i,:], u[:,t - 1]));

    # Output Equations
    @constraint(sub_model, y_eq_cons[i in 1:G.nY, t in 2:G.nT], y[i,t] == dot(G.Cm[m,i,:], x[:,t]) + dot(G.Dm[m,i,:], u[:,t - 1]));
    
    # Output bounds 
    @constraint(sub_model, y_min_cons[i in 1:G.nY, t in G.N1m[m]:G.N2m[m]], G.y_min[m,i] - y[i,t] <= 0.);
    @constraint(sub_model, y_max_cons[i in 1:G.nY, t in G.N1m[m]:G.N2m[m]], y[i,t] - G.y_max[m,i] <= 0.);
    
    # bounds on control variation
    @constraint(sub_model, Δu_min_cons[i in 1:G.nU, t in 1:G.Num[m] - 1], G.Du_min[m,i] - Du[i,t] <= 0.);
    @constraint(sub_model, Δu_max_cons[i in 1:G.nU, t in 1:G.Num[m] - 1], Du[i,t] - G.Du_max[m,i] <= 0.);

    # bounds on control signal
    @constraint(sub_model, u_min_cons[i in 1:G.nU, t in 1:G.Num[m] - 1], G.u_min[m,i] - u[i,t] <= 0.);
    @constraint(sub_model, u_max_cons[i in 1:G.nU, t in 1:G.Num[m] - 1], u[i,t] - G.u_max[m,i] <= 0.);

    # Resource constraints
    @constraint(sub_model, res_cons[r in 1:G.nR, t in 1:G.nT - 1], dot(G.rmu[m,r,:], u[:,t]) - ss[r,t,m] <= 0.);

    # Solve the problem
    tempo = @elapsed JuMP.optimize!(sub_model)

    if has_values(sub_model)
        sm = JuMP.objective_value(sub_model);
    else
        sm_output = (sm = 0.0, S = (x = 0.0, y = 0.0, u = 0.0, Du = 0.0), dual_nu_y1 = 0.0, dual_nu_y2 = 0.0, 
                        dual_nu_u1 = 0.0, dual_nu_u2 = 0.0,
                        dual_nu_du1 = 0.0, dual_nu_du2 = 0.0,
                        dual_mu = 0.0, flg = 1);
        return sm_output
    end        

    S = (x = JuMP.value.(x), y = JuMP.value.(y), u = JuMP.value.(u), Du = JuMP.value.(Du))
    
    if JuMP.termination_status(sub_model) == MOI.OPTIMAL || JuMP.termination_status(sub_model) == MOI.LOCALLY_SOLVED
        flg = 0;
    elseif JuMP.termination_status(sub_model) == MOI.INFEASIBLE
        flg = 1;
    else
        flg = 2;
    end

    dual_nu_y1 = zeros(G.nY, G.N2m[m]);
    dual_nu_y2 = zeros(G.nY, G.N2m[m]);

    dual_nu_u1 = zeros(G.nU, G.Num[m] - 1);
    dual_nu_u2 = zeros(G.nU, G.Num[m] - 1);

    dual_nu_du1 = zeros(G.nU, G.Num[m] - 1);
    dual_nu_du2 = zeros(G.nU, G.Num[m] - 1);

    dual_mu = zeros(G.nR, G.Num[m] - 1);

    ϵ = 1e-1;

    for t in G.N1m[m]:G.N2m[m]
        for i in 1:G.nY                
            if abs(G.y_min[m,i] - JuMP.value(y[i,t]) ) < ϵ
                dual_nu_y1[i,t] = - JuMP.dual(y_min_cons[i,t]);
            else
                dual_nu_y1[i,t] = 0.0;
            end
        end
    end    
    
    for t in G.N1m[m]:G.N2m[m]
        for i = 1:G.nY                
            if abs(JuMP.value(y[i,t]) - G.y_max[m,i] ) < ϵ
                dual_nu_y2[i,t] = - JuMP.dual(y_max_cons[i,t]);
            else
                dual_nu_y2[i,t] = 0.0;
            end
        end
    end

    for t in 1:G.Num[m] - 1
        for i in 1:G.nU                
            if abs(G.Du_min[m,i] - JuMP.value(Du[i,t]) ) < ϵ
                dual_nu_du1[i,t] = - JuMP.dual(Δu_min_cons[i,t]);
            else
                dual_nu_du1[i,t] = 0.0;
            end
        end
    end


    for t in 1:G.Num[m] - 1
        for i in 1:G.nU                
            if abs(JuMP.value(Du[i,t]) - G.Du_max[m,i] ) < ϵ
                dual_nu_du2[i,t] = - JuMP.dual(Δu_max_cons[i,t]);
            else
                dual_nu_du2[i,t] = 0.0;
            end
        end
    end

    for t in 1:G.Num[m] - 1
        for i in 1:G.nU                
            if abs(G.u_min[m,i] - JuMP.value(u[i,t]) ) < ϵ
                dual_nu_u1[i,t] = - JuMP.dual(u_min_cons[i,t]);
            else
                dual_nu_u1[i,t] = 0.0;
            end
        end
    end

    for t in 1:G.Num[m] - 1
        for i in 1:G.nU                
            if abs(JuMP.value(u[i,t]) -  G.u_max[m,i] ) < ϵ
                dual_nu_u2[i,t] = - JuMP.dual(u_max_cons[i,t]);
            else
                dual_nu_u2[i,t] = 0.0;
            end
        end
    end

    for t in 1:G.Num[m] - 1
        for r in 1:G.nR
            eq = dot(G.rmu[m,r,:], JuMP.value.(u[:,t]))
            if abs(eq - ss[r,t,m] ) < ϵ
                dual_mu[r,t] = - JuMP.dual(res_cons[r,t]);
            else
                dual_mu[r,t] = 0.0;
            end
        end
    end
    
    sm_output = (sm = sm, S = S, dual_nu_y1 = dual_nu_y1, dual_nu_y2 = dual_nu_y2, 
                dual_nu_u1 = dual_nu_u1, dual_nu_u2 = dual_nu_u2,
                dual_nu_du1 = dual_nu_du1, dual_nu_du2 = dual_nu_du2,
                dual_mu = dual_mu, flg = flg, tempo = tempo);

    return sm_output;
end