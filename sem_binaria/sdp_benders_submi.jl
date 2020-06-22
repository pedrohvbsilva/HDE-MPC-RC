#############
# Benders Infeasibility Subproblem BI_m(ss_m)


function sdp_benders_submi(ss, m)

    # fix_γ = 0
    # fix_γ = 1, then γ is fixed at γ_value

    sub_model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))
    # sub_model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "mumps", max_iter = 200))

    y = @variable(sub_model, y[1:nY[1], 1:nT[1]]);
    x = @variable(sub_model, x[1:nX[1], 1:nT[1]]);
    u = @variable(sub_model, u[1:nU[1], 1:nT[1] - 1]);
    Du = @variable(sub_model, Du[1:nU[1], 1:nT[1] - 1]);

    @variable(sub_model, γ);

    @objective(sub_model, Min, γ);

    # --------------------------
    # Constraints

    # if fix_γ == 0
    #     @constraint(sub_model, γ >= 0);
    # else
    #     @constraint(sub_model, γ >= γ_value);
    # end    #-----------------------

    # Initial states
    @constraint(sub_model, x0_cons[i in 1:nX[1]], x[i,1] == x0[m,i]);

    # Control variation
    @constraint(sub_model, Δu_cons[i in 1:nU[1], 1], Du[i,1] == u[i,1] - u0[m,i]);

    @constraint(sub_model, Δu_cons2[i in 1:nU[1], t in 2:nT[1] - 1], Du[i,t] == u[i,t] - u[i,t - 1]);

    # State Dynamic Equations
    @constraint(sub_model, x_eq_cons[i in 1:nX[1], t in 2:nT[1]], x[i,t] == dot(Am[m,i,:], x[:,t - 1]) + dot(Bm[m,i,:], u[:,t - 1]));

    # Output Equations
    @constraint(sub_model, y_eq_cons[i in 1:nY[1], t in 2:nT[1]], y[i,t] == dot(Cm[m,i,:], x[:,t]) + dot(Dm[m,i,:], u[:,t - 1]));

    # Output bounds 
    @constraint(sub_model, y_min_cons[i in 1:nY[1], t in N1m[m]:N2m[m]], y_min[m,i] - y[i,t] - γ <= 0.);
    @constraint(sub_model, y_max_cons[i in 1:nY[1], t in N1m[m]:N2m[m]], y[i,t] - y_max[m,i] - γ <= 0.);

    # bounds on control variation
    @constraint(sub_model, Δu_min_cons[i in 1:nU[1], t in 1:Num[m] - 1], Du_min[m,i] - Du[i,t] - γ <= 0.);
    @constraint(sub_model, Δu_max_cons[i in 1:nU[1], t in 1:Num[m] - 1], Du[i,t] - Du_max[m,i] - γ <= 0.);

    # bounds on control signal
    @constraint(sub_model, u_min_cons[i in 1:nU[1], t in 1:Num[m] - 1],  u_min[m,i] - u[i,t] - γ <= 0.);
    @constraint(sub_model, u_max_cons[i in 1:nU[1], t in 1:Num[m] - 1],  u[i,t] -  u_max[m,i] - γ <= 0.);

    # Resource constraints
    @constraint(sub_model, res_cons[r in 1:nR[1], t in 1:Num[m] - 1], dot(rmu[m,r,:], u[:,t]) - ss[r,t,m] - γ <= 0.);

    # Solve the problem
    tempo = @elapsed JuMP.optimize!(sub_model);

    if has_values(sub_model)
        sm = JuMP.objective_value(sub_model);
    else    
        return error("Infactivel")
    end  

    S = (x = JuMP.value.(x), y = JuMP.value.(y), u = JuMP.value.(u), Du = JuMP.value.(Du))

    if JuMP.termination_status(sub_model) == MOI.OPTIMAL || JuMP.termination_status(sub_model) == MOI.LOCALLY_SOLVED
        flg = 0;
    elseif JuMP.termination_status(sub_model) == MOI.INFEASIBLE
        flg = 1;
    elseif JuMP.termination_status(sub_model) == MOI.INFEASIBLE_OR_UNBOUNDED
        flg = 2;
    end

    #= ======================================== =#
    # Dual variables for INEQUALITIES
    #= ======================================== =#

    dual_nu_y1 = zeros(nY[1], N2m[m]);
    dual_nu_y2 = zeros(nY[1], N2m[m]);

    dual_nu_u1 = zeros(nU[1], Num[m] - 1);
    dual_nu_u2 = zeros(nU[1], Num[m] - 1);

    dual_nu_du1 = zeros(nU[1], Num[m] - 1);
    dual_nu_du2 = zeros(nU[1], Num[m] - 1);

    dual_mu = zeros(nR[1], Num[m] - 1);

    ϵ = 1e-8;

    for t in N1m[m]:N2m[m]
        for i in 1:nY[1]                
            if abs(y_min[m,i] - JuMP.value(y[i,t]) - JuMP.value(γ)) < ϵ
                dual_nu_y1[i,t] = - JuMP.dual(y_min_cons[i,t]);
            else
                dual_nu_y1[i,t] = 0.0;
            end
        end
    end    
    
    for t in N1m[m]:N2m[m]
        for i = 1:nY[1]                
            if abs(JuMP.value(y[i,t]) - y_max[m,i] - JuMP.value(γ)) < ϵ
                dual_nu_y2[i,t] = - JuMP.dual(y_max_cons[i,t]);
            else
                dual_nu_y2[i,t] = 0.0;
            end
        end
    end

    for t in 1:Num[m] - 1
        for i in 1:nU[1]                
            if abs(Du_min[m,i] - JuMP.value(Du[i,t]) - JuMP.value(γ)) < ϵ
                dual_nu_du1[i,t] = - JuMP.dual(Δu_min_cons[i,t]);
            else
                dual_nu_du1[i,t] = 0.0;
            end
        end
    end


    for t in 1:Num[m] - 1
        for i in 1:nU[1]                
            if abs(JuMP.value(Du[i,t]) - Du_max[m,i] - JuMP.value(γ)) < ϵ
                dual_nu_du2[i,t] = - JuMP.dual(Δu_max_cons[i,t]);
            else
                dual_nu_du2[i,t] = 0.0;
            end
        end
    end

    for t in 1:Num[m] - 1
        for i in 1:nU[1]                
            if abs(u_min[m,i] - JuMP.value(u[i,t]) - JuMP.value(γ)) < ϵ
                dual_nu_u1[i,t] = - JuMP.dual(u_min_cons[i,t]);
            else
                dual_nu_u1[i,t] = 0.0;
            end
        end
    end

    for t in 1:Num[m] - 1
        for i in 1:nU[1]                
            if abs(JuMP.value(u[i,t]) -  u_max[m,i] - JuMP.value(γ)) < ϵ
                dual_nu_u2[i,t] = - JuMP.dual(u_max_cons[i,t]);
            else
                dual_nu_u2[i,t] = 0.0;
            end
        end
    end

    for t in 1:Num[m] - 1
        for r in 1:nR[1]
            eq = dot(rmu[m,r,:], JuMP.value.(u[:,t]))
            if abs(eq - ss[r,t,m] - JuMP.value(γ)) < ϵ
                dual_mu[r,t] = - JuMP.dual(res_cons[r,t]);
            else
                dual_mu[r,t] = 0.0;
            end
        end
    end


    sm_output = (sm = sm, S = S, dual_nu_y1 = dual_nu_y1, dual_nu_y2 = dual_nu_y2, 
                dual_nu_u1 = dual_nu_u1, dual_nu_u2 = dual_nu_u2,
                dual_nu_du1 = dual_nu_du1, dual_nu_du2 = dual_nu_du2,
                dual_mu = dual_mu, flg = flg, tempo = tempo)

    if flg == 0 # 'solved optimally'
        return sm_output;
    else
        println("sdp_benders_submi: Something went wrong!!");
    end
end