#############
# OA Infeasibility Subproblem


function sdp_OA_submi(ss, δ, m)

    # fix_γ = 0
    # fix_γ = 1, then γ is fixed at γ_value

    submi_m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))
    # submi_m = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "mumps", max_iter = 200))

    @variable(submi_m, y[1:nY[1], 1:nT[1]]);
    @variable(submi_m, x[1:nX[1], 1:nT[1]]);
    @variable(submi_m, u[1:nU[1], 1:nT[1] - 1]);
    @variable(submi_m, Du[1:nU[1], 1:nT[1] - 1]);

    @variable(submi_m, γ);

    @objective(submi_m, Min, γ);
    
    # --------------------------
    # Constraints

    # Initial states
    @constraint(submi_m, x0_consA[1,i in 1:nX[1]], x[i,1] == x0[m,i]);
    # @constraint(submi_m, x0_consA[i in 1:nX[1]],   x[i,1] - x0[m,i]  <= γ);
    # @constraint(submi_m, x0_consB[i in 1:nX[1]], -(x[i,1] - x0[m,i]) <= γ);

    # Control variation
    @constraint(submi_m, Δu_cons[1,i in 1:nU[1]], Du[i,1] == u[i,1] - u0[m,i]);
    # @constraint(submi_m, Δu_consA[i in 1:nU[1], 1],   Du[i,1] - (u[i,1] - u0[m,i])  <= γ);
    # @constraint(submi_m, Δu_consB[i in 1:nU[1], 1], -(Du[i,1] - (u[i,1] - u0[m,i])) <= γ);

    @constraint(submi_m, Δu_cons2[t in 2:nT[1] - 1, i in 1:nU[1]], Du[i,t] == u[i,t] - u[i,t - 1]);
    # @constraint(submi_m, Δu_cons2A[i in 1:nU[1], t in 2:nT[1] - 1],   Du[i,t] - (u[i,t] - u[i,t - 1])  <= γ);
    # @constraint(submi_m, Δu_cons2B[i in 1:nU[1], t in 2:nT[1] - 1], -(Du[i,t] - (u[i,t] - u[i,t - 1])) <= γ);

    # State Dynamic Equations
    @constraint(submi_m, x_eq_cons[t in 2:nT[1], i in 1:nX[1]], x[i,t] == dot(Am[m,i,:], x[:,t - 1]) + dot(Bm[m,i,:], u[:,t - 1]));
    # @constraint(submi_m, x_eq_consA[i in 1:nX[1], t in 2:nT[1]],    x[i,t] - (dot(Am[m,i,:], x[:,t - 1]) + dot(Bm[m,i,:], u[:,t - 1])) <= γ);
    # @constraint(submi_m, x_eq_consB[i in 1:nX[1], t in 2:nT[1]], - (x[i,t] - (dot(Am[m,i,:], x[:,t - 1]) + dot(Bm[m,i,:], u[:,t - 1]))) <= γ);

    # Output Equations
    @constraint(submi_m, y_eq_cons[t in 2:nT[1], i in 1:nY[1]], y[i,t] == dot(Cm[m,i,:], x[:,t]) + dot(Dm[m,i,:], u[:,t - 1]));
    # @constraint(submi_m, y_eq_consA[i in 1:nY[1], t in 2:nT[1]],    y[i,t] - (dot(Cm[m,i,:], x[:,t]) + dot(Dm[m,i,:], u[:,t - 1])) <= γ);
    # @constraint(submi_m, y_eq_consB[i in 1:nY[1], t in 2:nT[1]], - (y[i,t] - (dot(Cm[m,i,:], x[:,t]) + dot(Dm[m,i,:], u[:,t - 1]))) <= γ);

    # Output bounds 
    @constraint(submi_m, y_min_cons[t in N1m[m]:N2m[m], i in 1:nY[1]], y_min[m,i] - y[i,t] <= γ);
    @constraint(submi_m, y_max_cons[t in N1m[m]:N2m[m], i in 1:nY[1]], y[i,t] - y_max[m,i] <= γ);

    # bounds on control variation
    @constraint(submi_m, Δu_min_cons[t in 1:Num[m] - 1, i in 1:nU[1]], Du_min[m,i] - Du[i,t] <= γ);
    @constraint(submi_m, Δu_max_cons[t in 1:Num[m] - 1, i in 1:nU[1]], Du[i,t] - Du_max[m,i] <= γ);

    # bounds on control signal
    @constraint(submi_m, u_min_cons[t in 1:Num[m] - 1, i in 1:nU[1]], (δ[m,t] * u_min[m,i]) - u[i,t] <= γ);
    @constraint(submi_m, u_max_cons[t in 1:Num[m] - 1, i in 1:nU[1]], u[i,t] - (δ[m,t] * u_max[m,i]) <= γ);

    # Resource constraints
    @constraint(submi_m, res_cons[t in 1:nT[1] - 1, r in 1:nR[1]], sum(rmu[m,r,j] * u[j,t] for j in 1:nU[1]) - ss[r,t,m] <= γ);

    # Solve the problem
    tempo = @elapsed optimize!(submi_m);

    if has_values(submi_m)
        sm = objective_value(submi_m);
        S = (x = value.(x), y = value.(y), u = value.(u), Du = value.(Du))
        flg = 0;
    else    
        return error("Infactivel")
    end  

    ϵ = 1e-8

    flgc1 = zeros(Num[m] - 1)
    flgc2 = zeros(Num[m] - 1)

    let flgc1 = flgc1, flgc2 = flgc2
        for t in 1:Num[m] - 1
            for i in 1:nU[1]
                if abs((δ[m,t] * u_min[m,i]) - value(u[i,t]) - value(γ)) < ϵ
                    flgc1[t] = 0
                else  
                    flgc1[t] = 1
                end
                if abs(value(u[i,t]) - (δ[m,t] * u_max[m,i]) - value(γ)) < ϵ
                    flgc2[t] = 0
                else  
                    flgc2[t] = 1
                end
            end
        end
    end

    sm_output = (sm = sm, S = S, flg = flg, flgc1 = flgc1, flgc2 = flgc2)

    if flg == 0 # 'solved optimally'
        return sm_output;
    else
        println("sdp_OA_submi: Something went wrong!!");
    end
end