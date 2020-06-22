#############
# OA subproblem

function sdp_OA_submf(ss, δ, m, G)
    
    submf_m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))
    # submf_m = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "mumps", max_iter = 800))

    @variable(submf_m, y[1:G.nY, 1:G.nT]);
    @variable(submf_m, x[1:G.nX, 1:G.nT]);
    @variable(submf_m, u[1:G.nU, 1:G.nT - 1]);
    @variable(submf_m, Du[1:G.nU, 1:G.nT - 1]);

    # -----------------------
    # QP OPtimization
    @objective(submf_m, Min, sum(G.rm[m] * Du[i,t]^2 for i in 1:G.nU, t in 1:G.Num[m] - 1) +
                             sum(G.qm[m] * (y[i,t] - G.wm[m,i,t])^2 for t in G.N1m[m]:G.N2m[m], i in 1:G.nY));

    # --------------------------
    # Constraints

    # Initial states
    @constraint(submf_m, x0_cons[1, i in 1:G.nX], x[i,1] == G.x0[m,i]);

    # Control variation
    @constraint(submf_m, Δu_cons[1, i in 1:G.nU], Du[i,1] == u[i,1] - G.u0[m,i]);

    @constraint(submf_m, Δu_cons2[t in 2:G.nT - 1, i in 1:G.nU], Du[i,t] == u[i,t] - u[i,t - 1]);

    # State Dynamic Equations
    @constraint(submf_m, x_eq_cons[i in 1:G.nX, t in 2:G.nT], x[i,t] == dot(G.Am[m,i,:], x[:,t - 1]) + dot(G.Bm[m,i,:], u[:,t - 1]));
    
    # Output Equations
    @constraint(submf_m, y_eq_cons[i in 1:G.nY, t in 2:G.nT], y[i,t] == dot(G.Cm[m,i,:], x[:,t]) + dot(G.Dm[m,i,:], u[:,t - 1]));

    # Output bounds 
    @constraint(submf_m, y_min_cons[t in G.N1m[m]:G.N2m[m], i in 1:G.nY,], G.y_min[m,i] - y[i,t] <= 0);
    @constraint(submf_m, y_max_cons[t in G.N1m[m]:G.N2m[m], i in 1:G.nY,], y[i,t] - G.y_max[m,i] <= 0);
    
    # bounds on control variation
    @constraint(submf_m, Δu_min_cons[t in 1:G.Num[m] - 1, i in 1:G.nU], G.Du_min[m,i] - Du[i,t] <= 0);
    @constraint(submf_m, Δu_max_cons[t in 1:G.Num[m] - 1, i in 1:G.nU], Du[i,t] - G.Du_max[m,i] <= 0);

    # bounds on control signal
    @constraint(submf_m, u_min_cons[t in 1:G.Num[m] - 1, i in 1:G.nU], (G.u_min[m,i] * δ[m,t]) - u[i,t] <= 0);
    @constraint(submf_m, u_max_cons[t in 1:G.Num[m] - 1, i in 1:G.nU], u[i,t] - (G.u_max[m,i] * δ[m,t]) <= 0);

    # Resource constraints
    @constraint(submf_m, res_cons[t in 1:G.nT - 1, r in 1:G.nR], sum(G.rmu[m,r,j] * u[j,t] for j in 1:G.nU) <= ss[r,t,m]);

    # Solve the problem
    tempo = @elapsed optimize!(submf_m)
    # println(submf_m)

    if has_values(submf_m)
        sm = objective_value(submf_m);

        flg = 0;
        flgc1 = zeros(G.Num[m] - 1)
        flgc2 = zeros(G.Num[m] - 1)

        S = (x = value.(x), y = value.(y), u = value.(u), Du = value.(Du))
        sm_output = (sm = sm, S = S, flg = flg, flgc1 = flgc1, flgc2 = flgc2);

        return sm_output;
    else
        sm_output = (sm = 0.0, S = (x = 0.0, y = 0.0, u = 0.0, Du = 0.0), 
                     flg = 1, flgc1 = ones(G.Num[m] - 1), flgc2 = ones(G.Num[m] - 1));
        return sm_output
    end         
end