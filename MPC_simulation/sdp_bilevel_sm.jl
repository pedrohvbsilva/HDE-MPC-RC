#Bilevel Optimization: 2nd level problem solver for unit "m"
function sdp_bilevel_sm(e, m, G)

    # sub_model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "ma97", tol = 1e-6, max_iter = 200));
    sub_model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0, "OptimalityTol" => 1e-6, "IterationLimit" => 200))
    # sub_model = Model(with_optimizer(CPLEX.Optimizer))

    y = @variable(sub_model, y[1:G.nY, 1:G.nT]);
    x = @variable(sub_model, x[1:G.nX, 1:G.nT]);
    u = @variable(sub_model, u[1:G.nU, 1:G.nT - 1]);
    Du = @variable(sub_model, Du[1:G.nU, 1:G.nT - 1]);

    rsm = @variable(sub_model, rsm[1:G.nR, 1:G.Num[m] - 1] >= 0);

    # QP OPtimization
    @objective(sub_model, Min, sum(G.rm[m] * Du[i,t]^2 for i in 1:G.nU for t in 1:G.Num[m] - 1) +
                           sum(G.qm[m] * (y[i,t] - G.wm[m,i,t])^2 for t in G.N1m[m]:G.N2m[m] for i in 1:G.nY));

    # Constraints
    Constraints = [];
    # Initial states
    for i in 1:G.nX
        push!(Constraints, @constraint(sub_model, x[i,1] == G.x0[m,i]));
    end
    
    # for t in 2:G.nT
    #     for i in 1:G.nX
    #         push!(Constraints, @constraint(sub_model, x[i,t] <= 0.99))
    #     end
    # end
    
    # Control variation
    for i in 1:G.nU
        push!(Constraints, @constraint(sub_model, Du[i,1] == u[i,1] - G.u0[m,i]));
    end

    for t in 2:G.nT - 1
        for i = 1:G.nU
            push!(Constraints, @constraint(sub_model, Du[i,t] == u[i,t] - u[i,t - 1]));
        end
    end

    # State Dynamic Equations
    for t in 2:G.nT
        for i in 1:G.nX
            eq = 0;
            for j in 1:G.nX
                eq = eq + G.Am[m,i,j] * x[j,t - 1];
            end
            for j in 1:G.nU
                eq = eq + G.Bm[m,i,j] * u[j,t - 1];
            end
            push!(Constraints, @constraint(sub_model, x[i,t] == eq))
        end
    end


    # Output Equations
    for t in 2:G.nT
        for i in 1:G.nY
            eq = 0;
            for j in 1:G.nX
                eq = eq + G.Cm[m,i,j] * x[j,t];
            end
            for j in 1:G.nU
                eq = eq + G.Dm[m,i,j] * u[j,t - 1];
            end
            push!(Constraints, @constraint(sub_model, y[i,t] == eq));
        end
    end

    # Output bounds
    for t in G.N1m[m]:G.N2m[m]
        for i in 1:G.nY
            push!(Constraints, @constraint(sub_model, G.y_min[m,i] <= y[i,t] <= G.y_max[m,i]));
        end
    end

    # bounds on control variation
    for t in 1:G.Num[m] - 1
        for i in 1:G.nU
            push!(Constraints, @constraint(sub_model, G.Du_min[m,i] <= Du[i,t] <= G.Du_max[m,i]));
        end
    end

    # bounds on control signal
    for t in 1:G.Num[m] - 1
        for i in 1:G.nU
            push!(Constraints, @constraint(sub_model, G.u_min[m,i] <= u[i,t] <= G.u_max[m,i]));
        end
    end

    # Resource constraints
    for t in 1:G.Num[m] - 1
        for r in 1:G.nR
            eq = 0;
            for j in 1:G.nU
                eq = eq + G.rmu[m,r,j] * u[j,t];
            end
            push!(Constraints, @constraint(sub_model, eq <= rsm[r,t]));
        end
    end

    # Sensitivity constraints
    for t = 1:G.Num[m] - 1
        for r = 1:G.nR
            push!(Constraints, @constraint(sub_model, rsm[r,t] - e[r,t,m] == 0));
        end
    end

    # Solve the problem
    tempo = @elapsed JuMP.optimize!.(sub_model);
    sm = JuMP.objective_value.(sub_model);
    S = (x = JuMP.value.(x), y = JuMP.value.(y), u = JuMP.value.(u), Du = JuMP.value.(Du))

    sConstraints = length(Constraints);
    it = sConstraints - G.nR * (G.Num[m] - 1);
    sensm =  zeros(G.nR, G.Num[m] - 1);
    let it = it
        for t in 1:G.Num[m] - 1
            for r in 1:G.nR
                it = it + 1;
                sensm[r,t] = JuMP.dual(Constraints[it]);
            end
        end
    end

    sm_output = (sm = sm, S = S, sensm = sensm, tempo = tempo)

    return sm_output
end
