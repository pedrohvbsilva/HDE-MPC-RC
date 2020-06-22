#Bilevel Optimization: 2nd level problem solver for unit "m"
function sdp_bilevel_sm(e,m)
    # sub_model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, tol = 1e-6, max_iter = 200))
    sub_model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))
    y = @variable(sub_model, y[1:nY[1], 1:nT[1]]);
    x = @variable(sub_model, x[1:nX[1], 1:nT[1]]);
    u = @variable(sub_model, u[1:nU[1], 1:nT[1]-1]);
    Du = @variable(sub_model, Du[1:nU[1], 1:nT[1]-1]);

    rsm = @variable(sub_model, rsm[1:nR[1], 1:Num[m]-1] >= 0);

    # QP OPtimization
    @objective(sub_model, Min, sum(rm1[m] * Du[i,t]^2 for i in 1:nU[1] for t in 1:Num[m] - 1) +
                           sum(qm[m] * (y[i,t]-wm[m,i,t])^2 for t in N1m[m]:N2m[m] for i in 1:nY[1]));

    # Constraints
    Constraints = [];
    # Initial states
    for i in 1:nX[1]
        push!(Constraints, @constraint(sub_model, x[i,1] == x0[m,i]));
    end

    # Control variation
    for i in 1:nU[1]
         push!(Constraints, @constraint(sub_model, Du[i,1] == u[i,1] - u0[m,i]));
    end

    for t in 2:nT[1]-1
        for i=1:nU[1]
            push!(Constraints, @constraint(sub_model, Du[i,t] == u[i,t] - u[i,t-1]));
        end
    end

    # State Dynamic Equations
    for t in 2:nT[1]
        for i in 1:nX[1]
            eq = 0;
            for j in 1:nX[1]
                eq = eq + Am[m,i,j]*x[j,t-1];
            end
            for j in 1:nU[1]
                eq = eq + Bm[m,i,j]*u[j,t-1];
            end
            push!(Constraints, @constraint(sub_model, x[i,t] == eq));
        end
    end

    # Output Equations
    for t in 2:nT[1]
        for i in 1:nY[1]
            eq = 0;
            for j in 1:nX[1]
                eq = eq + Cm[m,i,j]*x[j,t];
            end
            for j in 1:nU[1]
                eq = eq + Dm[m,i,j]*Du[j,t-1];
            end
            push!(Constraints, @constraint(sub_model, y[i,t] == eq));
        end
    end

    # Output bounds
    for t in N1m[m]:N2m[m]
        for i in 1:nY[1]
            push!(Constraints, @constraint(sub_model, y_min[m,i] <= y[i,t] <= y_max[m,i]));
        end
    end

    # bounds on control variation, and control signal
    for t in 1:Num[m]-1
        for i in 1:nU[1]
            push!(Constraints, @constraint(sub_model, Du_min[m,i] <= Du[i,t] <= Du_max[m,i]));
        end
    end

    for t in 1:Num[m]-1
        for i in 1:nU[1]
            push!(Constraints, @constraint(sub_model, u_min[m,i] <= u[i,t] <= u_max[m,i]));
        end
    end

    # Resource constraints
    for t in 1:Num[m]-1
        for r in 1:nR[1]
            eq = 0;
            for j in 1:nU[1]
                eq = eq + rmu[m,r,j]*u[j,t];
            end
            push!(Constraints, @constraint(sub_model, eq <= rsm[r,t]));
        end
    end

    # Sensitivity constraints
    for t=1:Num[m]-1
        for r=1:nR[1]
            push!(Constraints, @constraint(sub_model, rsm[r,t] - e[r,t,m] == 0));
        end
    end

    # Solve the problem
    tempo = @elapsed JuMP.optimize!.(sub_model);
    sm = JuMP.objective_value.(sub_model);
    S = (x = JuMP.value.(x), y = JuMP.value.(y), u = JuMP.value.(u), Du = JuMP.value.(Du))

    sConstraints = length(Constraints);
    it = sConstraints - nR[1]*(Num[m]-1);
    sensm =  zeros(nR[1], Num[m]-1);
    let it = it
        for t in 1:Num[m]-1
             for r in 1:nR[1]
                 it = it+1;
                 sensm[r,t] = JuMP.dual(Constraints[it]);
            end
        end
    end

    sm_output = (sm = sm,S = S,sensm = sensm, tempo = tempo)

return sm_output
end
