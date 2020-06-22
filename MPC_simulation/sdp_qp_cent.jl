function sdp_qp_cent(G)
        # model = Model(with_optimizer(Ipopt.Optimizer, tol = 1e-7))
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0))

    y = @variable(model, y[1:G.nM, 1:G.nY, 1:G.nT]);
    x = @variable(model, x[1:G.nM, 1:G.nX, 1:G.nT]);
    u = @variable(model, u[1:G.nM, 1:G.nU, 1:G.nT - 1]);
    Du = @variable(model, Du[1:G.nM, 1:G.nU, 1:G.nT - 1]);
    δ = @variable(model, δ[1:G.nM, 1:G.nT - 1], Bin);

    # QP OPtimization
    @objective(model, Min, sum(G.rm[m] * Du[m,i,t]^2 for m in 1:G.nM for i in 1:G.nU for t in 1:G.Num[m] - 1) +
                            sum(G.qm[m] * (y[m,i,t] - G.wm[m,i,t])^2 for m in 1:G.nM for t in G.N1m[m]:G.N2m[m] for i in 1:G.nY));

    # Constraints
    # Initial states
    for m in 1:G.nM
        for i in 1:G.nX
            @constraint(model, x[m,i,1] == G.x0[m,i]);
        end
    end

    # for m in 1:G.nM
    #     for t in 2:G.nT
    #         for i in 1:G.nX
    #             @constraint(model, x[m,i,t] == 0.99)
    #         end
    #     end
    # end

    # Control Variation
    for m in 1:G.nM
        for i in 1:G.nU
            @constraint(model, Du[m,i,1] == u[m,i,1] - G.u0[m,i])
        end
    end

    for m in 1:G.nM
        for t in 2:G.nT-1
            for i in 1:G.nU
                @constraint(model, Du[m,i,t] == u[m,i,t] - u[m,i,t-1])
            end
        end
    end

    #State Dynamic Equations
    for m in 1:G.nM
        for t in 2:G.nT
            for i in 1:G.nX
                eq = 0;
                for j in 1:G.nX
                    eq = eq + G.Am[m,i,j]*x[m,j,t-1];
                end
                for j in 1:G.nU
                    eq = eq + G.Bm[m,i,j]*u[m,j,t-1];
                end
                @constraint(model, x[m,i,t] == eq);
            end
        end
    end

    # Output Equations

    for m in 1:G.nM
        for t in 2:G.nT
            for i in 1:G.nY
                eq = 0;
                for j in 1:G.nX
                    eq = eq + G.Cm[m,i,j]*x[m,j,t];
                end
                for j in 1:G.nU
                    eq = eq + G.Dm[m,i,j]*u[m,j,t-1];
                end
                @constraint(model, y[m,i,t] == eq);
            end
        end
    end

    # Output bounds

    for m in 1:G.nM
        for t in G.N1m[m]:G.N2m[m]
            for i in 1:G.nY
                @constraint(model, G.y_min[m,i] <= y[m,i,t] <= G.y_max[m,i]);
            end
        end
    end

    # bounds on control variation, and control signal
    for m in 1:G.nM
        for t in 1:G.Num[m]-1
            for i in 1:G.nU
                @constraint(model, G.Du_min[m,i] <= Du[m,i,t] <= G.Du_max[m,i]);
            end
        end
    end

    for m in 1:G.nM
        for t in 1:G.Num[m]-1
            for i in 1:G.nU
                @constraint(model, G.u_min[m,i] <= u[m,i,t] <= G.u_max[m,i]);
            end
        end
    end

    # for m in 1:G.nM
    #     for t in 1:G.Num[m]-1
    #         for i in 1:G.nU
    #             @constraint(model, G.u_min[m,i] * δ[m,t] <= u[m,i,t]);
    #             @constraint(model, u[m,i,t] <= δ[m,t] * G.u_max[m,i]);
    #         end
    #     end
    # end

    # Resource constraints
    for t in 1:G.nT-1
        for r in 1:G.nR
            eq = 0;
            for m in 1:G.nM
                for j in 1:G.nU
                    eq = eq + G.rmu[m,r,j]*u[m,j,t];
                end
            end
            @constraint(model, eq <= G.r_max[r,t]);
        end
    end

    JuMP.optimize!(model)

    sol = JuMP.objective_value.(model);

    S = (x = JuMP.value.(x), y = JuMP.value.(y), u = JuMP.value.(u), Du = JuMP.value.(Du))

    return sol, S
end