using JuMP, Gurobi, LinearAlgebra
# include("sdp_setprob_sbai.jl")
# include("sdp_setpb_dissertation.jl")
# include("sdp_setpbcamps.jl")
# include("sdp_setpb.jl")
const GRB_ENV = Gurobi.Env()

model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), "Presolve" => 0))
# model = Model(Ipopt.Optimizer)

@variable(model, y[1:nM[1],  1:nY[1], 1:nT[1]]);
@variable(model, x[1:nM[1], 1:nX[1], 1:nT[1]]);
@variable(model, u[1:nM[1], 1:nU[1], 1:nT[1] - 1]);
@variable(model, Du[1:nM[1], 1:nU[1], 1:nT[1] - 1]);

# -----------------------
# QP OPtimization
@objective(model, Min, sum(rm1[m] * Du[m,i,t]^2 for m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1) + 
                                sum(qm[m] * (y[m,i,t] - wm[m,i,t])^2 for m in 1:nM[1], t in N1m[m]:N2m[m], i in 1:nY[1]));

# --------------------------
# Constraints

# Initial states
@constraint(model, x0_cons[m in 1:nM[1], i in 1:nX[1]], x[m,i,1] == x0[m,i]);

# Control variation
@constraint(model, Δu_cons[m in 1:nM[1], i in 1:nU[1], 1], Du[m,i,1] == u[m,i,1] - u0[m,i]);

@constraint(model, Δu_cons2[m in 1:nM[1], i in 1:nU[1], t in 2:nT[1] - 1], Du[m,i,t] == u[m,i,t] - u[m,i,t - 1]);

# # State Dynamic Equations
# @constraint(model, x_eq_cons[i in 1:nX[1], t in 2:nT[1]], x[i,t] == sum(Am[m,i,j] * x[j,t - 1] for j in 1:nX[1]) + sum(Bm[m,i,j] * u[j,t - 1] for j in 1:nU[1]));

# # Output Equations
# @constraint(model, y_eq_cons[i in 1:nY[1], t in 2:nT[1]], y[i,t] == sum(Cm[m,i,j] * x[j,t] for j in 1:nX[1]) + sum(Dm[m,i,j] * u[j,t-1] for j in 1:nU[1]));

# State Dynamic Equations
@constraint(model, x_eq_cons[m in 1:nM[1], i in 1:nX[1], t in 2:nT[1]], x[m,i,t] == dot(Am[m,i,:], x[m,:,t - 1]) + dot(Bm[m,i,:], Du[m,:,t - 1]));

# Output Equations
@constraint(model, y_eq_cons[m in 1:nM[1], i in 1:nY[1], t in 2:nT[1]], y[m,i,t] == dot(Cm[m,i,:], x[m,:,t]) + dot(Dm[m,i,:], Du[m,:,t - 1]));

# Output bounds 
@constraint(model, y_min_cons[m in 1:nM[1], i in 1:nY[1], t in N1m[m]:N2m[m]], y_min[m,i] - y[m,i,t] <= 0.);
@constraint(model, y_max_cons[m in 1:nM[1], i in 1:nY[1], t in N1m[m]:N2m[m]], y[m,i,t] - y_max[m,i] <= 0.);

# bounds on control variation
@constraint(model, Δu_min_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], Du_min[m,i] - Du[m,i,t] <= 0.);
@constraint(model, Δu_max_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], Du[m,i,t] - Du_max[m,i] <= 0.);

# bounds on control signal
@constraint(model, u_min_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], u_min[m,i] - u[m,i,t] <= 0.);
@constraint(model, u_max_cons[m in 1:nM[1], i in 1:nU[1], t in 1:Num[m] - 1], u[m,i,t] - u_max[m,i] <= 0.);

# Resource constraints
@constraint(model, res_cons[r in 1:nR[1], t in 1:nT[1] - 1], sum(rmu[m,r,j] * u[m,j,t] for j in 1:nU[1], m in 1:nM[1]) <= r_max[r,t]);

# Solve the problem
tempo = @elapsed JuMP.optimize!(model)

if has_values(model)
    sol = JuMP.objective_value(model)
    ts = JuMP.termination_status(model)
    println("Objective value: ", sol)
    open("cent.csv", "a") do io
        println(io, "$tempo, $sol, $ts")
    end;
else   
    open("cent.csv", "a") do io
        println(io, "Modelo Infactivel")
    end;
end