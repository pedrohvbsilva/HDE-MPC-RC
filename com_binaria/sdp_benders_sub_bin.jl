#############
# Solve Benders Subproblem

function sdp_benders_sub_bin(ss, δ)

    sm = zeros(nM[1]);
    Sm = Array{NamedTuple}(undef, nM[1]);
    flg = zeros(nM[1]);

    sm_output = Array{NamedTuple}(undef, nM[1]);

    println("Trying to Generate OPTIMALITY Cuts!!!");

    ar = []
    ard = []
    for m in 1:nM[1]
        push!(ar, ss)
        push!(ard, δ)
    end

    tempo1 = @elapsed sm_output = pmap(sdp_benders_submf_bin, ar, ard, 1:nM[1])

    for m in 1:nM[1]
        # sm_output[m] = sdp_benders_submf_bin(ss, δ, m);
        sm[m] = sm_output[m].sm;
        Sm[m] = (x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du,
                    dual_nu_y1 = sm_output[m].dual_nu_y1, dual_nu_y2 = sm_output[m].dual_nu_y2,
                    dual_nu_u1 = sm_output[m].dual_nu_u1, dual_nu_u2 = sm_output[m].dual_nu_u2,
                    dual_nu_du1 = sm_output[m].dual_nu_du1, dual_nu_du2 = sm_output[m].dual_nu_du2,
                    dual_mu = sm_output[m].dual_mu);
        # Sm[m] = (x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du,
        #             dual_nu_y1 = max.(sm_output[m].dual_nu_y1, 0), dual_nu_y2 = max.(sm_output[m].dual_nu_y2, 0),
        #             dual_nu_u1 = max.(sm_output[m].dual_nu_u1, 0), dual_nu_u2 = max.(sm_output[m].dual_nu_u2, 0),
        #             dual_nu_du1 = max.(sm_output[m].dual_nu_du1, 0), dual_nu_du2 = max.(sm_output[m].dual_nu_du2, 0),
        #             dual_mu = max.(sm_output[m].dual_mu, 0));

        flg[m] = sm_output[m].flg;
    end

    flg_all = maximum(flg);

    # All subproblems are feasible

    if flg_all == 0 # All subproblems were feasible

        # Compute cut
        cut_cte = sum(sm);
        cut_coef = zeros(nR[1], nT[1] - 1, nM[1]);
        cut_cte2 = 0.0;
        cut_coef2 = zeros(nT[1] - 1, nM[1]);

        println("Generating OPTIMALITY Cuts!!!");

        for m in 1:nM[1]
            # Output bounds
            v = sum(Sm[m].dual_nu_y1[i,t] * (y_min[m,i] - Sm[m].y[i,t]) for t in N1m[m]:N2m[m], i in 1:nY[1]);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_y2[i,t] * (Sm[m].y[i,t] - y_max[m,i]) for t in N1m[m]:N2m[m], i in 1:nY[1]);
            # v = max(0,v)
            cut_cte += v

            # bounds on control variation
            v = sum(Sm[m].dual_nu_du1[i,t] * (Du_min[m,i] - Sm[m].Du[i,t]) for t in 1:Num[m] - 1, i in 1:nU[1]);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_du2[i,t] * (Sm[m].Du[i,t] - Du_max[m,i]) for t in 1:Num[m] - 1, i in 1:nU[1]);
            # v = max(0,v)
            cut_cte += v
            
            v = sum(dot(rmu[m,r,:], Sm[m].u[:,t]) * Sm[m].dual_mu[r,t] for t in 1:nT[1] - 1, r in 1:nR[1]);
            # v = max(0,v) # This is valid if control signal u >= 0
            cut_cte += v;                    
            cut_coef[:,:,m] = - Sm[m].dual_mu;

            v = sum(dot(Sm[m].dual_nu_u2[:,t], Sm[m].u[:,t]) - dot(Sm[m].dual_nu_u1[:,t], Sm[m].u[:,t]) for t in 1:Num[m] - 1, r in 1:nR[1]);
            # v = sum((Sm[m].dual_nu_u2[j,t] * Sm[m].u[j,t]) - (Sm[m].dual_nu_u1[j,t] * Sm[m].u[j,t]) for t in 1:nT[1] - 1, r in 1:nR[1], j in 1:nU[1]);
            # v = max(0,v)              
            cut_cte2 += v;    

            # for t in 1:Num[m] - 1           
            #     for j in 1:nU[1]
            #         cut_coef2[t,m] = ((Sm[m].dual_nu_u1[j,t] * u_min[m,j]) - (Sm[m].dual_nu_u2[j,t] * u_max[m,j]));
            #     end
            # end
            for t in 1:Num[m] - 1           
                cut_coef2[t,m] = (dot(Sm[m].dual_nu_u1[:,t], u_min[m,:]) - dot(Sm[m].dual_nu_u2[:,t], u_max[m,:]));
            end  
        end
        # println("CTEf: ", cut_cte)
        # println("CTE2f: ", cut_cte2)
        # println("COEFf: ", cut_coef)
        # println("COEF2f: ", cut_coef2)
        
        cut_benders = Array{NamedTuple}(undef, 1);
        cut_benders[1] = (cut_cte = cut_cte, cut_coef = cut_coef, cut_cte2 = cut_cte2, cut_coef2 = cut_coef2);

        return (sm, Sm, flg_all, flg, cut_benders, tempo1);
    end

    # There are infeasible subproblems

    println("-----------------------------");
    println("sdp_benders_sub: Entering Feasibility Phase -- 1");
    sm = zeros(nM[1]);
    Sm = Array{NamedTuple}(undef, nM[1]);

    # k = findall(x->x == 1, vec(flg))

    flg = zeros(nM[1]);

    sm_output = Array{NamedTuple}(undef, nM[1]);
    
    ar3 = []
    ar4 = []
    for _ in 1:nM[1]
        push!(ar3, 0)
        push!(ar4, 0)
    end

    tempo1 = @elapsed sm_output = pmap(sdp_benders_submi_bin, ar, ar3, ar4, ard, 1:nM[1])

    # for m in k
    for m in 1:nM[1]
        # sm_output[m] = sdp_benders_submi_bin(ss, 0, 0, δ, m);
        sm[m] = sm_output[m].sm;
        Sm[m] = (sm = sm_output[m].sm, x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du,
                dual_nu_y1 = sm_output[m].dual_nu_y1, dual_nu_y2 = sm_output[m].dual_nu_y2,
                dual_nu_u1 = sm_output[m].dual_nu_u1, dual_nu_u2 = sm_output[m].dual_nu_u2,
                dual_nu_du1 = sm_output[m].dual_nu_du1, dual_nu_du2 = sm_output[m].dual_nu_du2,
                dual_mu = sm_output[m].dual_mu);
        # Sm[m] = (sm = sm_output[m].sm, x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du,
        #         dual_mu_x0_A = sm_output[m].dual_mu_x0_A, dual_mu_x0_B = sm_output[m].dual_mu_x0_B,     
        #         dual_mu_Du0_A = sm_output[m].dual_mu_Du0_A, dual_mu_Du0_B = sm_output[m].dual_mu_Du0_B,
        #         dual_mu_Du_A = sm_output[m].dual_mu_Du_A, dual_mu_Du_B = sm_output[m].dual_mu_Du_B,   
        #         dual_mu_X_A = sm_output[m].dual_mu_X_A, dual_mu_X_B = sm_output[m].dual_mu_X_B,
        #         dual_mu_Y_A = sm_output[m].dual_mu_Y_A, dual_mu_Y_B = sm_output[m].dual_mu_Y_B,
        #         dual_nu_y1 = sm_output[m].dual_nu_y1, dual_nu_y2 = sm_output[m].dual_nu_y2,
        #         dual_nu_u1 = sm_output[m].dual_nu_u1, dual_nu_u2 = sm_output[m].dual_nu_u2,
        #         dual_nu_du1 = sm_output[m].dual_nu_du1, dual_nu_du2 = sm_output[m].dual_nu_du2,
        #         dual_mu = sm_output[m].dual_mu);

        # Sm[m].solm = sm_output[m].solm;
        flg[m] = sm_output[m].flg;
    end

    gamma_value = max(maximum(sm), 1e-6);

    # Compute cut
    cut_cte = 0.0;
    cut_coef = zeros(nR[1], nT[1] - 1, nM[1]);
    cut_cte2 = 0.0;
    cut_coef2 = zeros(nT[1] - 1, nM[1]);

    # for m in 1:nM[1]

    #     if Sm[m].sm <= gamma_value
           
    #         cut_cte += sum(Sm[m].dual_mu_x0_A[i] * (Sm[m].x[i,1] - x0[m,i]) + 
    #                        Sm[m].dual_mu_x0_B[i] * (-(Sm[m].x[i,1] - x0[m,i])) for i in 1:nX[1]);

    #         cut_cte += sum(Sm[m].dual_mu_Du0_A[i] * (Sm[m].Du[i,1] - (Sm[m].u[i,1] - u0[m,i])) + 
    #                        Sm[m].dual_mu_Du0_B[i] * (-(Sm[m].Du[i,1] - (Sm[m].u[i,1] - u0[m,i]))) for i in 1:nU[1]);

    #         cut_cte += sum(Sm[m].dual_mu_Du_A[t,i] * (   Sm[m].Du[i,t] - (Sm[m].u[i,t] - Sm[m].u[i,t - 1])  ) + 
    #                        Sm[m].dual_mu_Du_B[t,i] * ( -(Sm[m].Du[i,t] - (Sm[m].u[i,t] - Sm[m].u[i,t - 1])) ) for t in 2:nT[1] - 1, i in 1:nU[1]);
                 
    #         cut_cte += sum(Sm[m].dual_mu_X_A[t,i] * (Sm[m].x[i,t] - (dot(Am[m,i,:], Sm[m].x[:,t - 1]) + dot(Bm[m,i,:], Sm[m].u[:,t - 1]))) + 
    #                        Sm[m].dual_mu_X_B[t,i] * - (Sm[m].x[i,t] - (dot(Am[m,i,:], Sm[m].x[:,t - 1]) + dot(Bm[m,i,:], Sm[m].u[:,t - 1]))) for t in 2:nT[1], i in 1:nX[1])

    #         cut_cte += sum(Sm[m].dual_mu_Y_A[t,i] * (Sm[m].y[i,t] - (dot(Cm[m,i,:], Sm[m].x[:,t]) + dot(Dm[m,i,:], Sm[m].u[:,t - 1]))) + 
    #                        Sm[m].dual_mu_Y_B[t,i] * - (Sm[m].y[i,t] - (dot(Cm[m,i,:], Sm[m].x[:,t]) + dot(Dm[m,i,:], Sm[m].u[:,t - 1]))) for t in 2:nT[1], i in 1:nY[1])
    #     end
    # end

    # for m in k
    for m in 1:nM[1]
        if Sm[m].sm > 0
        # if Sm[m].sm <= 0
        # if Sm[m].sm <= gamma_value
        # if Sm[m].sm >= gamma_value
       
            # Output bounds
            v = sum(Sm[m].dual_nu_y1[i,t] * (y_min[m,i] - Sm[m].y[i,t]) for t in N1m[m]:N2m[m], i in 1:nY[1]);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_y2[i,t] * (Sm[m].y[i,t] - y_max[m,i]) for t in N1m[m]:N2m[m], i in 1:nY[1]);
            # v = max(0,v)
            cut_cte += v

            # bounds on control variation
            v = sum(Sm[m].dual_nu_du1[i,t] * (Du_min[m,i] - Sm[m].Du[i,t]) for t in 1:Num[m] - 1, i in 1:nU[1]);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_du2[i,t] * (Sm[m].Du[i,t] - Du_max[m,i]) for t in 1:Num[m] - 1, i in 1:nU[1]);
            # v = max(0,v)
            cut_cte += v

            v = sum(dot(rmu[m,r,:], Sm[m].u[:,t]) * Sm[m].dual_mu[r,t] for t in 1:nT[1] - 1, r in 1:nR[1]);
            # v = max(0,v) # This is valid if control signal u >= 0
            cut_cte += v;                    
            cut_coef[:,:,m] = - Sm[m].dual_mu;
 
            # bounds on control signal
            v = sum(dot(Sm[m].dual_nu_u2[:,t], Sm[m].u[:,t]) - dot(Sm[m].dual_nu_u1[:,t], Sm[m].u[:,t]) for t in 1:Num[m] - 1, r in 1:nR[1]);
            # v = sum((Sm[m].dual_nu_u2[j,t] * Sm[m].u[j,t]) - (Sm[m].dual_nu_u1[j,t] * Sm[m].u[j,t]) for t in 1:nT[1] - 1, r in 1:nR[1], j in 1:nU[1]);
            # v = max(0,v)              
            cut_cte2 += v;    

            # for t in 1:Num[m] - 1           
            #     for j in 1:nU[1]
            #         cut_coef2[t,m] = ((Sm[m].dual_nu_u1[j,t] * u_min[m,j]) - (Sm[m].dual_nu_u2[j,t] * u_max[m,j]));
            #     end
            # end  

            for t in 1:Num[m] - 1           
                cut_coef2[t,m] = (dot(Sm[m].dual_nu_u1[:,t], u_min[m,:]) - dot(Sm[m].dual_nu_u2[:,t], u_max[m,:]));
            end  
            
        end
    end
    # println("CTEi: ", cut_cte)
    # println("CTE2i: ", cut_cte2)
    # println("COEFi: ", cut_coef)
    # println("COEF2i: ", cut_coef2)
    cut_benders = Array{NamedTuple}(undef, 1);
    cut_benders[1] = (cut_cte = cut_cte, cut_coef = cut_coef, cut_cte2 = cut_cte2, cut_coef2 = cut_coef2);
    
    return (sm, Sm, flg_all, flg, cut_benders, tempo1);
end