#############
# Solve Benders Subproblem

function sdp_benders_sub(ss,G)

    sm = zeros(G.nM);
    Sm = Array{NamedTuple}(undef, G.nM);
    flg = zeros(G.nM);

    sm_output = Array{NamedTuple}(undef, G.nM);

    # println("Trying to Generate OPTIMALITY Cuts!!!");

    ar = []
    ar2 = []
    for m in 1:G.nM
        push!(ar, ss)
        push!(ar2, G)
    end

    # sm_output = pmap(sdp_benders_submf, ar, ar2, 1:G.nM)

    for m in 1:G.nM
        # sm_output[m] = sdp_benders_submf_pmap(ss, δ, m);
        sm_output[m] = sdp_benders_submf(ss, G, m)
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
        cut_coef = zeros(G.nR, G.nT - 1, G.nM);
 
        # println("Generating OPTIMALITY Cuts!!!");

        for m in 1:G.nM
            # Output bounds
            v = sum(Sm[m].dual_nu_y1[i,t] * (G.y_min[m,i] - Sm[m].y[i,t]) for t in G.N1m[m]:G.N2m[m], i in 1:G.nY);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_y2[i,t] * (Sm[m].y[i,t] - G.y_max[m,i]) for t in G.N1m[m]:G.N2m[m], i in 1:G.nY);
            # v = max(0,v)
            cut_cte += v

            # bounds on control variation
            v = sum(Sm[m].dual_nu_du1[i,t] * (G.Du_min[m,i] - Sm[m].Du[i,t]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_du2[i,t] * (Sm[m].Du[i,t] - G.Du_max[m,i]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v

            # bounds on control
            v = sum(Sm[m].dual_nu_u1[i,t] * (G.u_min[m,i] - Sm[m].u[i,t]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_u2[i,t] * (Sm[m].u[i,t] - G.u_max[m,i]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v
            
            v = sum(dot(G.rmu[m,r,:], Sm[m].u[:,t]) * Sm[m].dual_mu[r,t] for t in 1:G.nT - 1, r in 1:G.nR);
            # v = max(0,v) # This is valid if control signal u >= 0
            cut_cte += v;                    
            cut_coef[:,:,m] = - Sm[m].dual_mu;

        end

        cut_benders = Array{NamedTuple}(undef, 1);
        cut_benders[1] = (cut_cte = cut_cte, cut_coef = cut_coef);

        return (sm, Sm, flg_all, flg, cut_benders);
    end

    # There are infeasible subproblems

    # println("-----------------------------");
    # println("sdp_benders_sub: Entering Feasibility Phase -- 1");
    sm = zeros(G.nM);
    Sm = Array{NamedTuple}(undef, G.nM);

    # k = findall(x->x == 1, vec(flg))
    
    flg = zeros(G.nM);

    sm_output = Array{NamedTuple}(undef, G.nM);

    # sm_output = pmap(sdp_benders_submi, ar, ar2, 1:G.nM)

    for m in 1:G.nM
        # sm_output[m] = sdp_benders_submi_pmap(ss,m);
        sm_output[m] = sdp_benders_submi(ss, G, m)
        sm[m] = sm_output[m].sm;
        Sm[m] = (sm = sm_output[m].sm, x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du,
                dual_nu_y1 = sm_output[m].dual_nu_y1, dual_nu_y2 = sm_output[m].dual_nu_y2,
                dual_nu_u1 = sm_output[m].dual_nu_u1, dual_nu_u2 = sm_output[m].dual_nu_u2,
                dual_nu_du1 = sm_output[m].dual_nu_du1, dual_nu_du2 = sm_output[m].dual_nu_du2,
                dual_mu = sm_output[m].dual_mu);

        flg[m] = sm_output[m].flg;
    end

    gamma_value = max(maximum(sm), 1e-6);

    # Compute cut
    cut_cte = 0.0;
    cut_coef = zeros(G.nR, G.nT - 1, G.nM);

    for m in 1:G.nM
        if Sm[m].sm > 0
            # if Sm[m].sm <= 0
            # if Sm[m].sm <= gamma_value
            # if Sm[m].sm >= gamma_value
           
            # Output bounds
            v = sum(Sm[m].dual_nu_y1[i,t] * (G.y_min[m,i] - Sm[m].y[i,t]) for t in G.N1m[m]:G.N2m[m], i in 1:G.nY);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_y2[i,t] * (Sm[m].y[i,t] - G.y_max[m,i]) for t in G.N1m[m]:G.N2m[m], i in 1:G.nY);
            # v = max(0,v)
            cut_cte += v
    
            # bounds on control variation   
            v = sum(Sm[m].dual_nu_du1[i,t] * (G.Du_min[m,i] - Sm[m].Du[i,t]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_du2[i,t] * (Sm[m].Du[i,t] - G.Du_max[m,i]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v
            
            # bounds on control
            v = sum(Sm[m].dual_nu_u1[i,t] * (G.u_min[m,i] - Sm[m].u[i,t]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v
            v = sum(Sm[m].dual_nu_u2[i,t] * (Sm[m].u[i,t] - G.u_max[m,i]) for t in 1:G.Num[m] - 1, i in 1:G.nU);
            # v = max(0,v)
            cut_cte += v
    
            v = sum(dot(G.rmu[m,r,:], Sm[m].u[:,t]) * Sm[m].dual_mu[r,t] for t in 1:G.nT - 1, r in 1:G.nR);
            # v = max(0,v) # This is valid if control signal u >= 0
            cut_cte += v;                    
            cut_coef[:,:,m] = - Sm[m].dual_mu;
        end
    end
    # println("CTEi: ", cut_cte)
    # println("COEFi: ", cut_coef)
    # println("COEFi: ", cut_coef)
    # println("COEF2i: ", cut_coef2)
    cut_benders = Array{NamedTuple}(undef, 1);
    cut_benders[1] = (cut_cte = cut_cte, cut_coef = cut_coef);

    return (sm, Sm, flg_all, flg, cut_benders);
end