#############
# Solve Benders Subproblem

function sdp_OA_sub(ss, δ)

    sm = zeros(nM[1], 1);
    Sm = Array{NamedTuple}(undef, nM[1]);
    flg = zeros(nM[1], 1);
    flgc1 = [];
    flgc2 = [];

    sm_output = Array{NamedTuple}(undef, nM[1]);

    println("Trying to Solve OPTIMALITY Problems!!!");

    ar = []
    ard = []
    for m in 1:nM[1]
        push!(ar, ss)
        push!(ard, δ)
    end

    tempo1 = @elapsed sm_output = pmap(sdp_OA_submf, ar, ard, 1:nM[1])

    for m in 1:nM[1]
        # sm_output[m] = sdp_OA_submf(ss, δ, m);
        Sm[m] = (sm = sm_output[m].sm, x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du);
        sm[m] = sm_output[m].sm;
              
        flg[m] = sm_output[m].flg;

        push!(flgc1, sm_output[m].flgc1);
        push!(flgc2, sm_output[m].flgc2);
    end
 
    flg_all = maximum(flg);

    # All subproblems are feasible

    if flg_all == 0 # All subproblems were feasible
    
        smh = sum(sm)
        
        cut_OA = Array{NamedTuple}(undef, 1);
        cut_OA[1] = (smh = smh, Smh = Sm)

        return (sm, Sm, flg_all, flgc1, flgc2, cut_OA, tempo1);
    end

    # There are infeasible subproblems

    println("-----------------------------");
    println("sdp_OA_sub: Entering Feasibility Phase -- 1");
    sm = zeros(nM[1]);
    Sm = Array{NamedTuple}(undef, nM[1]);
    
    # v = findall(x->x == 1, vec(flg))
    
    flg = zeros(nM[1]);
    flgc1 = [];
    flgc2 = [];

    sm_output = Array{NamedTuple}(undef, nM[1]);
    
    tempo1 = @elapsed sm_output = pmap(sdp_OA_submi, ar, ard, 1:nM[1])

    for m in 1:nM[1]
    # for m in v
        # sm_output[m] = sdp_OA_submi(ss, δ, m);
        Sm[m] = (sm = sm_output[m].sm, x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du);
        sm[m] = sm_output[m].sm;
        
        flg[m] = sm_output[m].flg;
        push!(flgc1, sm_output[m].flgc1);
        push!(flgc2, sm_output[m].flgc2);
    end

    smh = 0.0
    cut_OA = Array{NamedTuple}(undef, 1);
    cut_OA[1] = (smh = smh, Smh = Sm)

    return (sm, Sm, flg_all, flgc1, flgc2, cut_OA, tempo1);
end