function create_m_models(A, B, C, D, nM, nX, nY, nU, nR, nT, x0, u0)

    N1m =  2*ones(nM);  N1m = trunc.(Int, N1m);
    N2m = nT*ones(nM);  N2m = trunc.(Int, N2m);
    Num = nT*ones(nM);  Num = trunc.(Int, Num);

    Am = Array{Float64}(undef, nM, nX, nX);
    Bm = Array{Float64}(undef, nM, nX, nU);
    Cm = Array{Float64}(undef, nM, nY, nX);
    Dm = Array{Float64}(undef, nM, nY, nU);

    for m in 1:nM
        Am[m, 1, 1] = A;
        Bm[m, 1, 1] = B;
        Cm[m, 1, 1] = C
        Dm[m, 1, 1] = D;
    end

    y_min  =    0*ones(nM,nY);   # output bounds
    y_max  =  2*ones(nM,nY);     # output bounds

    u_min  =    20*ones(nM,nU);       # control signal lower bound
    u_max  =  50*ones(nM,nU);    # control signal lower bound

    Du_min = -50*ones(nM,nU);    # bound on control delta
    Du_max =  50*ones(nM,nU);    # bound on control delta

    Du_min = trunc.(Int, Du_min);
    Du_max = trunc.(Int, Du_max);

    wm = 1*ones(nM,nY,nT);    # reference trajectory

    r_max = 100*ones(nR,nT);
    rmu   = ones(nM,nR,nU);

    qm = ones(nM);              # state deviation cost
    rm = 0.000001*ones(nM);          # control signal variation cost

    G = (Am = Am, Bm = Bm, Cm = Cm, Dm = Dm, x0 = x0, u0 = u0,
        u_min = u_min, u_max = u_max, y_min = y_min, y_max = y_max, 
        Du_min = Du_min, Du_max = Du_max, wm = wm, r_max = r_max, 
        rmu = rmu, qm = qm, rm = rm, N1m = N1m, N2m = N2m, Num = Num, 
        nM = nM, nX = nX, nY = nY, nU = nU, nR = nR, nT = nT);

    return G
end