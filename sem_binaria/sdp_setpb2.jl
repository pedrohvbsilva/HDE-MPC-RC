using Random

function set_pb(m,t,r)
    Random.seed!(2)
    # M=20, T={4,6,8,10,12} rmax = 0.5 5?
    # M=40, T={4,6,8,10,12} rmax = 2  8?
    # M=80, T={4,6,8,10,12} rmax = 10 15?
    # M=120, T={4,6,8,10,12} rmax = 25?

    nM = m;
    nX = 2;
    nY = 1;
    nU = 2;
    nT = t;
    nR = 1;

    N1m =  2*ones(nM);
    N2m = nT*ones(nM);
    Num = nT*ones(nM);

    N1m = trunc.(Int, N1m);
    N2m = trunc.(Int, N2m);
    Num = trunc.(Int, Num);

    Am = rand(nM,nX,nX); # system matrices
    Bm = rand(nM,nX,nU);

    Cm = rand(nM,nY,nX); # system output matrices
    Dm = zeros(nM,nY,nU);

    y_min  = -25*ones(nM,nY);  # output bounds
    y_max  =  25*ones(nM,nY);  # output bounds

    u_min  =  0*ones(nM,nU);  # control signal lower bound
    u_max  =  3*ones(nM,nU);  # control signal lower bound

    Du_min = -3*ones(nM,nU);    # bound on control delta
    Du_max =  3*ones(nM,nU);    # bound on control delta

    Du_min = trunc.(Int, Du_min);
    Du_max = trunc.(Int, Du_max);

    wm = 0.2*ones(nM,nY,nT);  # reference trajectory

    r_max = r*ones(nR,nT);
    rmu   = ones(nM,nR,nU);

    qm = ones(nM);      # state deviation cost
    rm = 0.1*ones(nM);  # control signal variation cost
    x0 = rand(nM,nX); # initial state
    u0 = rand(nM,nU);  # previous control signal

    return nM, nX, nY, nU, nT, nR, N1m, N2m, Num, Am, Bm, Cm, Dm, y_min, y_max, u_min, u_max,
           Du_min, Du_max, wm, r_max, rmu, qm, rm, x0, u0
end
