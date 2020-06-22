include("sdp_setpb2.jl")

# M = [4, 8, 12, 20];
# T = [4, 6, 8, 10];
# M = [4, 8, 10, 20];
T = [4];
M = [4];

for m in M
    for t in T
        r = 2
        global nM, nX, nY, nU, nT, nR, N1m, N2m, Num, Am, Bm, Cm, Dm, y_min, y_max, u_min, u_max, Du_min, Du_max, wm, r_max, rmu, qm, rm, x0, u0 = set_pb(m, t, r)
        println("nM = ", nM, " nT = ", nT)
        # include("sdp_qp_cent_bin.jl")
        include("sdp_benders_master_bin.jl")
    end
end 
