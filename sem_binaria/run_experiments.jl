using Distributed
@everywhere include("sdp_setpb2.jl")
# include("sdp_setpbcamps.jl")

# M = [20, 40, 80];      
# T = [4, 6, 8, 10];

M = [4];
T= [10];


for m in M
    @eval @everywhere m = $m
    for t in T
        @eval @everywhere t = $t
        r = 2
        @eval @everywhere r = $r
        @everywhere nM, nX, nY, nU, nT, nR, N1m, N2m, Num, Am, Bm, Cm, Dm, y_min, y_max, u_min, u_max, Du_min, Du_max, wm, r_max, rmu, qm, rm1, x0, u0 = set_pb(m,t,r)
        
        println("nM = ", nM, " nT = ", nT)
        include("sdp_qp_cent.jl")
        # include("sdp_benders_master.jl")
        # include("sdp_bilevel_master_pmap.jl")
    end
end
