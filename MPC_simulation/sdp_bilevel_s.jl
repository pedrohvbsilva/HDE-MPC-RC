function sdp_bilevel_s(ex, G)
    e = zeros(G.nR, G.nT - 1, G.nM);
    it = 1;
    let it = it
        for r in 1:G.nR
            for t in 1:G.nT - 1
                for m in 1:G.nM
                    e[r,t,m] = ex[it];
                    it = it + 1;
                end
            end
        end
    end
    sm = zeros(G.nM, 1);
    tempo = zeros(G.nM, 1);
    flg = zeros(G.nM, 1);

    Sm =  Array{NamedTuple}(undef, G.nM, 1);
    sensm =  Array{NamedTuple}(undef, G.nM, 1);
    sm_output = Array{NamedTuple}(undef, G.nM, 1);
    
    let sm_output = sm_output, sensm = sensm, Sm = Sm
            for m in 1:G.nM
            sm_output[m] = sdp_bilevel_sm(e, m, G)
            sm[m] = sm_output[m].sm;
            tempo[m] = sm_output[m].tempo;
            Sm[m] = (x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du);
            sensm[m] = (grad = sm_output[m].sensm,);
        end
    end

    s = sum(sm);
    tempoo = sum(tempo);
    grad_ = zeros(G.nR, G.nT - 1, G.nM);
    for m in 1:G.nM
        it = 0;
        for t in 1:G.nT - 1
            for r in 1:G.nR
                it = it + 1
                grad_[r,t,m] = sensm[m].grad[r,t];
            end
        end
    end

    grad = [];
    for r in 1:G.nR
        for t in 1:G.nT - 1
            for m in 1:G.nM
                push!(grad, grad_[r,t,m])
            end
        end
    end

    return s, grad, sm, Sm, tempoo
end
