function sdp_bilevel_s(ex)
   
    e = zeros(nR[1],nT[1]-1,nM[1]);
    it = 1;
    let it = it
        for r in 1:nR[1]
            for t in 1:nT[1]-1
                for m in 1:nM[1]
                    e[r,t,m] = ex[it];
                    it = it+1;
                end
            end
        end
    end
   
    sm = zeros(nM[1],1);
    flg = zeros(nM[1],1);
    tempo = zeros(nM[1],1);

    Sm =  Array{NamedTuple}(undef, nM[1],nM[1]);
    sensm =  Array{NamedTuple}(undef, nM[1],nM[1]);
    sm_output = Array{NamedTuple}(undef, nM[1],nM[1]);

    teste = []
    for i in 1:nM[1]
        push!(teste, e)
    end
 
    sm_output = pmap(sdp_bilevel_sm,teste, 1:nM[1])

    let sm_output = sm_output
        for m in 1:nM[1]
            sm[m] = sm_output[m].sm;
            tempo[m] = sm_output[m].tempo
            Sm[m] = (x = sm_output[m].S.x, y = sm_output[m].S.y, u = sm_output[m].S.u, Du = sm_output[m].S.Du);
            sensm[m] = (grad = sm_output[m].sensm,);
        end
    end
    s = sum(sm);
    if nprocs() > 1
        tempo = sum(tempo)/nworkers()
    else
        tempo = sum(tempo)
    end

    grad_ = zeros(nR[1],nT[1]-1,nM[1]);
    for m in 1:nM[1]
        it = 0;
        for t in 1:nT[1]-1
            for r in 1:nR[1]
                it = it+1
                grad_[r,t,m] = sensm[m].grad[r,t];
            end
        end
    end

    grad = [];
    for r in 1:nR[1]
        for t in 1:nT[1]-1
            for m in 1:nM[1]
               push!(grad, grad_[r,t,m])
            end
        end
    end

    return s, grad, sm, Sm, tempo
end
