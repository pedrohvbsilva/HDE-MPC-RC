function read_objective_bilevel(nM,nT)
    file = readlines("bilevel_data/M$(nM)T$(nT).txt")

    aux = []
    for (line,content) in enumerate(file)
        if occursin("iter", content)
            append!(aux, file[line+1:line+10])
        end
    end 
    it = []
    for line in eachindex(aux)
        try
            push!(it,parse(Float64, split(aux[line])[2]))
            catch
                break
        end
    end
    return it
end