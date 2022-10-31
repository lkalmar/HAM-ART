using DelimitedFiles
using LinearAlgebra
using Printf

data = readdlm(ARGS[1])
global threshold = 0.12

sorted_ids = sort(unique(data[:, 1]))
n = length(sorted_ids)

global multi = ones(n)
global cladeHeight = zeros(n)
global is = indexin(data[:,1], sorted_ids)
global js = indexin(data[:,2], sorted_ids)

D = zeros(n, n)
for k = 1:n*n
    D[is[k], js[k]] = data[k, 3]
end
D += I

global ids = map(x -> [x], sorted_ids)
global clades = deepcopy(ids)
global DD = copy(D)
global pass = 0

while size(DD, 1) > 2
    global DD
    global pass
    global is
    global js
    global cladeHeight
    global clades
    all(isone, DD) && return
    i, j = argmin(DD).I
    aktmin = DD[i,j]
    if aktmin > threshold && pass == 0
        pass = 1
        outf = open("clade_text_annotation.txt", "w")
        println(outf, "DATASET_TEXT")
        println(outf, "SEPARATOR TAB")
        println(outf, "DATASET_LABEL\tClade labels")
        println(outf, "DATA")
        for k = 1:size(clades,1)
            for l = 1:size(clades[k],1)
                println(outf, clades[k][l], "\t", "clade_", @sprintf("%05d", k), "\t", "-1\t#000000\tnormal\t1\t0")
            end 
        end
        close(outf)
        println("Number of clades: ", string(size(ids,1)))
    end
    if multi[i] == 1
        h1 = aktmin/2
    else
        h1 = aktmin/2 - cladeHeight[i]
    end
    if multi[j] == 1
        h2 = aktmin/2
    else
        h2 = aktmin/2 - cladeHeight[j]
    end
    cladeHeight[i] += h1
    newStr = "(" * ids[i][1] * ":" * string(h1) * "," * ids[j][1] * ":" * string(h2) * ")"
    ave = (view(DD, :, i) .* multi[i] + view(DD, :, j) .* multi[j]) ./ (multi[i] + multi[j])
    DD[:, i] = ave
    DD[i, :] = ave
    DD[i, i] = 1.0
    idx = collect(1:size(DD, 1))
    deleteat!(idx, j)
    append!(clades[i], clades[j])
    
    DD = view(DD, idx, idx)
    ids[i][1] = newStr
    multi[i]+=1
    #append!(ids[i], ids[j])
    deleteat!(multi, j)
    deleteat!(cladeHeight, j)
    deleteat!(ids, j)
    deleteat!(clades, j)
end

s = open("pairs_upgma.tre", "w")
println(s, "(" * ids[1][1] * ":" * string(DD[1,2]/2 - cladeHeight[1]) * "," * ids[2][1] * ":" * string(DD[2,1]/2 - cladeHeight[2]) * ")")
close(s)
