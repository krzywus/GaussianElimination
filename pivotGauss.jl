include("fileOperations.jl")
using fileOperations

function pivotalGaussianElimination(A, b, n, l)
# - Gaussian Elimination -
    p = collect(1:n)
    s = zeros(n)
    count = 0
    for i in 1:n # rows iterator, this block computes the array of row maximal elements
        for j in getNonZeroElementsIndexes(i, n, l)
            s[i] = max(s[i], abs(A[i, j] ))
        end
    end
    for k in 1:(n-1)
        rmax = 0 # this block finds the largest scaled column entry
        j = -1
        for i in k:n
            r = abs(A[p[i], k] / s[p[i]])
            if r > rmax
                rmax = r
                j = i # row index of largest scaled entry
            end
        end
        p[k], p[j] = p[j], p[k] # exchange row pointers
        for i in (k+1):getLimit(k%l, k, l, n) # perform elimination on submatrix
            pi, pk = p[i], p[k]
            A[pi, k]  = A[pi, k] / A[pk, k]
            count += 1
            for j in (k+1):n
                if A[pk, j] != 0
                    A[pi, j] -= A[pi, k] * A[pk, j]
                    count += 1
                end
            end
        end
        # println("$k")
    end
    println("count:$count")
    count += forwardPivotal(A, b, p, n, l)
    println("count:$count")
    x, tmpCount = backwardPivotal(A, b, p, n, l)
    count += tmpCount
    println("count:$count")
    return x
end

function forwardPivotal(A, b, p, n, l)
# - Forward Elimination -
    count = 0
    for k in 1:(n-1)
        for i in (k+1):getLimit(k%l, k, l, n)
            b[p[i]] -= A[p[i], k] * b[p[k]]
            count += 1
        end
    end
    return count
end

function backwardPivotal(A, b, p, n, l)
# - Backward Solve -
    count = 0
    x = zeros(n)
    for i in n:-1:1
        for j in (i+1):min(convert(Int, i+l*1.5+1), n)
            b[p[i]] -= A[p[i], j] * x[j]
            count += 1
        end
        x[i] = b[p[i]] / A[p[i], i]
        # count += 1
    end
    return x, count
end


function getLimit(mod, i, l, n)
    if mod == l-1
        return min(i+1+l, n)
    elseif mod == 0
        return min(i+l, n)
    else
        return min(i-mod+l, n)
    end
end

function getNonZeroElementsIndexes(i, n, l)
    imodl = i%l
    if imodl == 0
        limit = collect(max(1, i-l-1):i)
    else
        limit = collect(max(1, i-imodl-1):min(n, i+l-imodl))
    end
    if n >= i+l
        push!(limit, i+l)
    end
    return limit
end

if size(ARGS)[1] > 0
    sizeArg = parse(Int, ARGS[1])
    sparseA, N, l = readMatrixFromFile("data/$(sizeArg)_A.txt")
    b = readVectorFromFile("data/$(sizeArg)_B.txt")
else
    sparseA, N, l = readMatrixFromFile("data/1000_A.txt")
    b = readVectorFromFile("data/1000_B.txt")
end
# writeMatrixToFile(sparseA, "starting.csv")
tic()
x = pivotalGaussianElimination(sparseA, b, N, l)
toc()
# writeMatrixToFile(sparseA, "comp.csv")
println(x[1:11])
println(x[size(x)[1]:-1:size(x)[1]-10])
