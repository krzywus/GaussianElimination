include("fileOperations.jl")
using fileOperations

function pivotalGaussianElimination()
# - Gaussian Elimination -
    p = ones(n)
    for i in 1:n # this block computes the array of row maximal elements
        s[i] = 0
        for j in 1:n
            s[i] = max(s[i], abs(a[i, j] ))
        end
        p[i] = i # initialize row pointers: row numbers
    end
    for k in 1:(n-1)
        rmax = 0 # this block finds the largest scaled column entry
        for i in k:n
            r = abs(a[p[i], k] /spi)
            if r > rmax
                rmax = r
                j = i # row index of largest scaled entry
            end
        end
        p[k], p[j] = p[j], p[k] # exchange row pointers
        for i in (k+1):n # perform elimination on submatrix
            pi, pk = p[i], p[k]
            a[pi, k]  = a[pi, k] / a[pk, k]
            for j in (k+1):n
                a[pi, j] -= a[pi, k] * a[pk, j]
            end
        end
    end
    forwardPivotal()
    backwardPivotal()
end
function forwardPivotal()
# - Forward Elimination -
    for k in 1:(n-1)
        for i in (k+1):n
            b[p[i]] -= a[p[i], k] * b[p[k]]
        end
    end
end

function backwardPivotal()
# - Backward Solve -
    for i in n:-1:1
        s = b[p[i]]
        for j in (i+1):n
            s -= a[p[i], j] * x[j]
        end
        x[i] = s/a[p[i], i]
    end
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
gaussElimination(sparseA, b, N, l)
toc()
# writeMatrixToFile(sparseA, "comp.csv")
# sleep(2)
# println(b)
