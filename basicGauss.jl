include("fileOperations.jl")
using fileOperations

function gaussElimination(A, b, n, l)
# A[wiersz, kolumna]
#     — Gaussian Elimination —
    count = 0
    for k in 1:(n-1) # iterator kolumn
        # println("gauss k: $k/$n")
        for i in (k+1):getLimit(k%l, k, l, n) # iterator wierszy, zaczynając od 'poniżej przekątnej' do /końca/wiadomych zer
            # if A[i, k] == 0
            #     println("CAUGHT A ZERO")
            # end
            A[i,k] = A[i,k]/A[k,k]
            count += 1
            for j in (k+1):n
                if A[k,j] != 0
                    A[i,j] -= A[i,k]*A[k,j]
                    count += 1
                end
            end
        end
    end
    println("count:$count")
    count += forward(A, b, n, l)
    println("count:$count")
    count += backward(A, b, n, l)
    println("count:$count")
end

function forward(A, b, n, l)
# — Forward Elimination —
    count = 0
    for i in 1:(n-1)
        # println("forward i: $i/$n")
        for j in (i+1):getLimit(i%l, i, l, n)
            b[j] -= A[j,i] * b[i]   # b[i] = x[i], mnożenie pozostałych przez ten współczynnik
            count += 1
        end
    end
    return count
end

function backward(A, b, n, l)
#     — Backward Solve —
    count = 0
    for i in n:-1:1
        # println("backward u: $(n-i)/$n")
        for j in (i+1):min(i+l, n)
            # if b[j] != 0 && A[i, j] != 0
                b[i] -= A[i,j]*b[j]
                count += 1
            # end
        end
        b[i] = b[i]/A[i,i]
    end
    return count
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

# A = [2.0    -2.0   0.0;
#      -2    0    2;
#      0   -2   0]
# N=3
# l=1
# b = [6, 0, -2]
# sparseA = sparse(copy(A))
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
