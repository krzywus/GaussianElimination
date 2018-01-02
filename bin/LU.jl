# module gauss
include("fileOperations.jl")
using fileOperations

# export LUdecompose


function LUdecompose(A, n, l)
    # A[wiersz, kolumna]
    count = 0
    for i in 1:(n-1)
        for j in (i+1):getLimit(i%l, i, l, n) # utworzenie współczyników L i przechowanie poniżej przekątnej
            A[j,i] = A[j,i]/A[i,i]
            count += 1
        end
        for k in (i+1):min(n, i+l)         # odjęcie od równań równań pomnożonych przez współczynnik
            for j in (i+1):min(i+l+1, n)
                A[j,k] = A[j,k] - A[j,i] * A[i,k]
                count += 2
            end
        end
    end
    println("count: $count")
end


# calculates Y and stores it in vector B. diagonal element of L are always ones.
function solveLyb(A, b, n)
    count = 0
    for i in 2:n
        Left = 0
        for j in i-1:-1:max(1, i-l-1)
            if A[i, j] == 0
                continue
            end
            Left += A[i, j] * b[j]
            count += 2
        end
        b[i] = (b[i]-Left)   # /A[i,i],  /1.0
    end
    println("count: $count")
end

# calculates X and stores it in vector Y
function solveUxy(A, y, n)
    count = 0
    y[n] = y[n]/A[n,n]
    for i in (n-1):-1:1
        Left = 0
        for j in (i+1):min(i+l, n)
            Left += A[i, j] * y[j]
            count += 2
        end
        y[i] = (y[i]-Left)/A[i,i]
    end
    println("count: $count")
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

# sparseA, N, l = readMatrixFromFile("data/100_A.txt")
# b = readVectorFromFile("data/100_B.txt")
if size(ARGS)[1] < 2
    println("Not enough arguments")
    exit(1)
end
sparseA, N, l = readMatrixFromFile(ARGS[1])
b = readVectorFromFile(ARGS[2])

println("LUdecompose")
tic()
LUdecompose(sparseA, N, l)
toc()
println("solveLyb")
tic()
solveLyb(sparseA, b, N)
toc()
println("solveUxy")
tic()
solveUxy(sparseA, b, N)
toc()
println(b[1:11])
println(b[size(b)[1]:-1:size(b)[1]-10])

# println(sparseA)
