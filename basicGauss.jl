function readMatrixFromFile(filename::String)
    println("Reading file $filename.")
    open(filename) do f
        lines = readlines(f)
        firstLine = lines[1]
        deleteat!(lines, 1)
        I = zeros(0)
        J = zeros(0)
        V = zeros(0)
#         k = 1 # debug
        n = parse(Int64, split(firstLine, r"\ |\t")[1])
        l = parse(Int64, split(firstLine, r"\ |\t")[2])  # debug
#         Q = zeros(Float64, n, n) # debug
        for line in lines
            splitLine = split(line, r"\ |\t")
#             println(splitLine[1], "|", splitLine[2], "|", splitLine[3])
            push!(I, parse(Int64, splitLine[1]))
            push!(J, parse(Int64, splitLine[2]))
            push!(V, parse(Float64, splitLine[3]))
#             println(I[k], "|", J[k], "|", V[k])
#             Q[Int(I[k]), Int(J[k])] = V[k]
#             k += 1
           end
        println("Done reading file $filename.")
        return sparse(I,J,V), n, l
    end
end


function readVectorFromFile(filename::String)
    println("Reading file $filename.")
    open(filename) do f
        lines = readlines(f)
        firstLine = lines[1]
        deleteat!(lines, 1)
        V = zeros(0)
#         n = parse(Int64, split(firstLine, r"\ |\t")[1])
        for line in lines
            push!(V, parse(Float64, line))
           end
        println("Done reading file $filename.")
        return V
    end
end

function writeMatrixToFile(matrix, filename, delim="\t")
    println("Writing file $filename.")
    open(filename, "w") do file
          writedlm(file, matrix, delim)
    end
    println("Done writing file $filename.")
end

function forward(A, b, n)
# — Forward Elimination —
    for i in 1:(n-1)
        for j in (i+1):n
            if b[i] != 0
                b[j] -= A[j,i] * b[i]
            end
        end
    end
end

function backward(A, b, n)
# remove x
#     — Backward Solve —
    x = zeros(n)
    for i in n:-1:1
        s = b[i]
        for j in (i+1):n
            if x[j] != 0
                s -= A[i,j]*x[j]
            end
        end
        x[i] = s/A[i,i]
    end
    return x
end

function gaussElimination(A, b, n, l)
# A[wiersz, kolumna]
#     — Gaussian Elimination —
    for k in 1:(n-1) # iterator kolumn
        println("gauss k: $k/$n")
        for i in (k+1):n # iterator wierszy
            # println("gauss i: $i/$n")
            if A[i, k] != 0
                A[i,k] = A[i,k]/A[k,k]
                for j in (k+1):n
                    if A[k,j] != 0
                        A[i,j] -= A[i,k]*A[k,j]
                    end
                end
            end
        end
    end
    forward(A, b, n)
    x = backward(A, b, n)
    return x
end
#
# A = [2.0    -2.0   0.0;
#      -2    0    2;
#      0   -2   0]
# N=3
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
tic()
x = gaussElimination(sparseA, b, N, l)
toc()
writeMatrixToFile(sparseA, "comp.csv")
println(x)
