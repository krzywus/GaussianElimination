# module gauss

# export readMatrixFromFile, writeMatrixToFile

# funkcja do zczytywania macierzy A
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


function LUdecompose(A, n)
    # A[wiersz, kolumna]
    for i in 1:(n-1)  # iterator kolumn
        println("LUdecompose i: $i")
        for j in (i+1):n # utworzenie współczyników L i przechowanie poniżej przekątnej
            if A[j, i] == 0
                continue
            end
            A[j,i] = A[j,i]/A[i,i]
        end
        for k in (i+1):n         # odjęcie od równań równań pomnożonych przez współczynnik
            for j in (i+1):n
                A[j,k] = A[j,k] - A[j,i] * A[i,k]
            end
        end
    end
end


# calculates Y and stores it in vector B. diagonal element of L are always ones.
function solveLyb(A, b, n)
    for i in 2:n
        println("solveLyb i: $i")
        Left = 0
        for j in i-1:-1:1
            Left += i==j ? b[j] : A[i, j] * b[j]
        end
        b[i] = (b[i]-Left)   # /A[i,i],  /1.0
    end
end

# calculates X and stores it in vector Y
function solveUxy(A, y, n)
    y[n] = y[n]/A[n,n]
    for i in (n-1):-1:1
        println("solveUxy i: $i")
        Left = 0
        for j in (i+1):n
            Left += A[i, j] * y[j]
        end
        y[i] = (y[i]-Left)/A[i,i]
    end
end
# end

#             if floor(j-1/n) == floor(i/j) # jesteśmy w A
#                 A[j,i] = A[j,i]/A[i,i]
#                continue
#             end


# sparseA, N, l = readMatrixFromFile("data/100_A.txt")
# b = readVectorFromFile("data/100_B.txt")
sparseA, N, l = readMatrixFromFile("data/medA.txt")
b = readVectorFromFile("data/medB.txt")

tic()
# writeMatrixToFile(sparseA, "starting.csv")
LUdecompose(sparseA, N)
# writeMatrixToFile(sparseA, "comp.csv")
solveLyb(sparseA, b, N)
# writeMatrixToFile(b, "calc1.csv")
# println("y: $b")
solveUxy(sparseA, b, N)
# writeMatrixToFile(b, "calc2.csv")
toc()
println("x: $b")

# println(sparseA)
