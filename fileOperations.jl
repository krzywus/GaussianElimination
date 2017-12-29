module fileOperations

export readMatrixFromFile, readVectorFromFile, writeMatrixToFile, writeVectorToFile

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

function writeVectorToFile(vector, filename, error)
    println("Writing file $filename.")
    open(filename, "w") do file
        if error != nothing
            write(file, error)
        end
        for x in vector
            write(file, "$x\n")
        end
    end
    println("Done writing file $filename.")
end

end #module
