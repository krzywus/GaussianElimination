
include("pivotGauss.jl")
include("fileOperations.jl")
include("basicGauss.jl")

using pivotGauss
using fileOperations
using basicGauss
# A[wiersz, kolumna]

function printUsage()
    println("\n  Usage: <method_name> <A_matrix_filename> [Optional]<b_vector_filename>")
    println("\t  <method_name> should be one of following:\n\t\t  -pivot/piv/p for gaussian elimination with partial pivoting.")
    println("\t\t  -basic/b for basic gaussian elimination.\n")
end

function checkArgs()
    if size(ARGS)[1] < 2
        printUsage()
        exit(0)
    end
    methods = ["pivot", "piv", "p", "basic", "b"]
    if !(ARGS[1] in methods)
        println("Wrong argument: method name must be one of $methods.")
        exit(0)
    end
    if !isfile(ARGS[2])
        println("Wrong argument: A matrix file does not exist.")
        exit(0)
    end
    if size(ARGS)[1] > 2
        if !isfile(ARGS[3])
            println("Wrong argument: B vector file does not exist.")
            exit(0)
        end
    end
end


checkArgs()

pivotAliases = ["pivot", "piv", "p"]
basicAliases = ["basic", "b"]

if ARGS[1] in pivotAliases
    println("Using pivotal gaussian elimination.")
    sparseA, n, l = readMatrixFromFile(ARGS[2])
    if size(ARGS)[1] > 2
        b = readVectorFromFile(ARGS[3])
        println("Starting gaussian elimination.")
        tic()
        x = pivotalGaussianElimination(sparseA, b, n, l)
        toc()

        writeVectorToFile(x, "x.txt", nothing)
        println(x[1:11])
        println(x[size(x)[1]:-1:size(x)[1]-10])
    else
        println("not implemented")
        x = ones(n)
    end
elseif ARGS[1] in basicAliases
    println("Using basic gaussian elimination.")
    sparseA, n, l = readMatrixFromFile(ARGS[2])
    if size(ARGS)[1] > 2
        b = readVectorFromFile(ARGS[3])
        println("Starting gaussian elimination.")
        tic()
        gaussElimination(sparseA, b, n, l)
        toc()

        writeVectorToFile(b, "x.txt", nothing)
        println(b[1:11])
        println(b[size(b)[1]:-1:size(b)[1]-10])
    else
        println("not implemented")
        x = ones(n)
    end
end
