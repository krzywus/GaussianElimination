
include("pivotGauss.jl")
include("fileOperations.jl")
include("basicGauss.jl")
include("LU.jl")

using fileOperations
using pivotGauss
using basicGauss
using gaussLU
# A[wiersz, kolumna]

function printUsage()
    println("\n  Usage: <method_name> <A_matrix_filename> [Optional]<b_vector_filename>")
    println("\t  <method_name> should be one of following:\n")
    println("\t\t  -pivot/piv/p for gaussian elimination with partial pivoting,\n")
    println("\t\t  -basic/b for basic gaussian elimination,\n")
    println("\t\t  -lu/LU for basic LU decomposition.\n")
end

function checkArgs()
    if size(ARGS)[1] < 2
        printUsage()
        exit(0)
    end
    methods = ["pivot", "piv", "p", "basic", "b", "lu", "LU"]
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

function getError(x, countError)
    if countError
        max_error = 0
        # max_diff = 0
        for X in x
            error = (1.0-X)/1.0
            # if(abs(1.0-X) > abs(max_diff))
            #     max_diff = X
            # end
            if(abs(error) > abs(max_error))
                max_error = error
            end
        end
        return max_error
    else
        return nothing
    end
end


checkArgs()

pivotAliases = ["pivot", "piv", "p"]
basicAliases = ["basic", "b"]
LUAliases = ["lu", "LU"]
xFilepath = "output/x.txt"

if ARGS[1] in pivotAliases
    println("Using pivotal gaussian elimination.")
    A, n, l = readMatrixFromFile(ARGS[2])
    if size(ARGS)[1] > 2
        b = readVectorFromFile(ARGS[3])
        countError = false
    else
        println("Starting calculating vector b from Ax=b.")
        tic()
        b = calculateVectorFromSolution(A, ones(n), n, l)
        toc()
        countError = true
    end
    println("Starting gaussian elimination.")
    tic()
    x = pivotalGaussianElimination(A, b, n, l)
    toc()

    mistake = getError(x, countError)
    writeVectorToFile(x, xFilepath, mistake)
    println(x[1:11])
    println(x[size(x)[1]:-1:size(x)[1]-10])
elseif ARGS[1] in basicAliases
    println("Using basic gaussian elimination.")
    A, n, l = readMatrixFromFile(ARGS[2])
    if size(ARGS)[1] > 2
        b = readVectorFromFile(ARGS[3])
        countError = false
    else
        println("Starting calculating vector b from Ax=b.")
        tic()
        b = calculateVectorFromSolution(A, ones(n), n, l)
        toc()
        countError = true
    end
    println("Starting gaussian elimination.")
    tic()
    gaussElimination(A, b, n, l)
    toc()

    mistake = getError(b, countError)
    writeVectorToFile(b, xFilepath, mistake)
    println(b[1:11])
    println(b[size(b)[1]:-1:size(b)[1]-10])
elseif ARGS[1] in LUAliases
    println("Using basic LU decomposition.")
    A, n, l = readMatrixFromFile(ARGS[2])
    if size(ARGS)[1] > 2
        b = readVectorFromFile(ARGS[3])
        countError = false
    else
        println("Starting calculating vector b from Ax=b.")
        tic()
        b = calculateVectorFromSolution(A, ones(n), n, l)
        toc()
        countError = true
    end

    println("Starting basic LU decompostion.")
    tic()
    LUdecompose(A, n, l)
    toc()
    println("Started solving Ly=b")
    tic()
    solveLyb(A, b, n, l)
    toc()
    println("Started solving Ux=y")
    tic()
    solveUxy(A, b, n, l)
    toc()

    mistake = getError(b, countError)
    writeVectorToFile(b, xFilepath, mistake)
    println(b[1:11])
    println(b[size(b)[1]:-1:size(b)[1]-10])
end
