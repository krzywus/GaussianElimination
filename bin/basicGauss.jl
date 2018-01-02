module basicGauss

export gaussElimination

function gaussElimination(A, b, n, l)
    count = 0
    for k in 1:(n-1)
        for i in (k+1):getLimit(k%l, k, l, n)
            A[i,k] = A[i,k]/A[k,k]
            count += 1
            for j in (k+1):n
                if A[k,j] != 0
                    A[i,j] -= A[i,k]*A[k,j]
                    count += 2
                end
            end
        end
    end
    # println("count:$count")
    count += forward(A, b, n, l)
    # println("count:$count")
    count += backward(A, b, n, l)
    println("Operations count: $count")
end

function forward(A, b, n, l)
    count = 0
    for i in 1:(n-1)
        for j in (i+1):getLimit(i%l, i, l, n)
            b[j] -= A[j,i] * b[i]   # b[i] = x[i], mnożenie pozostałych przez ten współczynnik
            count += 2
        end
    end
    return count
end

function backward(A, b, n, l)
    count = 0
    for i in n:-1:1
        for j in (i+1):min(i+l, n)
            b[i] -= A[i,j]*b[j]
            count += 2
        end
        b[i] = b[i]/A[i,i]
        count += 1
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

end
