#=======================

=======================#

function bullet(U::Array, V::Array)
    sum(U.*V)
end

function isConverged(ϵ, x, X, Y)
    residuals = 1
    if residuals <= ϵ && bullet(X,Y)/n <= ϵ # μ = bullet(X,Y)/n 
        return true
    else
        return false
    end
end

function update(m, n, c, F0, F, β, x, X, Y)

    #Step2: Search Direction
    B = zeros(m,m)
    invX = inv(X)
    for q=1:m
        B[:,q] = [bullet(invX*F[p]*Y, F[q]) for p=1:m]
    end

    μ = bullet(X,Y)/n
    R_c = β*μ*eye(n) - X*Y
    R_p = sum(F[k]*x[k] for k=1:m) - F0 - X
    d = c - [bullet(F[k],Y) for k=1:m]
    r = -d + [bullet(F[k], inv(X)*(R_c-R_p*Y)) for k=1:m]

    dx = inv(B)*r
    dX = R_p + sum(F[k]*dx[k] for k=1:m)
    dYtil = inv(X)*(R_c-dX*Y)
    dY = (dYtil + dYtil')/2

    #Step3: Step Length
    #=
    X+α*dX => 0 を満たすαを求める
    X = LL' とCholesky分解する。Lは下三角行列
    LL'+α*dX => 0
    I + α*inv(L)*dX*inv(L') => 0
    inv(L)*dX*inv(L') = V*D*V' と固有値分解する。D:固有値対角行列、V:固有値ベクトル行列
    I + α*V*D*V' => 0
    V'*V + α*D => 0 (V*V' = I を利用)
    I + α*D => 0
    I + α*D は(i,i)成分が 1+α*λ_i の対角行列
    任意のvに対し v*(I + α*D)*v => 0 が成立
    任意のiに対し (1+α*λ_i)*v*v => 0
    =#
    α_p = 0.9
    α_d = 0.9

    #Step4: Update
    x += α_p*dx
    X += α_p*dX
    Y += α_d*dY
    println(x)
    println(X)
    println(Y)
    return (x, X, Y)
end

function sdpa_newton(m, n, c, F0, F)
    # parameters
    ϵ = 1.0e-6
    γ = 0.9
    β = 0.1
    λ = 1.0e+4
    maxIter = 30

    #Step0: Initialization
    iter = 0
    X = λ*eye(n)
    Y = λ*eye(n)
    x = zeros(m)
    B = zeros(m,m) #schur complement matrix(SCM)
    μ = bullet(X,Y)/n

    while iter <= maxIter
        iter += 1
        if isConverged(ϵ, x, X, Y)
            println("Converged!")
            break
        end
        (x, X, Y) = update(m, n, c, F0, F, β, x, X, Y)
    end
    return (x, X, Y)
end




function sdpa_newton2(m, n, c, F0, F)
    #=
        We choose an initial point (x_0,X_0,Y_0), with X_0,
        Y_0 to be definite.
        Let μ_0 = (X_0 bullet Y_0)/n and h = 0.
        We set the parameters 0<β<1, 0<γ<1
    Step1: Checking Convergence
        if μ_h is suffeciently small and (x_h,X_h,Y_h)
        approximately satisfies the feasibility,
        we stop the iteration and output the current point
        (x_h,X_h,Y_h) as an approximate optimal solution.
    Step2: Search Direction
       We compute a search direction (dx,dX,dY) toward a
        target point (x(β*μ^h),X(β*μ^h),Y(β*μ^h))on the
        central path with μ = β*μ^h.
    Step3: Step Length
        We compute step length α_p and α_d such that
        X^h + α_p*dX is semidefinite and
        Y^h + α_d*dY is semidefinite.
    Step4: Update
        We update the current point by 
        (x^(h+1),X^(h+1),Y^h+1) = 
            (x^h+γ*α_p*dx,γ*α_p*dX,Y*α_d*dX).
        Let μ^(h+1) = (X^(h+1) bullet Y^(h+1))/n
        and h <- h+1. Go to Step1.
    =#

    # parameters
    ε = 1.0e-6
    γ = 0.9
    β = 0.1
    λ = 1.0e+4
    maxIter = 30

    #Step0: Initialization
    iter = 0
    X = λ*eye(n)
    Y = λ*eye(n)
    x = zeros(m)
    B = zeros(m,m) #schur complement matrix(SCM)
    μ = bullet(X,Y)/n

    while iter <= maxIter
        iter += 1
        
        #Step1: Checking Convergence
        residuals = 0
        if residuals <= ϵ
            break
        end

        #Step2: Search Direction
        invX = inv(X)
        for q=1:m
            B[:,q] = [bullet(invX*F[:;:;p]*Y, F[:;:;q]) for p=1:m]
        end

        R_c = β*μ*eye(n) - X*Y
        R_p = sum(F[:;:;k]*x for k=1:m) - F_0 - X
        d = c - [bullet(F[:;:;k],Y) for k=1:m]
        r = -d + [bullet(F[:;:;k], inv(X)*(R_c-R_p*Y)) for k=1:m]

        dx = inv(B)*r
        dX = R_p + sum(F[:;:;k]*dx[k] for k=1:m)
        dYtil = inv(X)*(R_c-dX*Y)
        dY = (dYtil + dYtil')/2

        #Step3: Step Length
        #=
        X+α*dX => 0 を満たすαを求める
        X = LL' とCholesky分解する。Lは下三角行列
        LL'+α*dX => 0
        I + α*inv(L)*dX*inv(L') => 0
        inv(L)*dX*inv(L') = V*D*V' と固有値分解する。D:固有値対角行列、V:固有値ベクトル行列
        I + α*V*D*V' => 0
        V'*V + α*D => 0 (V*V' = I を利用)
        I + α*D => 0
        I + α*D は(i,i)成分が 1+α*λ_i の対角行列
        任意のvに対し v*(I + α*D)*v => 0 が成立
        任意のiに対し (1+α*λ_i)*v*v => 0

        =#

        #Step4: Update
        x += α_p*dx
        X += α_p*dX
        Y += α_d*dY
    end

    return (x, X, Y, )
end