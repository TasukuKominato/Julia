function sdpa_newton(m, nBlock, blockStruct, b, A0, A)

    start_All = now()
    time_All = 0
    time_Schur = 0
    time_Chol = 0

    # parameters
    β = 0.1
    ε = 1.0e-6
    γ = 0.9
    λ = 1.0e+4
    maxIter = 100

    n = sum(abs.(blockStruct))

    # initialization
    iter = 0
    X = Array{Array{Float64,2}}(nBlock)
    Z = Array{Array{Float64,2}}(nBlock)
    y = zeros(m,1)
    B = zeros(m,m) # Schur complement
    dX = Array{Array{Float64,2}}(nBlock)

    for l = 1:nBlock
	if blockStruct[l]<0
	    X[l] = λ*(ones(-blockStruct[l],1))
	else
	    X[l] = λ*eye(blockStruct[l])
	end
	
	if blockStruct[l]<0
	    Z[l] = λ*(ones(-blockStruct[l],1))
	else
	    Z[l] = λ*eye(blockStruct[l])
	end
    end

    @printf "iter time     r_p     r_d      mu  alphaP  alphaD     objP     objD\n"

    # main loop
    while iter < maxIter
		iter += 1
		# compute residual
		AX = zeros(m)

		for i=1:m
		    for l=1:nBlock
			if blockStruct[l]<0
			    AX[i] += sum(A[i,l].*X[l])
			else
			    AX[i] += sum(A[i,l].*X[l])
			end
		    end
		end

		rp = b - AX
		Rd = [A0[l] - sum([A[k,l]*y[k] for k=1:m]) - Z[l] for l=1:nBlock]
		max_rp = maximum(abs.(rp))
		max_Rd = maximum(maximum(abs.(Rd[l])) for l=1:nBlock)
		μ = sum(sum(X[l].*Z[l]) for l=1:nBlock)/n
		residuals = max(max_rp, max_Rd, μ)
		if  residuals <= ε
		    #@printf "iter = %3d, r_p = %.2e, r_d = %.2e, μ = %.2e, α = %.2e\n" iter max_rp max_Rd μ α
		    println("Converge !!")
		    break
		end

		# build Schur complement B
		invX = Array{Array}(nBlock)
		for l=1:nBlock
		    if blockStruct[l]<0
			invX[l] = inv.(X[l])
		    else
			invX[l] = inv(X[l])
		    end
		end

		invZ = Array{Array}(nBlock)
		for l=1:nBlock
		    if blockStruct[l]<0
			invZ[l] = inv.(Z[l])
		    else
			invZ[l] = inv(Z[l])
		    end
		end
	        isSparse = true
	        start_Schur = now()

	        if false
	            B = zeros(m,m)
	        else
	        B = SharedArray{Float64}(m,m)
	        @sync @parallel for j=1:m
	            for i=1:m
	                B[i,j] = 0
	            end
	        end
	        end
	        # 列ごとに計算
	        if isSparse
	            for l=1:nBlock
	                if blockStruct[l] > 0
	                    @sync @parallel for j=1:m
	                        (R2, S2, W2) = findnz(A[j,l])
	                        for i=1:j
	                            (P1, Q1, V1) = findnz(A[i,l])
	                            sum1 = 0.0
	                            for i1 = 1:length(P1)
	                                for i2 = 1:length(R2)
	                                    sum1 += V1[i1]*X[l][P1[i1],R2[i2]]*W2[i2]*invZ[l][S2[i2],Q1[i1]]
	                                end # i2
	                            end # i1
	                            B[i,j] += sum1
	                        end # i
	                    end # j
	                else # blockStruct[l] < 0
	                    @sync @parallel for j=1:m
			        		G = X[l].*A[j,l].*invZ[l]
			        		for i=1:j
			            		B[i,j] += sum(A[i,l].*G)
			        		end
	                    end
	                end
	            end
	        else # dense
		    for l=1:nBlock
		        for j=1:m
			    G = blockStruct[l]<0 ? X[l].*A[j,l].*invZ[l] : X[l]*A[j,l]*invZ[l]
			    for i=1:j
			        B[i,j] += sum(A[i,l].*G)
			    end
		        end
		    end
	        end

	        B = B + B' - spdiagm(diag(B))
	        time_Schur += Dates.value(now()-start_Schur)

		Rc = Array{Array{Float64,2}}(nBlock)
		for l=1:nBlock
		    if blockStruct[l]<0
			Rc[l] = β*μ*ones(-blockStruct[l]) - X[l].*Z[l]
		    else
			Rc[l] = β*μ*eye(blockStruct[l]) - X[l]*Z[l]
		    end
		end
		
		W = Array{Array}(nBlock)
		for l=1:nBlock
		    if blockStruct[l]<0
			W[l] = (Rc[l]-X[l].*Rd[l]).*invZ[l]
		    else
			W[l] = (Rc[l]-X[l]*Rd[l])*invZ[l]
		    end
		end

		s = [rp[i] - sum(sum(A[i,l].*W[l]) for l=1:nBlock) for i=1:m]

	        # B = B + (1.0e-3)*speye(m)
	        start_Chol = now()
		dy = B\s
	        time_Chol = Dates.value(now()-start_Chol)

		dZ = [Rd[l] - sum(A[i,l]*dy[i] for i=1:m) for l=1:nBlock]
		
		for l=1:nBlock
		    if blockStruct[l]<0
			# println(size(X[l]))
			dX[l] = (Rc[l] - X[l].*dZ[l]).*invZ[l]
		    else
			dX[l] = (Rc[l] - X[l]*dZ[l])*invZ[l]
		        dX[l] = (dX[l] + dX[l]')/2
		    end
		end


	    # predictor-corrector
	    if false

		Rc = Array{Array{Float64,2}}(nBlock)
		for l=1:nBlock
		    if blockStruct[l]<0
			Rc[l] = β*μ*ones(-blockStruct[l]) - X[l].*Z[l] -dX[l].*Z[l]
		    else
			Rc[l] = β*μ*eye(blockStruct[l]) - X[l]*Z[l] - dX[l]*Z[l]
		    end
		end
		
		W = Array{Array}(nBlock)
		for l=1:nBlock
		    if blockStruct[l]<0
			W[l] = (Rc[l]-X[l].*Rd[l]).*invZ[l]
		    else
			W[l] = (Rc[l]-X[l]*Rd[l])*invZ[l]
		    end
		end

		s = [rp[i] - sum(sum(A[i,l].*W[l]) for l=1:nBlock) for i=1:m]

	        # B = B + (1.0e-3)*speye(m)
	        start_Chol = now()
		dy = B\s
	        time_Chol = Dates.value(now()-start_Chol)

		dZ = [Rd[l] - sum(A[i,l]*dy[i] for i=1:m) for l=1:nBlock]
		
		for l=1:nBlock
		    if blockStruct[l]<0
			# println(size(X[l]))
			dX[l] = (Rc[l] - X[l].*dZ[l]).*invZ[l]
		    else
			dX[l] = (Rc[l] - X[l]*dZ[l])*invZ[l]
		        dX[l] = (dX[l] + dX[l]')/2
		    end
		end
		end 
		# end of predictor-corrector
	        

        #--
        # The following values should be 0
        # execpt SS (blockStuct[l]>0)
        #--

        if false
            s0 = [rp[i] - sum(sum(A[i,l].*dX[l]) for l=1:nBlock) for i=1:m]
            @show norm(s0)
	    
	    S0 = Rd
            for l=1:nBlock
                S0[l] -= dZ[l]
                for k=1:m
                    S0[l] -= A[k,l]*dy[k]
                end
                @show norm(S0[l])
            end
            for l=1:nBlock
	        if blockStruct[l]<0
                    SS = X[l].*dZ[l] + dX[l].*Z[l] - Rc[l]
	        else
                    SS = X[l]*dZ[l] + dX[l]*Z[l] - Rc[l]
	        end
                @show norm(SS)
            end
        end
        
        #--

		# step length
		αp = 1.0
		αd = 1.0
		for l=1:nBlock
		    if blockStruct[l]<0
			eigx = minimum(X[l])
			eigz = minimum(Z[l])
		    else
		    # chol(X) = U, X = U'U: the cholesky factorization of X, U: UpperTriangular matrix
		    # Lx = U': LowerTriangular matrix
			Lx = ctranspose(chol(X[l]+(1.0e-5)*speye(blockStruct[l])))
			Lz = ctranspose(chol(Z[l]+(1.0e-5)*speye(blockStruct[l])))

			LXL = inv(Lx)*dX[l]*inv(Lx)'
			LZL = inv(Lz)*dZ[l]*inv(Lz)'

			(Dx, Vx) = eig((LXL + LXL')/2)
			(Dz, Vz) = eig((LZL + LZL')/2)
			eigx = minimum(Dx)
			eigz = minimum(Dz)
		    end
		    αpbar = eigx < 0 ? -1.0/eigx : 100.0
		    αdbar = eigz < 0 ? -1.0/eigz : 100.0
		    αp = min(αp, γ*αpbar)
		    αd = min(αd, γ*αdbar)
		end

        # mu を 0 に近づけるよりも早く、実行可能性を達成するようにする
        if max_rp > 1.0e-6 || max_Rd > 1.0e-6

            xMatvMat = sum(sum(X[l].*dZ[l]) for l=1:nBlock)/n
            uMatzMat = sum(sum(dX[l].*Z[l]) for l=1:nBlock)/n
            uMatvMat = sum(sum(dX[l].*dZ[l]) for l=1:nBlock)/n
            μNew = μ + αp*uMatzMat + αd*xMatvMat + αp*αd*uMatvMat
            thetaMax = max(1.0-αp, 1.0-αd)


            while μNew / μ < thetaMax
                # @show μNew, μ
                αMax = 0.95*max(αp, αd)
                αp = min(αp, αMax)
                αd = min(αd, αMax)
                if αp < 1.0e-6 && αd < 1.0e-6
                    break
                end
                μNew = μ + αp*uMatzMat + αd*xMatvMat + αp*αd*uMatvMat
                thetaMax = max(1.0-αp, 1.0-αd)
            end
        end
	    
		# update
		X += αp*dX
		y += αd*dy
		Z += αd*dZ
		α = min(αp, αd)

        if α < 1.0e-6
            println("Step length is too small")
            break
        end

		obj_p = sum(sum(A0[l].*X[l]) for l=1:nBlock)
		obj_d = b'y[:,1]

		@printf "%4d %4.1f %2.1e %2.1e %2.1e %2.1e %2.1e %+2.1e %+2.1e\n" iter Dates.value(now()-start_All)/1000 max_rp max_Rd μ αp αd obj_p obj_d
		# @printf " %d %f" iter Dates.value(now()-start_All)/1000
		# @printf "iter = %3d, time = %.3f, r_p = %.2e, r_d = %.2e, μ = %.2e, αp = %.2e, αd = %.2e, obj_p = %.2e, obj_d = %.2e\n" iter Dates.value(now()-start_All)/1000 max_rp max_Rd μ αp αd obj_p obj_d
    end
	time_All = Dates.value(now()-start_All)
	result = (nprocs(), nworkers(), time_All, time_Schur)
	result
	# @show time_All, time_Schur, time_Chol

end
