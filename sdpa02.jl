"""
function sdpa_read
    (mDim, nBlock, blockStruct, c, F0, F) = sdpa_read(filename)
Read SDPA sparse format file

    primal: min  \sum_{i=1}^m c_i x_i
            s.t. X = \sum_{i=1}^m F_i x_i - F_0 \succeq O

    dual  : max  <F_0, Y>
            s.t. <F_i, Y> = c_i (i=1,..,m), Y \succeq O

"""
function sdpa_read(filename::String)
    if filename[end-4:end] != "dat-s"
        error("read_sdpa now supports only the files that end with dat-s.")
    end
    
    fio = open(filename,"r")
    lines = readlines(fio)
    close(fio)
    lineNo = 0
    # read comments
    while lineNo < length(lines)
        lineNo += 1
        str = lines[lineNo]
        if str[1] == '*' || str[1] =='"'
            print("Comment: Line $lineNo : $str")
            continue
        end
        break
    end
    # read m # this part should be better
    str = lines[lineNo]
    mDim = parse(Int, str[search(str, r"\ *[0-9]*")])
    println(" mDim = $mDim");
    lineNo += 1
    str = lines[lineNo]
    nBlock = parse(Int, str[search(str, r"\ *[0-9]*")])
    println(" nBlock = $nBlock")
    blockStruct = zeros(Int,nBlock)
    k = 1
    while k <= nBlock
        str = ""
        while length(str) == 0
            lineNo += 1
            str = lines[lineNo]
        end
        strIndex = 0
        while k <= nBlock && strIndex <= length(str)
            substr = str[strIndex+1:end]
            # @show substr
            searchIndex = search(substr, r"\ *[+-]*[0-9]*")
            # @show substr[searchIndex]
            blockStruct[k] = parse(Int,substr[searchIndex])
            k += 1
            strIndex += searchIndex[end]
        end
    end
    for k = 1:nBlock
        if blockStruct[k] == 1
            blockStruct[k] = -1
        end
    end
    println(" blockStruct = $blockStruct")

    
    lineNo += 1
    c = zeros(mDim)
    str = lines[lineNo]
    while str[1] == '*' || str[1] == '\n'
        lineNo += 1
        str = lines[LineNo]
    end
    str = replace(str, ',', ' ')
    str = replace(str, '{', ' ')
    str = replace(str, '}', ' ')
    strs = split(str, ' ', keep = false)
    # @show strs
    for k = 1:mDim
        # @show k, strs[k]
        c[k] = parse(Float64, strs[k])
    end
    # @show c
    
    tmpItemNo = length(lines) - lineNo
    item_ijkl = zeros(Int, tmpItemNo,4)
    item_v = zeros(Float64, tmpItemNo)
    item_line = zeros(Int, tmpItemNo)
    current_pos = 0
    itemNo = 0
    # @show itemNo
    while current_pos < tmpItemNo
        current_pos += 1
        str = lines[lineNo+current_pos]
        # @show str
        if str[1] == '*' || str[1] == '\n'
            # Note that ijkl remain zeros
            continue
        end
        str = replace(str, ',', ' ')
        str = replace(str, '{', ' ')
        str = replace(str, '}', ' ')
        itemNo += 1
        strIndex = 0
        for ijkl = 1:4
            substr = str[strIndex+1:end]
            searchIndex = search(substr, r"\ *[+-]*[0-9]*")
            item_ijkl[itemNo,ijkl] = parse(Int,substr[searchIndex])
            strIndex += searchIndex[end]
            # @show item_ijkl[itemNo,ijkl]
        end
        
        substr = str[strIndex+1:end]
        item_v[itemNo] = parse(Float64,substr)
        # @show item_v[itemNo]
        item_line[itemNo] = current_pos
    end

    println("Finished reading data, lines = $(lineNo+itemNo)")
    # checking the data consistency
    for current_pos = 1:itemNo
        k = item_ijkl[current_pos,1]
        l = item_ijkl[current_pos,2]
        i = item_ijkl[current_pos,3]
        j = item_ijkl[current_pos,4]
        v = item_v[current_pos]
        current_line = item_line[current_pos]+lineNo

        if l <= 0 ||  l > nBlock
            error("The block number l = $l is out of the range at line $current_line")
        end
        if k < 0 || k > mDim
            error("The constraint number k = $k is out of the range at line $current_line")
        end
        if blockStruct[l] < 0
            if i != j
                error("The LP part should have i = j, but i = $i and j = $j at line $current_line")
            end
        end
        if i < 0 || i > abs(blockStruct[l])
            error("The row index i = $i is out of the range at line $current_line")
        end
        if j < 0 || j > abs(blockStruct[l])
            error("The column index j = $j is out of the range at line $current_line")
        end
    end
    println("The input data is consistent with the diagonal block structure.")

    # @show item_ijkl
    # @show item_v

    if false
        NNZ = zeros(Int, m+1, nBlock)
        for current_pos = 1:itemNo
            k = item_ijkl[current_pos,1]
            l = item_ijkl[current_pos,2]
            # i = item_ijkl[current_pos,3]
            # j = item_ijkl[current_pos,4]
            NNZ[k+1,l] += 1
        end
    end
    
    F0 = Array{SparseMatrixCSC}(nBlock)
    F = Array{SparseMatrixCSC}(mDim, nBlock)
    for l = 1:nBlock
        if blockStruct[l] < 0
            F0[l] = spzeros(-blockStruct[l],1)
            for k = 1:mDim
                F[k,l] = spzeros(-blockStruct[l],1)
            end
        else
            F0[l] = spzeros(Float64, blockStruct[l],blockStruct[l])
            for k = 1:mDim
                F[k,l] = spzeros(Float64, blockStruct[l],blockStruct[l])
            end
        end
    end
    for current_pos = 1:itemNo
        k = item_ijkl[current_pos,1]
        l = item_ijkl[current_pos,2]
        i = item_ijkl[current_pos,3]
        j = item_ijkl[current_pos,4]
        v = item_v[current_pos]
        if abs(v) < 1.0e-10 # remove small elements
            continue
        end
        if k == 0
            if blockStruct[l] < 0
                F0[l][i] = v
            else
                F0[l][i,j] = v
            end
        else
            if blockStruct[l] < 0
                F[k,l][i] = v
            else
                F[k,l][i,j] = v
            end
        end
    end

    # symmetrization
    for l=1:nBlock
        if blockStruct[l] > 0
            F0[l] = F0[l] + F0[l] - diagm(diag(F0[l]))
            for k = 1:mDim
                F[k,l] = F[k,l] + F[k,l]' - diagm(diag(F[k,l]))
            end
        end
    end
    # return
    return mDim, nBlock, blockStruct, c, F0, F
end


function combine(A, nBlock, blockStruct)
	n = sum(abs.(blockStruct))
	Icomb = Array{Int64,1}(0)
	Jcomb = Array{Int64,1}(0)
	Vcomb = Array{Float64,1}(0)
	i=0
	for k=1:nBlock
		(I,J,V) = findnz(A[k])
		Icomb = append!(Icomb, I.+i)
		Jcomb = blockStruct[k]<0 ? append!(Jcomb, I.+i) : append!(Jcomb, J.+i)
		Vcomb = append!(Vcomb, V)
		i += blockStruct[k]
	end
	return sparse(Icomb,Jcomb,Vcomb,n,n)
end


function sdpa_newton(m, n, A, b, C)
	# parameters
	ε = 1.0e-6
	γ = 0.9
	β = 0.1
	λ = 1.0e+4
	maxIter = 50

	# initialization
	iter = 0
	X = λ*eye(n);
	Y = λ*eye(n);
	z = zeros(m,1);

	M = zeros(m,m); # Schur complement
	# main loop
	while iter < maxIter
		iter += 1
		# compute residual
		p = b - Float64[sum(A[:,:,k].*X) for k=1:m] # .* 要素ごとの積
		# p = b - Array{Float64}([sum(A[:,:,k].*X) for k=1:m])
		# float is necessary here
	    D = C - Y - sum([ A[:,:,k]*z[k] for k=1:m ])
		max_p = maximum(abs.(p)) # f.(x) = (f(x1), ... , f(xn)) 関数のベクトル化
		max_D = maximum(abs.(D))
		μ = bullet(X,Y)/n
		residuals = max(max_p, max_D, μ);
		if  residuals <= ε
			@printf "iter = %3d, r_p = %.2e, r_d = %.2e, μ = %.2e, α = %.2e\n" iter max_p max_D μ α
			println("Converge !!")
			break
		end
		
		Yinv = inv(Y)
		# build Schur complement
		for j=1:m
			XAY = X*A[:,:,j]*Yinv
			M[:,j] = [sum(A[:,:,i].*XAY) for i=1:m]
		end
		
		# search direction
		dz = M \ (p - Float64[ sum(A[:,:,k].* (β*μ*Yinv - X - X*D*Yinv)) for k=1:m ] )
		dY = D - sum([ A[:,:,k]*dz[k] for k=1:m ]);
		dX = β*μ*Yinv - X - X*dY*Yinv
		dX = (dX + dX')/2;

		# step length
	        Lx = ctranspose(chol(X))
	        Ly = ctranspose(chol(Y))
		
		(Dx, Vx) = eig(inv(Lx)*dX*inv(Lx)')
		eigx = minimum(Dx)
		αp = eigx < 0 ? -1.0/eigx : 100.0
		(Dy, Vy) = eig(inv(Ly)*dY*inv(Ly)')
		eigy = minimum(Dy)
		αd = eigy < 0 ? -1.0/eigy : 100.0
		
		α = min(γ*αp, γ*αd, 1.0)
		# update
		X += α*dX
		Y += α*dY
		z += α*dz
		@printf "iter = %3d, r_p = %.2e, r_d = %.2e, μ = %.2e, α = %.2e\n" iter max_p max_D μ α
	end

	if iter <= maxIter
		println("--Optimal Objective Value is ", bullet(C,X))
		#println("--Optimal Solution--")
		#println("X = ", X)
		#println("Y = ", Y)
		#println("z = ", z)
	else
		println("Failed to converge during ", iter, " iterations.")
	end

	return X, Y, z 
end

function sdpa_main(filename)
	mDim, nBlock, blockStruct, c, F0, F = sdpa_read(filename)

	m = mDim
	n = sum(abs.(blockStruct))
	
	C = combine(F0, nBlock, blockStruct)
	A = zeros(n,n,m)
	for k=1:m
		A[:,:,k] = combine(F[k,:], nBlock, blockStruct)
	end
	b = c

	sdpa_newton(m,n,A,b,C)
end
