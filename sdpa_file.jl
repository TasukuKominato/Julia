"""
function sdpa_read
    (mDim, nBlock, blockStruct, c, F0, F) = read_sdpa(filename)
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
    
    # X = Array{Array{Float64,2}}(nBlock)
    F0 = Array{SparseMatrixCSC}(nBlock)
    F = Array{SparseMatrixCSC}(mDim, nBlock)
    for l = 1:nBlock
        if blockStruct[l] < 0
            F0[l] = spzeros(-blockStruct[l])
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


"""
function sedumi2sdpa
    (mDim, nBlock, blockStruct, c2, F0, F) = sedumi2sdpa(At, b, c, K_l, K_s)

Convert SeDuMi format to SDPA format
Note that K.l and K.s should be input separately as K_l::Int and K_s::Array{Int}
    The conversion is as follow
    -c   -> F0
    -A_i -> F_i
    -b_i -> c_i
"""
function sedumi2sdpa(At::SparseMatrixCSC, b::Array, c::SparseMatrixCSC, K_l::Int, K_s::Array{Int})
    #  b = float(b) # b can be Array{Int} and Array{Float64}
    (n, m) = size(At)
    mDim = m
    nBlock = length(K_s)
    blockStruct = K_s
    if size(blockStruct,1) < size(blockStruct,2)
        # convert to a vertical vector
        blockStruct = blockStruct'
    end
    if K_l > 0
        nBlock += 1
        blockStruct = [-K_l; K_s]
    end
    F0 = Array{SparseMatrixCSC}(nBlock)
    F = Array{SparseMatrixCSC}(mDim, nBlock)
    for l = 1:nBlock
        if blockStruct[l] < 0
            F0[l] = spzeros(-blockStruct[l])
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

    si = 0 # startIndex
    for l = 1:nBlock
        bs = blockStruct[l] 
        if bs < 0
            F0[l] = -c[si+1:si+abs(bs)]
            si += abs(bs)
        else
            F0[l] = -sparse(reshape(c[si+1:si+bs*bs], bs, bs))
            si += bs*bs
        end
    end

    
    for k = 1:m
        si = 0 # startIndex
        for l = 1:nBlock
            bs = blockStruct[l] 
            if bs < 0
                F[k,l] = -At[si+1:si+abs(bs),k]
                si += abs(bs)
            else
                F[k,l] = -sparse(reshape(At[si+1:si+bs*bs,k], bs, bs))
                si += bs*bs
            end
        end
    end
    c2 = -b
    return (mDim, nBlock, blockStruct, c2, F0, F)
end

