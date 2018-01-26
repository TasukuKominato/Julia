function print_tree(typ)
    println(typ)
    print_tree(typ, [])
end

function print_tree(typ, level)
    stypes = subtypes(typ)
    for stype in stypes
        if stype !== stypes[end]
            println(join(level) * "├─── $stype")
            push!(level, "|    ")
        else
            println(join(level) * "└─── $stype")
            push!(level, "     ")
        end
        print_tree(stype, level)
        pop!(level)
    end
end