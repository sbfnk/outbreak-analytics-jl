# Pretty printing for ContactMatrix

function Base.show(io::IO, ::MIME"text/plain", cm::ContactMatrix{T, N}) where {T, N}
    n = ndimgroups(cm.groupings)
    dimstr = join(size(cm), "×")
    println(io, "ContactMatrix{$T} $dimstr (setting: $(cm.setting))")

    if n == 1
        _show_2d(io, cm)
    else
        names_str = join(cm.groupings.names, " × ")
        println(io, "  Groupings: $names_str")
        for (nm, lb) in zip(cm.groupings.names, cm.groupings.labels)
            println(io, "    $nm: ", join(lb, ", "))
        end
        println(io, "  Data: ", summary(cm.data))
    end
end

function _show_2d(io::IO, cm::ContactMatrix)
    labels = cm.groupings.labels[1]
    println(io, "  Grouping: ", join(labels, ", "))

    # Compute column widths
    ncols = length(labels)
    col_widths = [max(length(labels[j]), maximum(_numwidth(cm.data[i, j]) for i in 1:ncols)) for j in 1:ncols]
    row_width = maximum(length, labels)

    # Header
    print(io, "  ", lpad("", row_width))
    for j in 1:ncols
        print(io, "  ", lpad(labels[j], col_widths[j]))
    end
    println(io)

    # Rows
    for i in 1:length(labels)
        print(io, "  ", lpad(labels[i], row_width))
        for j in 1:ncols
            print(io, "  ", lpad(_fmtnum(cm.data[i, j]), col_widths[j]))
        end
        if i < length(labels)
            println(io)
        end
    end
end

_fmtnum(x::AbstractFloat) = string(round(x; digits=2))
_fmtnum(x::Integer) = string(x)
_fmtnum(x::Real) = string(round(Float64(x); digits=2))

_numwidth(x) = length(_fmtnum(x))
