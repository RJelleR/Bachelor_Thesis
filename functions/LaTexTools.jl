module LaTexTools

using DataFrames, PrettyTables
export save_latex_table

function save_latex_table(df::DataFrame, title::String, file_path::String; transpose=false, row_index=true)
    base_dir = normpath(joinpath(@__DIR__, ".."))
    latex_dir = joinpath(base_dir, "latex_project")

    if !isdir(latex_dir)
        error("LaTeX project directory not found at: $latex_dir")
    end

    path = joinpath(latex_dir, file_path * ".tex")

    if transpose
        df = DataFrame(permutedims(Matrix(df), ["col_" * string(i) for i in 1:size(df, 1)]), :auto)
    end

    open(path, "w") do io
        pretty_table(io, df;
            backend = Val(:latex),
            title = title,
            alignment = :c,
            header = names(df),
            tf = tf_latex_booktabs
        )
    end

    return path
end
end 