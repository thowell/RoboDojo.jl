function module_dir()
    joinpath(@__DIR__, "../..")
end

using Pkg
Pkg.activate(module_dir())
using RoboDojo
