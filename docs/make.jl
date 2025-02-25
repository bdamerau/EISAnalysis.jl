using Documenter
using EISAnalysis

push!(LOAD_PATH,"../src/")
makedocs(sitename="EISAnalysis.jl",
        doctest = false,
        warnonly = true,
         pages = [
            "Index" => "index.md",
            "Tutorial"=> "Tutorial.md",
            "DRT" => "DRT.md",
            "API" => "api.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/bdamerau/EISAnalysis.jl.git",
    devbranch = "main"
)