using Documenter
using EISAnalysis

makedocs(
    sitename = "EISAnalysis",
    format = Documenter.HTML(),
    modules = [EISAnalysis]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
