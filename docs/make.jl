using Documenter
using LShaped

makedocs(
    sitename = "LShaped",
    format = Documenter.HTML(),
    modules = [LShaped]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
