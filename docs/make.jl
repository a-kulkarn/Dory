push!(LOAD_PATH,"../src/")
using Documenter, DocumenterMarkdown, Dory, Pkg
makedocs(sitename="Dory.jl Documentation")


# makedocs(
#     doctest= false,
#     modules = [Dory],
#     format = Markdown(),
# )

# docsdir = joinpath(@__DIR__, "build/")

# for (root, dirs, files) in walkdir(docsdir)
#   for file in files
#     filename = joinpath(root, file) # path to files
#     run(`sed -i 's/.*dash; \*Method.*/---/g' $filename`)
#     run(`sed -i 's/.*dash; \*Type.*/---/g' $filename`)
#     run(`sed -i 's/.*dash; \*Function.*/---/g' $filename`)
#     run(`sed -i '/>source<\/a>/d' $filename`)
#     run(`sed -i '/>\#<\/a>/d' $filename`)
#   end
# end

deploydocs(
  repo = "github.com/a-kulkarn/Dory.jl.git",
  deps = Deps.pip("pymdown-extensions", "pygments", "mkdocs", "python-markdown-math", "mkdocs-material", "mkdocs-cinder"),
  target = "site",
  make = () -> run(`mkdocs build`),
)

#makedocs(
#    modules = [Hecke, Nemo, AbstractAlgebra],
#    clean   = true,
#    format = :html,
#    sitename = "Hecke",
#    doctest = !false,
#    pages = [
#      "index.md",
#      "number_fields/intro.md",
#      "Orders" => [ "orders/introduction.md",
#                    "orders/orders.md",
#                    "orders/elements.md",
#                    "orders/ideals.md",
#                    "orders/frac_ideals.md"
#                  ],
##      "Maximal Orders" => [ "MaximalOrders/Introduction.md",
##                            "MaximalOrders/Creation.md",
##                            "MaximalOrders/Elements.md",
##                            "MaximalOrders/Ideals.md"
##                          ],
#      "abelian/introduction.md",
#      "class_fields/intro.md",
#      "sparse/intro.md",
#      "FacElem.md"
#      ]
#)
#
## Hack around to get syntax highlighting working
##cd(joinpath(dirname(pathof(Hecke)), "..", "docs"))
##
##cp("application-f78e5cb881.palette.css", "build/application-f78e5cb881.palette.css", force = true)
##cp("application-e2807e330f.css", "build/application-e2807e330f.css", force = true)
#
#deploydocs(
#     deps = nothing,
##    deps = Deps.pip()#="pygments",
##                    "mkdocs==0.16.3",
##                    "python-markdown-math",
##                    "mkdocs-cinder")=#,
#    repo = "github.com/thofma/Hecke.jl.git",
#    target = "build",
#    make = nothing,
#    osname = "linux",
#    julia = "1.0",
#)
#
## Try out the following deps = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4" ,"python-markdown-math", "pygments", "pymdown-extensions")
