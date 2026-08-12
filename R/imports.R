#' @keywords internal
"_PACKAGE"

## Namespace directives ------------------------------------------------------

# nimble's model DSL exposes helpers (inprod, nimNumeric, nimPrint, the in-model
# `rnorm`, etc.) that R CMD check cannot resolve statically; importing the whole
# nimble namespace makes them visible and silences the false positives.
#' @import nimble
#' @importFrom graphics par plot
#' @importFrom stats density rnorm var kmeans
NULL

# `Z0` is a node name that only exists inside the nimbleCode() model block, not
# an R object, so declare it as a known global to avoid a spurious NOTE.
utils::globalVariables(c("Z0"))
