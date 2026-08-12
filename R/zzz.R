.onLoad <- function(libname, pkgname) {
  # Register the SAD distribution eagerly when the user calls library(BPFC).
  # Use getFromNamespace so the package name is not hard-coded.
  reg <- utils::getFromNamespace(".mcmcgraph_register_sad", pkgname)
  try(reg(), silent = TRUE)
}
