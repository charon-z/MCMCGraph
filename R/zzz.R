.onLoad <- function(libname, pkgname) {
  # ✅ 双保险之二：用户 library(MCMCGraph) 时就提前注册
  # 用 getFromNamespace 避免写死包名
  reg <- utils::getFromNamespace(".mcmcgraph_register_sad", pkgname)
  try(reg(), silent = TRUE)
}
