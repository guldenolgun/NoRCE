.onLoad <- function(libname, pkgname) {
  #pkgVersion <- packageDescription(pkg = pkgname, fields="Version")

  msg <- paste0("Welcome to Noncoding RNA Gene Set Enrichment package, ", pkgname,"\n")

  cit <- paste0("If you use ", pkgname, " this tool, please cite:\n")

  packageStartupMessage(paste0(msg, cit))
}
