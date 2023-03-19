.onAttach <- function(libname, pkgname) {
  packageStartupMessage("I see you're loading the Gumbo... nice.
      This is an evolving, often awkwardly implemented
      collection of data analysis tools.
      If you're not in the TYP lab and can't freely ask me
      what on earth I was thinking for a particular function
      you may not want this package. If you are in the lab, I
      hope this makes some of our standard analyses a touch easier.
      -Richard")
}
