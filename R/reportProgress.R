#' Report Progress
#'
#' For use within long running scripts. Will report progress in 10% intervals when running iterations.
#' Must start script by defining t1 <- Sys.time()
#'
#' @param n current iteration
#' @param iter total number of iterations to be run
#' @param t1 starting time
#'
#' @return
#' @export
#'
reportProgress <- function(n, iter, t1) {
  if(n == .1 * iter) {
    t.1 <- (Sys.time() - t1)
    cat(t.1, attr(t.1, "units"),"have elapsed: 10% iterations complete","\n")}

  if(n == .2 * iter) {
    t.2 <- (Sys.time() - t1)
    cat(t.2, attr(t.2, "units"),"have elapsed: 20% iterations complete","\n")}

  if(n == .3 * iter) {
    t.3 <- (Sys.time() - t1)
    cat(t.3, attr(t.3, "units"),"have elapsed: 30% iterations complete","\n")}

  if(n == .4 * iter) {
    t.4 <- (Sys.time() - t1)
    cat(t.4, attr(t.4, "units"),"have elapsed: 40% iterations complete","\n")}

  if(n == .5 * iter) {
    t.5 <- (Sys.time() - t1)
    cat(t.5, attr(t.5, "units"), "have elapsed: 50% iterations complete","\n")}

  if(n == .6 * iter) {
    t.6 <- (Sys.time() - t1)
    cat(t.6, attr(t.6, "units"),"have elapsed: 60% iterations complete","\n")}

  if(n == .7 * iter) {
    t.7 <- (Sys.time() - t1)
    cat(t.7, attr(t.7, "units"), "have elapsed: 70% iterations complete","\n")}

  if(n == .8 * iter) {
    t.8 <- (Sys.time() - t1)
    cat(t.8, attr(t.8, "units"), "have elapsed: 80% iterations complete","\n")}

  if(n == .9 * iter) {
    t.9 <- (Sys.time() - t1)
    cat(t.9, attr(t.9, "units"), "have elapsed: 90% iterations complete","\n")}
}
