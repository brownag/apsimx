# suppress messages and warnings based on logical `verbose` argument
# captures expression `expr`, optionally surrounded with suppressWarnings and suppressMessages
.suppressExpression <- function(expr, verbose = FALSE) {
    expr <- substitute(expr)
    .FUN0 <- function(x) x
    .FUN1 <- .FUN0
    .FUN2 <- .FUN0
    if (isFALSE(verbose)) {
        .FUN1 <- suppressMessages
        .FUN2 <- suppressWarnings
    }
    eval(substitute(function() {
        FUN1(FUN2(expr))
    }, list(
        FUN1 = .FUN1,
        FUN2 = .FUN2,
        expr = expr
    )), envir = parent.frame())()
}
