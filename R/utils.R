.rleLocalMaxima <- function(x) {
    # see https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
        y <- y[-1]
    }
    y
}
