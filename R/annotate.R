annotate <- function(fft) {

    fft_track <- fft[['fft']]

    ## Map gaps in nucleosome covearge track from FFT track
    pos_signal <- as(fft_track, 'GRanges') |> 
        filter(score >= 0.05) |> 
        reduce() |> 
        filter(width >= 90)
    gaps <- gaps(pos_signal, ignore.strand = TRUE)

    ## Filter NFRs in gaps
    NFRs <- filter(gaps, width >= 70, width <= 1000)

    ## Map dyads from FFT track
    dyads <- lapply(seq_along(fft_track), function(K) {
        chr <- names(fft_track)[K]
        pos <- fft_track[[K]] |> 
            as.vector() |> 
            .rleLocalMaxima()
        gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start = pos, width = 1))
        seqlevels(gr) <- seqlevels(fft_track)
        gr
    }) |> as("GRangesList") |> unlist()

    ## Filter out dyads mapping over gaps
    dyads <- filter_by_non_overlaps(dyads, gaps)

    ## Extend nucleosomes
    nucs <- resize(dyads, fix = 'center', width = 147)

    return(list(
        fft = fft_track, 
        NFRs = NFRs, 
        dyads = dyads, 
        calls = nucs
    ))
}
