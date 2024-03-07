filter_fragments <- function(mnase_frags) {

    ## -- Filter nucleosomal frags
    nuc_frags <- filter(mnase_frags, width <= 200, width >= 130)

    ## -- Get nucleosomal dyads
    nuc_dyads <- resize(nuc_frags, fix = 'center', width = 3)

    ## -- Get FFT track
    nuc_dyad_cov <- GRangesList(nuc_dyads) |> unlist() |> coverage()
    fft <- nucleR::filterFFT(nuc_dyad_cov, pcKeepComp=0.01, mc.cores = 16)
    fft <- map(fft, ~ Rle(.x)) |> RleList()

    return(list(
        dyads = nuc_dyads, 
        fft = fft
    ))
}
