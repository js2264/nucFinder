parse <- function(BamFile, ScanBamParam = NULL) {
    if (!is.null(ScanBamParam)) {
        param <- ScanBamParam
    }
    else {
        param <- ScanBamParam(
            flag=scanBamFlag(
                isPaired = TRUE,
                isProperPair = TRUE,
                isDuplicate = FALSE,
                isSecondaryAlignment = FALSE
            ), 
            mapqFilter = 1, 
            what = c("mapq", "isize")
        )
    }
    GenomicAlignments::readGAlignmentPairs(BamFile, param = param) |> 
        as('GRanges')
}
