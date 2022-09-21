library(gUtils)
library(GxG)

## idk testing this is slightly beyond my pay grade
## GENOME = system.file("extdata", "human_g1k_v37.no.extra.chrom.sizes", package = "GxG")
## if (nchar(GENOME) == 0) {
##     download.file("http://mskilab.com/gUtils/hg19/hg19.chrom.sizes", "~/hg19.chrom.sizes")
##     GENOME="~/hg19.chrom.sizes"
## }
GENOME="human_g1k_v37.no.extra.chrom.sizes"
Sys.setenv(DEFAULT_GENOME = GENOME)

test_that("test .hic", {
    suppressWarnings({
        hic.file = system.file('extdata', "test.hic", package = "GxG")

        ## can load in specific chromosomes via character or numeric vector
        ## and the other normalizations
        gm = straw(hic.file, 1:2, res = 5e4)
        gm = straw(hic.file, c(1:3, 'X'), norm = "NONE", res = 1e5)

        ## we can use any GRanges as input to straw and use alternate norms
        ## like KR and VC, though the default is NONE
        ## (which are not provided with the small .hic matrix bundled with
        ## the package but are available to most Juicer outputs
        gm = straw(hic.file, GRanges('1:1-250e6'), norm = 'NONE', res = 5e5)
        expect_true(inherits(gm, "gMatrix"))

        ## check that gTrack can be created
        gt = gm$gtrack(colormap = c('white', 'green', 'blue'), clim = c(0, 50))
        expect_true(inherits(gt, "gTrack"))

        ## test disjoin
        new.ranges = gr.tile(gm$footprint, 2e5)
        gmd = gm$disjoin(new.ranges)
        expect_true(inherits(gm, "gMatrix"))
        expect_true(length(gmd$gr) >= length(new.ranges))
    })
})

