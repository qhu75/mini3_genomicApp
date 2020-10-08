# init.R
#
# Example R code to install packages if not already installed
#

## my_packages = c("httr", "jsonlite", "shiny", "shinycssloaders", "DT",
##                 "GenomicRanges", "VariantAnnotation",
##                 "BSgenome.Hsapiens.UCSC.hg38",
##                 "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
##                 "g3viz")

my_packages = c("httr", "jsonlite", "shiny", "shinycssloaders", "DT",
                "org.Hs.eg.db", "g3viz")

install_if_missing = function(p) {
    if(!"BiocManager" %in% rownames(installed.packages())){
        install.packages("BiocManager")
    }
    if (!p %in% rownames(installed.packages())) {
        BiocManager::install(p)
    }
}

invisible(sapply(my_packages, install_if_missing))
