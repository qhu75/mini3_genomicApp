# example from http://shiny.rstudio.com/gallery/kmeans-example.html

library(shiny)

SearchRegion <- function(chr, start, end, gbuild = "hg38"){
    ## RESTful API
    r2 <- GET(paste0("http://hapi.fhir.org/baseR4/MolecularSequence?chromosome=", chr,
                     "&variant-start=ge", start, "&variant-end=le", end))
    j <- i <- 0
    vars <- c()
    while(TRUE){
        j <- j + 1
        i <- i + 20
        r2c <- fromJSON(httr::content(r2, as = "text"))
        if(is.null(r2c$entry)){
            break
        }
        v1 <- r2c$entry$resource$variant
        v1 <- lapply(v1, function(x)x[,c("start", "end", "referenceAllele", "observedAllele")])
        patient <- r2c$entry$resource$patient$reference
        purl <- paste0("http://hapi.fhir.org/baseR4/", patient)
        pt <- paste0("<a href='",purl,"' target='_blank'>",patient,"</a>")
        gbuild <- r2c$entry$resource$referenceSeq$genomeBuild
        chr1 <- r2c$entry$resource$referenceSeq$chromosome$coding
        vars <- cbind(do.call(rbind, chr1), do.call(rbind, v1))
        vars <- data.frame(pt, gbuild, vars)
        link2 <- paste0("http://hapi.fhir.org/baseR4?_getpages=", r2c$id,
                        "&_getpagesoffset=",i,
                        "&_count=20&_pretty=true&_bundletype=searchset")
        r2 <- GET(link2)
    }
    ## vars <- vars[vars$start > start & vars$end <= end,]
    if(nrow(vars) > 0){
        vars <- vars[vars$gbuild == gbuild & vars$code == chr &
                     vars$start >= start & vars$end <= end,,drop=F]
        if("system" %in% colnames(vars)){
            vars <- vars[,-match("system", colnames(vars))]
        }
    }
    rownames(vars) <- NULL
    return(vars)
}

gene2reg <- function(gene){
    gid <- get(gene, org.Hs.egSYMBOL2EG)
    start <- get(gid, org.Hs.egCHRLOC)[1]
    end <- tail(get(gid, org.Hs.egCHRLOCEND), 1)
    c(start, end)
}

VarAnno <- function(variants){
    v1 <- VRanges(paste0("chr", variants$code),
                  IRanges(as.integer(variants$start), as.integer(variants$end)),
                  ref = variants$referenceAllele,
                  alt = variants$observedAllele)
    
    v1a <- suppressWarnings(predictCoding(v1, TxDb.Hsapiens.UCSC.hg38.knownGene, seqSource=Hsapiens))
    v1a <- v1a[!is.na(v1a$GENEID),]
    idx <- findOverlaps(v1, v1a, type = "equal")
    idx <- na.omit(subjectHits(idx)[match(seq(v1), queryHits(idx))])
    v1a <- v1a[idx,]
    
    ## to MAF
    ploc <- unlist(lapply(v1a$PROTEINLOC, "[[", 1))
    muts <- data.frame(Chromosome = as.character(seqnames(v1a)),
                       Start_Position = start(v1a),
                       Tumor_Seq_Allele2 = v1a$varAllele,
                       Hugo_Symbol = unlist(mget(v1a$GENEID, org.Hs.egSYMBOL)),
                       Variant_Classification = v1a$CONSEQUENCE,
                       amino_acid_change = paste0("p.", v1a$REFAA, ploc, v1a$VARAA),
                       AA_Position = ploc)
    return(muts)
}

SearchDGIdb <- function(gene){
    gdb <- GET(paste0("https://dgidb.org/api/v2/interactions.json?genes=", gene))
    gres <- fromJSON(httr::content(gdb, as = "text"))
    intact <- gres$matchedTerms$interactions[[1]]
    return(intact)
}

server <- function(input, output){
    datasetInput <- eventReactive(input$search, {
        gene <- input$gene
        region <- input$region

        if(gene != ""){
            regs <- gene2reg(gene)
            SearchRegion(names(regs)[1], regs[1], regs[2])
        }else if(region != ""){
            regs <- strsplit(region, split = "[:-]")[[1]]
            SearchRegion(regs[1], as.integer(regs[2]), as.integer(regs[3]))
        }
    })

    mutInput <- eventReactive(input$annot, {
        vars <- datasetInput()
        if(!is.null(vars)){
            VarAnno(vars)
        }
    })

    plotInput <- eventReactive(input$lplot, {
        muts <- mutInput()
        g3Lollipop(muts, gene.symbol = muts$Hugo_Symbol[1],
                   protein.change.col = "amino_acid_change",
                   factor.col = "Variant_Classification")
    })

    dgiInput <- eventReactive(input$dgi, {
        muts <- mutInput()
        SearchDGIdb(muts$Hugo_Symbol[1])
    })
    
    output$vars <- DT::renderDataTable({datasetInput()}, escape = FALSE)
    output$muts <- DT::renderDataTable({mutInput()})
    output$llplot <- renderG3Lollipop({plotInput()})
    output$drug <- DT::renderDataTable({dgiInput()}, options = list(scrollX = T))
}
