# example from http://shiny.rstudio.com/gallery/kmeans-example.html

library(shiny)
library(httr)
library(jsonlite)
library(shinycssloaders)
library(DT)
##library(GenomicRanges)
##library(VariantAnnotation)
##library(BSgenome.Hsapiens.UCSC.hg38)
##library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(g3viz)

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
    v1 <- VariantAnnotation::VRanges(paste0("chr", variants$code),
                  IRanges(as.integer(variants$start), as.integer(variants$end)),
                  ref = variants$referenceAllele,
                  alt = variants$observedAllele)
    
    v1a <- suppressWarnings(
        VariantAnnotation::predictCoding(v1,
                                         TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                                         seqSource=BSgenome.Hsapiens.UCSC.hg38::Hsapiens))
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

VEP <- function(chr, start, end, alt){
    url <- paste0("https://rest.ensembl.org/vep/human/region/",
                  chr, ":", start, ":", end, "/", alt)
    if(Sys.info()[1]=="Darwin"){
        res <- GET(url)
    }else{
        httr_config <- config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
        res <- with_config(config = httr_config, GET(url)) 
    }
    a1 <- fromJSON(httr::content(res, as = "text", encoding = "utf8"))
    t1 <- a1$transcript_consequences[[1]]
    if(all(c("gene_symbol", "consequence_terms", "protein_start", "amino_acids")
           %in% colnames(t1))){
        aa1 <- t1[1, c("gene_symbol", "consequence_terms", "protein_start", "amino_acids")]
        pp <- strsplit(aa1$amino_acid, split = "/")[[1]]
        aa1$amino_acid_change <- paste0("p.", pp[1], aa1$protein_start, pp[2])
    }else{
        aa1 <- NULL
    }
    return(aa1)
}

VarAnnoVEP <- function(variants){
    anno <- c()
    for(i in 1:nrow(variants)){
        print(i)
        x <- variants[i,,drop=F]
        a1 <- VEP(x$code, x$start, x$end, x$observedAllele)
        if(!is.null(a1)){
            anno <- rbind(
                anno,
                cbind(variants[i,c("code", "start", "observedAllele")],
                      a1))
        }
    }
    anno <- anno[,c("code", "start", "observedAllele", "gene_symbol",
                    "consequence_terms", "amino_acid_change", "protein_start")]
    colnames(anno) <- c("Chromosome", "Start_Position", "Tumor_Seq_Allele2",
                        "Hugo_Symbol", "Variant_Classification", "amino_acid_change",
                        "AA_Position")
    return(anno)
}

SearchDGIdb <- function(gene){
    gdb <- GET(paste0("https://dgidb.org/api/v2/interactions.json?genes=", gene))
    gres <- fromJSON(httr::content(gdb, as = "text"))
    intact <- gres$matchedTerms$interactions[[1]]
    return(intact)
}

server <- function(input, output, session){
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
        vars_s <- vars[input$vars_rows_selected,]
        if(nrow(vars_s)>0){
            VarAnnoVEP(vars_s)
        }else{
            VarAnnoVEP(vars)
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
    output$muts <- DT::renderDataTable({mutInput()}, options = list(scrollX = T))
    output$llplot <- renderG3Lollipop({plotInput()})
    output$drug <- DT::renderDataTable({dgiInput()}, options = list(scrollX = T))
}
