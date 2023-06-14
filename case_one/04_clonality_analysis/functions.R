facets_run = function(WD, cval1) {
    library("facets")
    library("future")
    plan(multicore)
    set.seed(666)
    setwd(WD)
    out = list()
    for (i in seq(1, length(list.files("pileup")))) {
        out[[i]] <- future({
            out_frame = data.frame(matrix(nrow = 1, ncol = 5))
            j = 1
            ea = list.files("pileup")[i]
            snp <- readSnpMatrix(paste0("pileup/", ea))
            xx <- preProcSample(snp, ndepth = 10, het.thresh = 0.1, ndepthmax = 10000, cval = cval1)
            pdf(paste0("facets/", gsub("\\..*", "", ea), ".pdf"))
            for (cval in seq(floor((cval1 / 100) + 1) * 100, 500, 50)) {
                name = paste(gsub("\\..*", "", ea), "cval", cval, sep = ".")
                oo <- procSample(xx, cval = cval, dipLogR = 0)
                fit <- emcncf(oo)
                out_frame[j, ] = data.frame(cval = cval, loglike = fit$loglik, purity = fit$purity, ploidy = fit$ploidy, dipLogR = fit$dipLogR)
                j = j + 1
                plotSample(oo, fit, sname = name, plot.type = "em")
                fit$cncf$`Major Copy Number` <- fit$cncf$tcn.em - fit$cncf$lcn.em
                colnames(fit$cncf)[is.element(colnames(fit$cncf), c("cf.em", "tcn.em", "lcn.em"))] <- c("Cellular Fraction", "Total Copy Number", "Minor Copy Number")
                fwrite(fit$cncf, file = paste0("facets/", name, ".csv"))
            }
            dev.off()
            colnames(out_frame) = c("cval", "logLike", "purity", "ploidy", "dipLogR")
            fwrite(out_frame, file = paste0("facets/", gsub("\\..*", "", ea), ".stats.csv"))
        })
    }
}
cnAssign = function(cnv, name) {
    out = data.frame(matrix(nrow = 1, ncol = 3))
    colnames(out) = c("Cellular.Fraction", "Major.Copy.Number", "Minor.Copy.Number")
    for (i in seq(1, nrow(gtf))) {
        temp = (cnv[intersect(intersect(which(cnv$chrom == gtf$chr[i]), which(cnv$start < gtf$mean[i])), which(cnv$end > gtf$mean[i])), ])
        if (length(temp$Cellular.Fraction) > 0) {
            out = rbind(out, temp[, c(12, 15, 14)])
        } else {
            out = rbind(out, rep(NA, 3))
        }
    }
    out = out[-1, ]
    colnames(out) = paste(name, colnames(out), sep = ".")
    gtf = cbind(gtf, out)
    gtf
}
citupToTree = function(citDir, snpDir, outDir) {
    library("treemap")
    library("data.tree")
    library("ape")
    library("ggtree")
    library("DiagrammeRsvg")
    library("DiagrammeR")
    setwd(citDir)
    objValue = c()
    for (ea in list.files("tmp/tmp/tree/")) {
        setwd(citDir)
        setwd(paste0("tmp/tmp/tree/", ea))
        if (length(grep("results", list.files())) == 0) {
            next
        }
        print(ea)
        objValue = c(objValue, as.numeric(readLines("results", n = 2)[2]))
    }

    setwd(snpDir)
    loci = read.csv("loci.tsv", sep = "\t", header = T)
    loci$mutation_id = gsub("C10orf54", "VSIR", loci$mutation_id)
    loci$mutation_id = gsub("C14orf37", "ARMH4", loci$mutation_id)
    loci$mutation_id = gsub("MLL4", "KMT2B", loci$mutation_id)
    cluster = read.csv("cluster.tsv", sep = "\t", header = T)
    x = aggregate(loci$variant_allele_frequency, by = list(loci$cluster_id), FUN = mean)
    x$ord = rank(x$x)
    x$ord = rank(1 - x$x)
    loci$cluster_id = x$ord[match(loci$cluster_id, x$Group.1)]
    setwd(citDir)
    l = 1
    for (ea in list.files("tmp/tmp/tree/")) {
        setwd(citDir)
        setwd(paste0("tmp/tmp/tree/", ea))
        if (length(grep("results", list.files())) == 0) {
            next
        }
        # 	if(as.numeric(readLines('results', n=2)[2])>33.55){next}
        x = readLines("results", n = 15)[-seq(1, 3)]
        x = x[-seq(grep("gamma_matrix", x), length(x))]
        x = as.data.frame(t(sapply(strsplit(x, " "), unlist)))
        pcTable = data.frame(matrix(nrow = length(unique(unlist(x))), ncol = 3))
        i = 1
        for (ea in unique(unlist(x))) {
            parent = paste(x[which(x[, 2] == ea), 1], collapse = ",")
            child = paste(x[which(x[, 1] == ea), 2], collapse = ",")
            if (length(parent) == 0) {
                parent = NA
            }
            if (length(child) == 0) {
                child = NA
            }
            pcTable[i, ] = c(ea, parent, child)
            i = i + 1
        }
        colnames(pcTable) = c("clone", "parent", "child")
        parent = pcTable$clone[which(nchar(pcTable$parent) == 0)]
        pPrime = 99
        outF = c()
        for (tail in pcTable$clone[which(nchar(pcTable$child) == 0)]) {
            pPrime = i = 99
            out = c()
            while (pPrime != parent) {
                pPrime = pcTable$parent[which(pcTable$clone == tail)]
                if (i == 99) {
                    out = c(pPrime, tail, out)
                } else {
                    out = c(pPrime, out)
                }
                tail = i = pPrime
            }
            outF = c(outF, paste(out, collapse = "/"))
            print(outF)
        }
        outF = data.frame(pathString = outF)
        if (nrow(outF) == 1) {
            next
        }
        outNode = as.Node(outF)
        labelVec = data.frame(matrix(nrow = 7, ncol = 1))
        for (ea in as.numeric(unlist(strsplit(as.character(paste(unlist(outF), collapse = "/")), "/")))) {
            if (ea == 0) {
                next
            }
            ea = as.character(ea)
            temp = loci[which(loci$cluster_id == ea), ]
            temp = temp[which(temp$variant_allele_frequency > 0), ]
            temp2 = aggregate(temp[, c(1, 3)], by = list(temp$mutation_id), FUN = mean)
            temp2 = temp2[order(temp2$variant_allele_frequency, decreasing = T), ]
            if (length(grep("del|amp", temp2$Group.1)) > 0) {
                temp2 = temp2[c(grep("del|amp", temp2$Group.1), seq(1, nrow(temp2))[-grep("del|amp", temp2$Group.1)]), ]
            }
            labelVec[[ea]] = temp2$Group.1[1:7]
        }
        labelVecL = lapply(labelVec[, -1], paste, collapse = "\n")
        labelVecL = lapply(labelVecL, function(x) gsub("\nNA", "", x))
        outNode$Do(function(node) SetEdgeStyle(node, inherit = FALSE, arrowhead = "vee", color = "grey35", penwidth = 2, label = labelVecL[as.numeric(node$name)]))
        setwd(outDir)
        treeAsSVG <- export_svg(render_graph(ToDiagrammeRGraph(outNode)))
        writeLines(treeAsSVG, paste0("file", l, ".svg"))
        l = l + 1
        asdf
    }
    setwd(outDir)
    system("for ea in file*svg; do sem -j 15 inkscape $ea --export-pdf=$ea.pdf; done")
    system("qpdf --empty --pages *.pdf -- out.pdf")
}
