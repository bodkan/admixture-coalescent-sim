suppressMessages(library(docopt))

"Detect introgressed segments in a given simulated genome.

Usage:
    find_introgressed_segments.R [options]

Options:
    --tree_file=<FILENAME>  file containing genealogies
    --n_afr=<COUNT>     number of simulated Africans
    --n_eur=<COUNT>     number of simulated Europeans
    --n_nea=<COUNT>     number of simulated Neandertals
    --split_time=<TIME>     time of the split between Africans and Eurasians
    --output=<FILENAME> name of the output PDF file
    -h --help               show this help" -> doc

opts <- docopt(doc)

if (with(opts, is.null(tree_file) | is.null(n_afr) | is.null(n_eur) | is.null(n_nea) | is.null(split_time) | is.null(output))) {
    writeLines(doc)
    stop("All arguments have to be specified!")
}
opts$n_afr <- as.integer(opts$n_afr)
opts$n_eur <- as.integer(opts$n_eur)
opts$n_nea <- as.integer(opts$n_nea)
opts$split_time <- as.numeric(opts$split_time)


library(magrittr)
library(ape)
library(phangorn)

# read a list of trees from a given file
trees <- read.tree(opts$tree_file)

pdf(opts$output, 10, 14)

for (i in 1 : length(trees)) {
    tree <- trees[[i]]

    # get IDs of tree nodes/leaves belonging to each simulated individual
    # (the names of the sequences do not necessarily correspond to their IDs
    # within the tree)
    afr_taxon_ids <- which(tree$tip.label %in% ((1                       : opts$n_afr)              %>% as.character))
    eur_taxon_ids <- which(tree$tip.label %in% ((opts$n_afr + 1          : opts$n_eur)              %>% as.character))
    nea_taxon_ids <- which(tree$tip.label %in% ((opts$n_afr + opts$n_eur : opts$n_eur + opts$n_nea) %>% as.character))

    # get index of an immediate parent node of the Neandertal sequence
    parent_node <- Ancestors(tree, node = nea_taxon_ids, type = "parent")
    # get indices of sequences that first coalesce with the Neandertal
    introgressed_seq <-
        Descendants(tree, parent_node, type = "tips") %>%
        unlist %>%
        .[! (. %in% nea_taxon_ids)]
    # get a time of coalescence (i.e. distance within the tree) with Neandertal
    coalesc_time <- dist.nodes(tree)[nea_taxon_ids, parent_node]

    tree_count <- paste(i, "/", length(trees))
    if (coalesc_time < opts$split_time)
        cat("tree ", tree_count, " -- introgression detected!\n")

    plot(tree, edge.width = 2, label.offset = 0.0003, main = paste0(tree_count, " - ", coalesc_time < opts$split_time))
    tiplabels(tip = afr_taxon_ids, bg = "yellow", pch = 24, cex = 1.4)
    tiplabels(tip = eur_taxon_ids, bg = "green", pch = 22, cex = 1.7)
    tiplabels(tip = nea_taxon_ids, bg = "red", pch = 21, cex = 1.7)
    legend("bottomleft", legend = c("AFR", "EUR", "NEA"), pt.bg = c("yellow", "green", "red"),
        pch = c(24, 22, 21), cex = 1.5)
}

dev.off()