suppressMessages(library(docopt))

"Detect introgressed segments in a given simulated genome.

Usage:
    find_introgressed_segments.R [options]

Options:
    --fasta=<FILENAME>  input FASTA filename
    --n_afr=<COUNT>     number of simulated Africans
    --n_eur=<COUNT>     number of simulated Europeans
    --n_nea=<COUNT>     number of simulated Neandertals
    --output=<FILENAME> name of the output PDF file
    -h --help           show this help" -> doc

opts <- docopt(doc)

if (with(opts, is.null(fasta) | is.null(n_afr) | is.null(n_eur) | is.null(n_nea) | is.null(output))) {
    writeLines(doc)
    stop("All arguments have to be specified!")
}
opts$n_afr <- as.integer(opts$n_afr)
opts$n_eur <- as.integer(opts$n_eur)
opts$n_nea <- as.integer(opts$n_nea)


library(magrittr)
library(ape)
library(phytools)

# load the simulated sequences and build a NJ tree
sequences <- read.dna(opts$fasta, format = "fasta")
tree <- dist.dna(sequences) %>% nj

# get IDs of tree nodes/leaves belonging to each simulated individual
# (the names of the sequences do not necessarily correspond to their IDs
# within the tree)
afr_taxon_ids <- which(tree$tip.label %in% ((1                       : opts$n_afr)              %>% as.character))
eur_taxon_ids <- which(tree$tip.label %in% ((opts$n_afr + 1          : opts$n_eur)              %>% as.character))
nea_taxon_ids <- which(tree$tip.label %in% ((opts$n_afr + opts$n_eur : opts$n_eur + opts$n_nea) %>% as.character))

# root the tree according to the outgroup (Neandertal), placing the root on
# the branch leading to it
nea_branch_len <- tree$edge.length[ tree$edge[, 2] == nea_taxon_ids ]
rooted_tree <- reroot(tree, nea_taxon_ids, position = 0.2 * nea_branch_len)

pdf(opts$output, 10, 14)

plot(rooted_tree, edge.width = 2, label.offset = 0.0003, main = "Overall phylogeny of simulated sequences")
tiplabels(tip = afr_taxon_ids, bg = "yellow", pch = 24, cex = 1.4)
tiplabels(tip = eur_taxon_ids, bg = "green", pch = 22, cex = 1.7)
tiplabels(tip = nea_taxon_ids, bg = "red", pch = 21, cex = 1.7)
legend("bottomleft", legend = c("AFR", "EUR", "NEA"), pt.bg = c("yellow", "green", "red"),
    pch = c(24, 22, 21))

dev.off()
