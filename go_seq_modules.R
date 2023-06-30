library(dplyr)
library(splitstackshape)
library(mgsub)
library(httr)
library(jsonlite)
library(xml2)
library(purrr)
library(forcats)
library(ggplot2)
library(svglite)
library(goseq)
library(data.table)
# Calculate genelength for PVIT genes
data <- read.table("path/to/plasmopara_viticola_INRA-PV221.geneAnnotation.v1.gff3", h=F)
test <- data %>% filter(V3 %in% "gene") %>% mutate(Length= V5-V4, Gene= mgsub(V9,c("^ID=",";Name=.*"), c("","")),Transcripts= mgsub(V9,c("^ID=",";Name=.*"), c("","")) ) %>% select(Gene,Transcripts,Length)
genelength <- read.table(file = "path/to/Genome_Vitis/GO/gene_length/full_gene_length.txt", sep="\t", header = TRUE)
genelength <- rbind(genelength, test)
# Retrieve GO terms for PVIT
go1 <- fread("path/to/plasmopara_viticola_INRA-PV221.functionalAnnot.blast2go.tsv", select= c("SeqName", "GO IDs"))
names(go1) <- c("Gene", "GO")
go1$Gene <- gsub(".T\\d", "", go1$Gene)
go1$GO <- mgsub(go1$GO, c("C:","F:", "P:"), c("","", ""))
go2 <-  as.data.frame(cSplit(go1, 2, ";", "long"))
GO_terms = read.table(file = "path/to/Genome_Vitis/GO/full_GO_clean.txt", sep="\t", header = TRUE)
GO_terms = rbind(GO_terms, go2)
# Get background info (genes that pass threshold for Deseq2 matrix)
dat <- read.csv("path/to/Deseq2/CI_vs_C_24h-Deseq2_matrix.txt")
background = rownames(dat)
View(dat)
background <- dat[,1]
genelength_lim <- genelength %>% filter(Gene %in% background)
GO_terms  <- GO_terms %>% filter(Gene %in% background)
# Retrieve genes in WGCNA calculated modules
mod <- read.csv("path/to/C_geneInfo.csv")
#MM <- lapply(split(mod, mod$clust), function(x) unname(unlist(x$gene)))
### Module of interest
cc <- c("black", "green", "magenta")
MM_6 <- list()
MM_8 <- list()
MM <- list()
genes = list()
for ( i in 1:length(cc))
{
MM_8[[i]] <- unname(unlist(mod %>% filter(abs(get(paste0("MM.", cc[i]))) > 0.8) %>% select(Gene)))
MM_6[[i]] <- unname(unlist(mod %>% filter(abs(get(paste0("MM.", cc[i]))) > 0.6) %>% select(Gene)))
MM[[i]] <- unname(unlist(mod %>% filter(moduleColor %in% cc[i]) %>% select(Gene)))
}

genes <- lapply(MM, function(x) as.integer(background %in% x))
genes <- lapply(genes, setNames, background)
# Vector genes length
genl <- as.integer(genelength_lim[,3])
#Analysis
pwf <- lapply(genes, function(x) nullp(x, bias.data = genl))
go.wall <- lapply(pwf, function(x) goseq(x, gene2cat = GO_terms))
#Plotting
#Plot MF
go.wall <- lapply(go.wall, function(x) 
				{rownames(x) <- 1:length(rownames(x))
				x})
reduced_go <- lapply(go.wall, function(x) x %>% filter(over_represented_pvalue < 1))
line_NA <- lapply(reduced_go, function(x) which(is.na(x$term)))
miss_GO <- lapply(seq_along(reduced_go), function(x) reduced_go[[x]][line_NA[[x]], -c(6,7)])
reduced_go <- lapply(seq_along(reduced_go), function(x) reduced_go[[x]][-c(line_NA[[x]]),])
#lapply(seq_along(miss_GO), function(x) writeLines(as.vector(miss_GO[[x]]$category), paste0("miss_GO_", cc[x], ".txt"), sep = ","))
#Retrieve missing GO
# Send request to QuickGo API
test = lapply(miss_GO, function(x) na.omit(gsub("GO:", "", x$category)))
requestURL <- lapply(test, function(y) lapply(y, function(x) paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=GO%3A", x, "&limit=5&page=1")))
miss <- lapply(requestURL, function(y) lapply(y, function(x) GET(x, accept("application/json"))))
miss <- lapply(miss, function(y) lapply(y, stop_for_status))
miss <- lapply(miss, function(y) lapply(y, function(x) fromJSON(toJSON(content(x)))))
miss <-lapply(miss, function(y) lapply(y, function(x) as.data.frame(sapply(x$results, unlist))))
# FromJSON retrieve corresponding GO
# Need to fix bug with filter and list
fGO <-lapply(seq_along(miss), function(x) lapply(seq_along(miss[[x]]), function(y) miss[[x]][[y]] %>% dplyr::filter(id %in% miss_GO[[x]]$category[y]) %>% dplyr::select(id, name, aspect)%>% dplyr::rename("category" = "id" , "term" = "name", "ontology" = "aspect")))
# Create dataframe containing missing GOs for each module
missed <- lapply(fGO, function(x) Reduce(rbind,x))
# Formatting correctly the ontology column to match with ontology column from reduced_go
missed <- missed %>% map(~ .x %>% mutate(ontology = mgsub(ontology, c("molecular_function", "biological_process", "cellular_component"), c("MF", "BP", "CC"))))
missed <- lapply(missed , function(x) lapply(x, function(y) unlist(y)))
# Gather miss_GO and retrieve missing GOs
miss_GO <- lapply(seq_along(miss_GO), function(x) merge(miss_GO[[x]], missed[[x]], by = "category")) 
# Merge together with reduced_go
merged <- lapply(seq_along(reduced_go), function(x) rbind(reduced_go[[x]], miss_GO[[x]]))
merged <- lapply(merged, function(x )x[order(x$ontology),])
merged <- merged %>% map( ~.x %>% mutate(hitsPerc=numDEInCat*100/numInCat))
# Create top 20 GOs terms by ontology
final <- merged %>% map( ~.x %>% group_by(ontology) %>%  filter(over_represented_pvalue < 0.2) %>% top_n(20, wt=-over_represented_pvalue))
# If need to reduce term name do it here
# Check if term name too long
too_long <- lapply(final, function(x) which(nchar(x$term) > 80))
#Manual modify name example : final[[3]]$term[too_long[[3]]] <- "..."
final <- final %>% map( ~.x%>% mutate(term = factor(term, levels= term), category = factor(category, levels= category))) 
# Coordinate for dash line in plot (separate ontology)
loc_pro <- lapply(final, function(x) as.vector(table(x$ontology)))
### Required to built the second y-axis
# solution from https://github.com/tidyverse/ggplot2/issues/3171

guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}
#usage :
#guides(y.sec = guide_axis_label_trans(~labels_desired_for_y_axis))

### Loop not working for the moment in Rstudio
for (i in seq_along(final))
{
name_ana = paste0("CH_", cc[[i]])
print(ggplot(final[[i]], aes(x=numDEInCat, y=fct_rev(term), colour=over_represented_pvalue, size=hitsPerc)) +
    geom_point(alpha = 0.7) +
    expand_limits(x=0) +
    labs(x="Number of genes", y="GO term", colour="p value", size="Enrichment (%)") +
    scale_colour_gradient2(midpoint = 0.1, mid= "grey", high = "lightgrey", low = "red") +
    scale_size(range = c(1,15)) +
    theme_bw() + theme(legend.position="bottom", legend.box = "vertical") +
    guides(y.sec = guide_axis_label_trans(~rev(final[[i]]$category))) +
    geom_hline(yintercept = length(final[[i]]$term) - (loc_pro[[i]][1] -0.5), linetype = "dotted") +
    geom_hline(yintercept = length(final[[i]]$term) - (loc_pro[[i]][1] + loc_pro[[i]][2] -0.5), linetype = "dotted"))
    ggsave(paste0(name_ana, "_top20_GO_MM.svg"), height = 12, width = 8)

}
