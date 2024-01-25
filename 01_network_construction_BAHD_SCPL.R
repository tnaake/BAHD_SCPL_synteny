## setwd
setwd("~/AG-Fernie/Thomas/Data/synteny")

library(igraph)

## get TPS genes from output of Orthofinder or MCL
genes_table <- read.table("./Results_Nov05/family_orthogroups.txt", sep = "\t", 
    header = FALSE, stringsAsFactors = FALSE)
bs_genes <- sort(genes_table[genes_table[, 2] %in% ## for OrthoFinder 
    c("OG0000346", "OG0000212", "OG0002767", "OG0000119", "OG0001133", ## BAHD-ATs
      "OG0001868", "OG0000161", "OG0000365", "OG0001959", "OG0002199", ## BAHD-ATs
      "OG0000185", "OG0000121", "OG0001568", "OG0001444", "OG0003193", ## SCP&SCPL-ATs
      "OG0003286"), 1]) ## SCP&SCPL-ATs                 
  
genes_table_mcl <- read.table("./Results_Nov05/family_orthogroups_mcl15.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
bs_genes_mcl <- sort(genes_table_mcl[genes_table_mcl[, 2] %in% ## for mcl
    c("group_414", "group_298", "group_1084", "group_179", "group_56", ## BAHD-ATs
      "group_50"), 1]) ## SCP-SCPL-ATs

## vector with all TPS genes identified by orthofinder and mcl, use this to create the bin_mat_... matrices
bs_genes_all <- sort(unique(c(bs_genes, bs_genes_mcl)))

## make the genes unique
bs_genes <- unique(bs_genes)
bs_genes_mcl <- unique(bs_genes_mcl)
bs_genes_all <- unique(bs_genes_all)

# ## load information on pks genes
# supp <- openxlsx::read.xlsx("~/winhome/Documents/03_Syntenic_linkage/01_Data/synteny_network_results/pks_genes_tree_id_type.xlsx", 
#             sheet = "pks_genes_tree_id_type")
# supp <- supp[, 1:27]

## save the pks_genes
setwd("~/AG-Fernie/Federico/BAHD_SCPL_synteny/")
save(bs_genes, bs_genes_mcl, bs_genes_all, file = "./bs_genes.RData")

## i-ADHore
## function to find synteny by finding pks_genes in i-ADHore output
find_synteny <- function(files, bin_mat, genes) {
    ## loop through files
    for (i in seq_along(files)) {
        files_i <- paste0(files[i], "/anchorpoints.txt")
        mult_pairs <- tryCatch(read.csv(files_i, header = FALSE, sep = "\t", 
            skip = 1, fill = TRUE, comment.char = "#", 
            stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(mult_pairs)) {
            gene_1_ind <- which(mult_pairs[, 4] %in% genes)
            gene_2_ind <- which(mult_pairs[, 5] %in% genes)
            gene_ind <- unique(c(gene_1_ind, gene_2_ind))
            
            links <- mult_pairs[gene_ind, 4:5]
            links <- links[!(!(links[, 1] %in% genes) | !(links[, 2] %in% genes)), ]
            
            ## write 1 to collinear pair when there is synteny reported
            for (j in seq_len(nrow(links))) {
                bin_mat[links[j, 1], links[j, 2]] <- bin_mat[links[j, 2], links[j, 1]] <- 1
            }
            
        } else {print(files_i)}
    }
    return(bin_mat)
}

## function to find tandem by finding pks_genes in i-ADHoRe output
find_tandem <- function(files, bin_mat, genes) {
    
    tandem_list <- vector("list", length(files))
    ## loop through files
    for (i in seq_along(files)) {
        tandem <- read.csv(paste(files[i], "genes.txt", sep = "/"), 
            header = TRUE, sep = "\t", fill = TRUE, comment.char = "#", 
            stringsAsFactors = FALSE)
        tandem <- tandem[tandem[, "id"] %in% genes, c("id", "tandem_representative")]
        tandem_ind <- which(tandem[, "tandem_representative"] != "")
        tandem <- tandem[tandem_ind,]
        
        tandem_list[[i]] <- tandem
        ## write 1 to tandem when there is a tandem gene reported
        for (j in seq_len(nrow(tandem))) {
            bin_mat[tandem[j, 1], tandem[j, 2]] <- bin_mat[tandem[j, 2], tandem[j, 1]] <- 1
        }
    }
    return(list(bin_mat, tandem_list))
}

## load i-ADHore results from Orthofinder input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore/output/")
output_files <- list.files()[grep(list.files(), pattern = "output_")]

## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem <- bin_mat <- matrix(data = 0, nrow = length(bs_genes_all), 
    ncol = length(bs_genes_all))
colnames(bin_mat) <- rownames(bin_mat) <- bs_genes_all
colnames(bin_mat_tandem) <- rownames(bin_mat_tandem) <- bs_genes_all
bin_mat <- find_synteny(output_files, bin_mat, bs_genes) ## use here bs_genes 
tandem <- find_tandem(output_files, bin_mat_tandem, bs_genes)

## load i-ADHore results from MCL input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore/output_mcl")
output_files <- list.files()[grep(list.files(), pattern = "output_")]
## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem_mcl <- bin_mat_mcl <- matrix(data = 0, nrow = length(bs_genes_all), 
    ncol = length(bs_genes_all))
colnames(bin_mat_mcl) <- rownames(bin_mat_mcl) <- bs_genes_all
colnames(bin_mat_tandem_mcl) <- rownames(bin_mat_tandem_mcl) <- bs_genes_all
bin_mat_mcl <- find_synteny(output_files, bin_mat_mcl, bs_genes_mcl) ## use here bs_genes_mcl
tandem_mcl <- find_tandem(output_files, bin_mat_tandem_mcl, bs_genes_mcl)

## MCScanX
## function to find synteny by finding pks_genes in MSScanX output
find_synteny_mcscanx <- function(files, bin_mat, genes) {
    ## loop through files
    for (i in 1:length(files)) {
        collinear_i <- NULL
        collinear_i <- tryCatch(read.csv(files[i], sep = "\t", fill = TRUE, 
            stringsAsFactors = FALSE, header = FALSE, comment.char = "#"), 
            error = function(x) NULL)
        if (!is.null(collinear_i)) {
            ## get index
            gene_1_ind <- which(collinear_i[, 2] %in% genes)
            gene_2_ind <- which(collinear_i[, 3] %in% genes)
            gene_ind <- unique(c(gene_1_ind, gene_2_ind))
            
            links <- collinear_i[gene_ind, 2:3]
            links <- links[!(!(links[, 1] %in% genes) | !(links[, 2] %in% genes)), ]
            
            ## write 1 to collinear pair when there is synteny reported
            for (j in seq_len(nrow(links))) {
                bin_mat[links[j,1], links[j, 2]] <- bin_mat[links[j, 2], links[j, 1]] <- 1
            }    
        } else (print(files[i]))
    }
    return(bin_mat)
}

## function to find tandems by finding pks_genes in MCScanX output
find_tandem_mcscanx <- function(files, bin_mat, genes) {
    ## loop through files
    tandem_list <- list()
    for (i in 1:length(files)) {
        file_i <- paste0(strsplit(files[i], split = "[.]collinearity")[[1]], ".tandem")
        tandem_i <- tryCatch(read.csv(file_i, sep = ",", 
            stringsAsFactors = FALSE, header = FALSE, comment.char = "#"), 
            error = function(x) NULL)
        if (!is.null(tandem_i)) {
            gene_1_ind <- which(tandem_i[,1] %in% genes)
            gene_2_ind <- which(tandem_i[,2] %in% genes)
            gene_ind <- unique(c(gene_1_ind, gene_2_ind)) 
            tandem <- tandem_i[gene_ind, ]
            
            tandem_list[[i]] <- tandem    
            ## write 1 to tandem pair when there is synteny reported
            for (j in seq_len(nrow(tandem))) {
                bin_mat[tandem[j, 1], tandem[j, 2]] <- bin_mat[tandem[j, 2], tandem[j, 1]] <- 1
            }
        } else (print(file_i))
    }
    return(list(bin_mat, tandem_list))
} ## replaces .collineary by .tandem

## load MCScanX results from Orthofinder input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/MCScanX_gff")
collinearity_files <- list.files()[grep(list.files(), pattern = "[.]collinearity")]

## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem_mcscanx <- bin_mat_mcscanx <- matrix(0, 
    nrow = length(bs_genes_all), ncol = length(bs_genes_all))
colnames(bin_mat_mcscanx) <- rownames(bin_mat_mcscanx) <- bs_genes_all
colnames(bin_mat_tandem_mcscanx) <- rownames(bin_mat_tandem_mcscanx) <- bs_genes_all
bin_mat_mcscanx <- find_synteny_mcscanx(collinearity_files, bin_mat_mcscanx, bs_genes)
tandem_mcscanx <- find_tandem_mcscanx(collinearity_files, bin_mat_tandem_mcscanx, bs_genes)


## load MScanX results from MCL input
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_MCScanX/MCScanX_gff_mcl")
collinearity_files <- list.files()[grep(list.files(), pattern = "[.]collinearity")]

## create binary matrix that will store if there exists a collinear or tandem
## relationship or not
bin_mat_tandem_mcscanx_mcl <- bin_mat_mcscanx_mcl <- matrix(0, 
    nrow = length(bs_genes_all), ncol = length(bs_genes_all))
colnames(bin_mat_mcscanx_mcl) <- rownames(bin_mat_mcscanx_mcl) <- bs_genes_all
colnames(bin_mat_tandem_mcscanx_mcl) <- rownames(bin_mat_tandem_mcscanx_mcl) <- bs_genes_all
bin_mat_mcscanx_mcl  <- find_synteny_mcscanx(collinearity_files, bin_mat_mcscanx_mcl, bs_genes_mcl)
tandem_mcscanx_mcl <- find_tandem_mcscanx(collinearity_files, bin_mat_tandem_mcscanx_mcl, bs_genes_mcl)

## save and load tandem and bin_mat files
setwd("~/AG-Fernie/Federico/BAHD_SCPL_synteny/")
save(bin_mat, bin_mat_mcl, bin_mat_mcscanx, bin_mat_mcscanx_mcl, 
    file = "./synteny_bin_mat.RData")
save(tandem, tandem_mcl, tandem_mcscanx, tandem_mcscanx_mcl, 
    file = "./synteny_tandem_mat.RData")
load("./synteny_bin_mat.RData")
load("./synteny_tandem_mat.RData")

##### TANDEM #####
## how to proceed with tandems? get tandems for all species from the four 
## methods and use the lowest level union of tandem affiliation to 
## concatenate/paste protein names, i.e. for
## method 1: x1 and x2 and x4 are a tandem
## method 2: x1 and x2 and x3 are a tandem 
## then x1 and x2 and x3 and x4 are a tandem
all(rownames(tandem[[1]]) == rownames(tandem_mcl[[1]]))
all(rownames(tandem[[1]]) == rownames(tandem_mcscanx[[1]]))
all(rownames(tandem[[1]]) == rownames(tandem_mcscanx_mcl[[1]]))
tandem_sum <- tandem[[1]] + tandem_mcl[[1]] + tandem_mcscanx[[1]] + tandem_mcscanx_mcl[[1]]
rownames(tandem_sum) <- rownames(tandem[[1]])
ind_keep <- apply(tandem_sum, 1, sum) > 0

tandem_sum <- tandem_sum[ind_keep, ind_keep]
##tandem_sum[which(tandem_sum>0)] <- 1
diag(tandem_sum) <- 0

## create a matrix that stores the source of tandem, i.e. if the information
## comes from i-ADHore + Orthofinder, i-ADHore + MCL, MCScanX + Orthofinder 
## and/or MCScanX + MCL
tandem_type <- matrix("", ncol = ncol(tandem[[1]]), nrow = nrow(tandem[[1]]))
rownames(tandem_type) <- rownames(tandem[[1]])
tandem_type[tandem[[1]] == 1] <- "iadhore/"
tandem_type[tandem_mcl[[1]] == 1] <- paste(tandem_type[tandem_mcl[[1]] == 1], "iadhore_mcl/", sep = "")
tandem_type[tandem_mcscanx[[1]] == 1] <- paste(tandem_type[tandem_mcscanx[[1]] == 1], "mcscanx/", sep = "")
tandem_type[tandem_mcscanx_mcl[[1]] == 1] <- paste(tandem_type[tandem_mcscanx_mcl[[1]] == 1], "mcscanx_mcl/", sep = "")

## remove the final slash ("/")
tandem_type_vec <- lapply(as.vector(tandem_type), function(x) {
    paste(strsplit(x, split = "/")[[1]], collapse = "/")})
tandem_type_vec <- unlist(tandem_type_vec)

## write back to a matrix
tandem_type <- matrix(tandem_type_vec, ncol = 
    ncol(tandem_type), nrow = nrow(tandem_type))
rownames(tandem_type) <- colnames(tandem_type) <- rownames(tandem[[1]])

## plot tandem network
g <- igraph::graph_from_adjacency_matrix(tandem_sum, diag = FALSE, weighted = TRUE)
plot(g, vertex.label.cex = 0.1, vertex.size = 0.1, edge.width = 1, 
     edge.arrow.size = 0.1)

## get components of complete graphs, i.e. members of each graphs, in this case
## all tandem genes per region
components_g <- igraph::components(g)$membership
components_g_unique <- unique(components_g)

## check reliability of connection by checking number of neighbouring genes
res <- vector("list", length(bs_genes_all))
setwd("~/AG-Fernie/Thomas/Data/synteny/synteny_iadhore")
gff_files <- list.files()
gff_files <- gff_files[grep(gff_files, pattern = "gff")]

## iterate through bs_genes_all
for(i in 1:length(bs_genes_all) ) {
    species <- unlist(lapply(strsplit(bs_genes_all[i], split = "_"), "[", 1))
    ## get chromosome of gene
    gene <- unlist(lapply(strsplit(
        bs_genes_all[i], split = paste(species, "_", sep = "")), "[", 2))
    if (species == "quero") gene <- gsub(gene, pattern = "_P", replacement = "_T")
    gff <- read.table(gff_files[grep(gff_files, pattern = species)], 
        sep = "\t", stringsAsFactors = FALSE, quote = "")
    chr <- unique(gff[grep(gff[, 9], pattern = gene), 1])
    
    ## for species covsu and selmo strsplit the gene, since the genes were 
    ## pasted from different names (not directly from the gff files), get the 
    ## identifier that can be found in the gff
    if (species == "amahy") chr <- gsub(chr, pattern = "[|]", replacement = "_")
    if (species == "carpa") chr <- unique(gff[grep(gff[, 9],
        pattern = paste(gene, ";", sep = "")), 1])
    if (species=="covsu") chr <- unique(gff[grep(gff[, 9], 
        pattern = unlist(lapply(strsplit(gene, split = "_"), "[", 2))), 1])
    if (species=="selmo") chr <- unique(gff[grep(gff[, 9], 
        pattern = unlist(lapply(strsplit(gene, split = "_"), "[", 2))), 1])
    if (species == "salmi") chr <- unique(gff[grep(gff[, 9],
        pattern = paste(gene, ";", sep = "")), 1])
    if (species == "tripr") chr <- unique(gff[grep(gff[, 9],
        pattern = paste(gene, ".v2", sep = "")), 1])
    if (species == "vacco") chr <- unique(gff[grep(gff[, 9],
        pattern = paste("=", gene, ";", sep = "")), 1])
    ## load chromosome file lst and get position of gene in file
    res[[i]] <- list()
    if (!length(chr) == 0) {
        lst <- read.table(paste("lst_files/", species, "/", chr, ".lst", sep = ""), 
            stringsAsFactors = FALSE, quote="")
        lst <- unique(lst[, 1])
        genes_lst <- substring(lst, 1, nchar(lst)-1)
        if (species == "quero") genes_lst <- gsub(genes_lst, pattern = "_P", replacement = "_T")
        position <- which(genes_lst == paste(species, gene, sep = "_"))
        end <- length(lst)
        ## write the following information to the list res: start end position, name of chromosome
        res[[i]][[1]] <- c(position-1, end-position)
        res[[i]][[2]] <- chr
    } else {
        res[[i]][[1]] <- NULL
        res[[i]][[2]] <- chr
    }
    print(c(bs_genes_all[i], res[[i]][[1]], res[[i]][[2]]))
}

names(res) <- bs_genes_all

## save res to neighbours_on_chromosomes.RData
setwd("~/AG-Fernie/Federico/BAHD_SCPL_synteny/")
save(res, file = "./neighbours_on_chromosomes.RData")
load("neighbours_on_chromosomes.RData")
## end check reliability

## some plots
## calculate length of chromosome/scaffold where PKS gene is located on
length_chr <- unlist(lapply(res, function(x) sum(x[[1]]) + 1)) 
hist(log(length_chr))
## calculate minimum to end to chromosome/scaffold
min_chr <- unlist(lapply(res, function(x) min(x[[1]]))) 
hist(log(min_chr))

## check tandem components if they are reported correctly --> tandems should 
## be on the same chromosome and within distance of 20 genes
for (i in 1:length(components_g_unique)) {
    name_i <- names(which(components_g == i))
    chrs <- unlist(lapply(res[name_i], "[", 2))

    if (length(unique(chrs)) != 1) {print(i)} else { 
        ## print i when genes are on different chr
        position <- lapply(res[name_i], "[", 1)
        
        ## get positions and calculate differences 
        position <- unlist(position)[c(TRUE, FALSE)] 
        position <- sort(position)
        position_diff <- numeric(length(position) -1)
        for (j in 1:length(position_diff)) position_diff[j] <- position[j + 1] - position[j]
        if (any(position_diff > 20)) print(i)
    }
}


## name comb will be the vector that stores the pasted names of all tandem genes
name_comb <- vector("numeric", igraph::components(g)$no)
for (i in seq_along(components_g_unique)) {
    names_i_ind <- components_g == components_g_unique[i]
    names_i <- names(components_g)[names_i_ind]
    name_comb[names_i_ind] <- paste(names_i, collapse = "/")
}

##### Syntenic links #####
## add all bin_mat* matrices

## some checks for bin_mat* matrices, they should all have the same rownames
all(rownames(bin_mat) == rownames(bin_mat_mcl))
all(rownames(bin_mat) == rownames(bin_mat_mcscanx))
all(rownames(bin_mat) == rownames(bin_mat_mcscanx_mcl))
bin_mat_complete <- 0.25 * bin_mat + 0.25 * bin_mat_mcl + 
    0.25 * bin_mat_mcscanx + 0.25 * bin_mat_mcscanx_mcl

## get source and write the source to mat_type
mat_type <- matrix("", ncol = ncol(bin_mat), nrow = nrow(bin_mat))
mat_type[bin_mat == 1] <- "iadhore/"
mat_type[which(bin_mat_mcl == 1)] <- paste(
    mat_type[which(bin_mat_mcl == 1)], "iadhore_mcl/", sep = "")
mat_type[which(bin_mat_mcscanx == 1)] <- paste(
    mat_type[which(bin_mat_mcscanx == 1)], "mcscanx/", sep = "")
mat_type[which(bin_mat_mcscanx_mcl == 1)] <- paste(
    mat_type[which(bin_mat_mcscanx_mcl == 1)], "mcscanx_mcl/", sep = "")

## remove the final slash ("/")
mat_type_vec <- lapply(as.vector(mat_type), 
    function(x) paste(strsplit(x, split = "/")[[1]], collapse="/"))
mat_type_vec <- unlist(mat_type_vec)

## write back a to matrix
mat_type <- matrix(mat_type_vec, ncol = ncol(mat_type), nrow = nrow(mat_type))

## assign rownames of bin_mat to col- and rownames of bin_mat_complete and 
## mat_type
colnames(bin_mat_complete) <- rownames(bin_mat_complete) <- rownames(bin_mat)
rownames(mat_type) <- colnames(mat_type) <- rownames(bin_mat)

## plot bin_mat_complete network
g <- igraph::graph_from_adjacency_matrix(bin_mat_complete, diag=FALSE, weighted = TRUE)
plot(g, vertex.label.cex=0.1, vertex.size=0.1, edge.width=1, edge.arrow.size=0.1)

## rename tandem genes to pasted names
name <- names(components_g) ## name of tandem genes
name_comb_unique <- unique(name_comb)

## iterate through the unique name_comb
for (i in 1:length(name_comb_unique)) {
    
    name_i_ind <- name_comb_unique[[i]] == name_comb
    name_i <- name[name_i_ind]
    
    ## calculate the sum from all tandem genes 
    connection_sum <- apply(mat_type[name_i, ], 2, 
        function(x) paste(x, collapse = "/"))
    connection_sum <- paste(connection_sum, apply(mat_type[, name_i], 1, 
        function(x) paste(x, collapse = "/")), sep = "/")
    connection_sum <- lapply(strsplit(unlist(connection_sum), split = "/"), unique)                    
    connection_sum <- lapply(connection_sum, function(x) x[x != ""])
    connection_sum <- unlist(lapply(connection_sum, length))
    
    if (any(connection_sum>4)) stop(i)
    if (sum(connection_sum > 0) > 0) {
        ## combine all sources from all name_i and assign to mat_type
        if (sum(connection_sum>0) > 1) {
            mat_type_comb <- apply(mat_type[name_i, connection_sum > 0], 2, 
                function(x) paste(x, collapse = "/"))
            mat_type_comb <- paste(mat_type_comb, 
                apply(mat_type[connection_sum > 0, name_i], 1, 
                function(x) paste(x, collapse = "/")), sep = "/")
        } else {
            mat_type_comb <- paste(mat_type[name_i, connection_sum > 0], 
                collapse = "/")
            mat_type_comb <- paste(mat_type_comb, paste(
                mat_type[connection_sum > 0, name_i], collapse = "/"), sep = "/")
        }
        ## split by /
        mat_type_comb_l <- strsplit(mat_type_comb, split = "/") 
        ## make unique and sort
        mat_type_comb_l <- lapply(mat_type_comb_l, function(x) sort(unique(x))) 
        ## remove ""
        mat_type_comb_l <- lapply(mat_type_comb_l, function(x) x[x != ""])
        ## paste again
        mat_type_comb_l <- lapply(mat_type_comb_l, function(x) 
            paste(x, collapse = "/"))
        mat_type_comb <- as.vector(unlist(mat_type_comb_l))
        mat_type[name_i, connection_sum > 0] <- matrix(
            rep(mat_type_comb, times = length(name_i)), 
            ncol = sum(connection_sum > 0), byrow = T)
        mat_type[connection_sum > 0, name_i] <- matrix(
            rep(mat_type_comb, times = length(name_i)), 
            nrow = sum(connection_sum > 0), byrow = F)
        
        ## assign connection_sum to the first element in name i
        bin_mat_complete[name_i, ] <- matrix(
            rep(connection_sum / 4, times = length(name_i)), 
            ncol = length(connection_sum), byrow = T)
        bin_mat_complete[, name_i] <- matrix(
            rep(connection_sum / 4, times = length(name_i)), 
            nrow = length(connection_sum), byrow = F)
    }
    
    ## rename the first element to the combined name
    rownames(bin_mat_complete)[rownames(bin_mat_complete) %in% name_i] <- name_comb_unique[[i]] 
    colnames(bin_mat_complete)[colnames(bin_mat_complete) %in% name_i] <- name_comb_unique[[i]]
    rownames(mat_type)[rownames(mat_type) %in% name_i] <- name_comb_unique[[i]]
    colnames(mat_type)[colnames(mat_type) %in% name_i] <- name_comb_unique[[i]]
}

## remove all other elements from bin_mat_complete
ind_remove <- duplicated(rownames(bin_mat_complete)) 
## remove the duplicated rownames
bin_mat_complete <- bin_mat_complete[!ind_remove, !ind_remove]
mat_type <- mat_type[!ind_remove, !ind_remove]

## remove type of connection that are "iadhore_mcl/mcscanx", 
## "iadhore/mxscanx_mcl" (not compatible techniques and clustering)
table(mat_type)
bin_mat_complete[which(
    mat_type %in% c( "iadhore_mcl/mcscanx", "iadhore/mcscanx_mcl"))] <- 0
mat_type[which(
    mat_type %in% c("iadhore_mcl/mcscanx", "iadhore/mcscanx_mcl"))] <- ""


setwd("~/AG-Fernie/Federico/BAHD_SCPL_synteny")
save(bin_mat_complete, mat_type, 
    file = "BAHD_SCPL_bin_mat_complete_mat_type.RData")

## plot the bin_mat_complete network
g <- graph_from_adjacency_matrix(bin_mat_complete, diag = FALSE, mode = "directed", weighted = TRUE)
plot(g, vertex.label.cex = 0.1, vertex.size = 0.1, edge.width = 1, 
     edge.arrow.size = 0.1)


## remove the proteins that do not link to others
inds_keep <- apply(bin_mat_complete, 1, sum) > 0
## keep also duplicated genes
inds_keep[rownames(bin_mat_complete) %in% name_comb_unique ] <- TRUE

## check how many genes are on the chromosome/scaffold for (not) linking ones 
## not linking ones
res_nlink <- res[rownames(bin_mat_complete)[!inds_keep]]
res_nlink <- lapply(1:length(res_nlink), function(x) res_nlink[[x]][[1]])
res_nlink_sum <- lapply(res_nlink, sum)
names(res_nlink_sum) <- rownames(bin_mat_complete)[!inds_keep]

## add 1 since the pks_genes is not counted
res_nlink_sum <- unlist(res_nlink_sum) + 1 

## create a vector that stores the colour depending on how many genes there 
## are on the chromosome/scaffold
reliability <- vector("character", dim(bin_mat_complete)[1])
names(reliability) <- rownames(bin_mat_complete)
tmp <- unlist(res_nlink_sum)
reliability[names(res_nlink_sum[tmp < 15])] <- "red"
reliability[names(res_nlink_sum[tmp >= 15 & tmp < 25])] <- "yellow"
reliability[names(res_nlink_sum[tmp >= 25])] <- "green"

## plot network that contains features that do not link to other
net <- graph_from_adjacency_matrix(
    bin_mat_complete[!inds_keep, !inds_keep], mode = "directed", diag = FALSE, weighted = TRUE)
plot(net, vertex.label.cex = 0.3, vertex.size = 5, 
    vertex.color = reliability[!inds_keep],edge.arrow.size = 0.1)

## build network for linking features
res_link <- res[names(res) %in% unlist(
    strsplit(rownames(bin_mat_complete)[inds_keep], split = "/"))]
res_link <- lapply(1:length(res_link), function(x) res_link[[x]][[1]])
res_link_sum <- lapply(res_link, sum)
names(res_link_sum) <- names(res[names(res) %in% unlist(
    strsplit(rownames(bin_mat_complete)[inds_keep], split = "/"))])

## rename res_link_sum and truncate with names used in bin_mat_complete
for (i in 1:length(components_g_unique)) {
    new_name <- paste(names(components_g)[components_g == components_g_unique[i]], collapse="/")
    inds <- grep(names(res_link_sum), pattern = gsub("/", "|", new_name))
    names(res_link_sum)[inds[1]] <- new_name
    res_link_sum <- res_link_sum[ -inds[ 2:length(inds) ]]
}
## add 1 since the pks_genes is not counted
res_link_sum <- unlist(res_link_sum) + 1 


write.table(
    data.frame(names = names(unlist(res_link_sum)), value = unlist(res_link_sum)), 
    file = "./figure_synteny_network_quality/sum_genes_scaffold_linking.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(
    data.frame(names = names(unlist(res_nlink_sum)), value = unlist(res_nlink_sum)), 
    file = "./figure_synteny_network_quality/sum_genes_scaffold_notlinking.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE)

## plot for figure synteny_network_quality
df <- unlist(lapply(strsplit(names(res_nlink_sum), split = "_"), "[", 1))
df <- sort(table(df))
df <- data.frame(species = names(df), number = as.vector(df))
df$species <- factor(df$species, levels=df$species[order(df$number)])
g <- ggplot(df) + geom_bar(aes(x = species, y = number), stat = "identity") + 
    ylim(0,575) + theme_bw() + 
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(g, 
    file = "./figure_synteny_network_quality/synteny_network_quality_notlinking_species_hist.pdf")

g <- ggplot(data.frame(value = unlist(res_link_sum)), aes(x = value)) + 
    ylim(-100, 1000) + xlim(0, 15000) + theme_bw() + geom_histogram(binwidth = 100) + 
    theme(panel.grid = element_blank())
ggsave(g, file = "./figure_synteny_network_quality/sum_genes_scaffold_linking_hist.pdf")

g <- ggplot(data.frame(value = unlist(res_nlink_sum)), aes(x = value)) + 
    ylim(-1,1000) + xlim(0, 15000) + theme_bw() + geom_histogram(binwidth = 100) + 
    theme(panel.grid = element_blank())
ggsave(g, 
       file = "./figure_synteny_network_quality/sum_genes_scaffold_notlinking_hist.pdf")

## set colour according to the number of neightbour genes
tmp <- unlist(res_link_sum)
reliability[names(tmp[tmp < 15])] <- "red"
reliability[names(tmp[tmp >= 15 & tmp < 25])] <- "yellow"
reliability[names(tmp[tmp >= 25])] <- "green"

## plot network that contains features that link to other
bin_mat_complete_cut <- bin_mat_complete[inds_keep, inds_keep]
mat_type_cut <- mat_type[inds_keep, inds_keep]
net <- graph_from_adjacency_matrix(bin_mat_complete_cut, weighted = TRUE, 
    mode = "undirected", diag = FALSE)
plot(net, vertex.label.cex = 0.1, vertex.size = 5, 
     vertex.color = reliability[inds_keep], edge.arrow.size = 0.1)

## plot type of links 
df <- data.frame(names = names(unlist(sort(table(mat_type_cut[lower.tri(mat_type_cut)])[-1]))), 
    value = as.vector(sort(table(mat_type_cut[lower.tri(mat_type_cut)])[-1])))
df$names <- factor(df$names, levels = df$names[order(df$value)])
g <- ggplot(df) + geom_bar(aes(x = names, y = value), stat = "identity") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(0, 120000)
ggsave(g, file = "./figure_synteny_network_quality/type_link_synteny.pdf")

## export the network
setwd("~/AG-Fernie/Federico/")
save(bin_mat_complete, bin_mat_complete_cut, mat_type_cut, file = "bin_mat_complete_cut_tandemgenes.RData")
load("bin_mat_complete_cut_tandemgenes.RData")

