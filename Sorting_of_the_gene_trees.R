# Cesta
setwd("/Users/marcosqueiroz/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/0_trees")

# Knihovny
library("ape")
library("phangorn")
library("stringr")
library("readr")
library("apex")
library("phylotools")
library("msa")
library("seqinr")



## List only trees
trees <- Sys.glob('*.contree')

#For loop to get all trees and save them in a single file
#for each file tree
for(i in 1:length(trees)){
  #get the name of the file, and remove the extensions (i.e., the gene name)
  tree_name <- str_replace(trees[i], '\\.fas.contree', '')
  #read the tree as a text file
  my_tree <- read_file(trees[i])
  #paste the name of the tree and the tree itself to a variable
  name_and_tree <- paste0(tree_name, my_tree)
  #write (append) the tree name and the tree itself to a file
  write(name_and_tree, file = 'all_trees.txt', append = TRUE)
}

#import the tree
all_trees <- read.tree(file = 'all_trees.txt')

#drop samples:
#'Curcuma-sp_Z54', Curcuma-sp_S662, Curcuma-manggaAndamans_C005 
all_trees <- lapply(all_trees, drop.tip, tip = 'Curcuma-sp_Z54')
all_trees <- lapply(all_trees, drop.tip, tip = 'Curcuma-sp_S662')
all_trees <- lapply(all_trees, drop.tip, tip = 'Curcuma-manggaAndamans_C005')
class(all_trees) <- 'multiPhylo'

#Parental tree 1 - all trees where curcuma and vamana are monophyletic

#tips with curcuma and vamana
tips_c_v <- c('Curcuma-mangga_S669', 'Curcuma-mangga_S670', 'Curcuma-amada_S673', 'Curcuma-amada_S580', 'Curcuma-sp_S659', 'Curcuma-sp_S660', 'Curcuma-sulphurea_S576', 'Curcuma-plicata_C343', 'Curcuma-mangga_S677', 'Curcuma-amada_S577', 'Curcuma-amada_S672', 'Curcuma-mangga_S581', 'Curcuma-amada_S671', 'Curcuma-aromatica_S547', 'Curcuma-aromatica_S561', 'Curcuma-aromatica_S528', 'Curcuma-picta_S514', 'Curcuma-rubescens_S661', 'Curcuma-repens_S516', 'Curcuma-rubescens_S578', 'Curcuma-montana_S544', 'Curcuma-montana_S530', 'Curcuma-montana_S531', 'Curcuma-viridiflora_S572', 'Curcuma-prakasha_S560', 'Curcuma-petiolata_S664', 'Curcuma-sp_S666', 'Curcuma-prakasha_S527', 'Curcuma-prakasha_S511', 'Curcuma-prakasha_S512', 'Curcuma-sp_S563', 'Curcuma-sp_C342', 'Curcuma-plicata_C295', 'Curcuma-plicata_S204', 'Curcuma-plicata_S657', 'Curcuma-rubrobracteata_S579', 'Curcuma-rubrobracteata_S559', 'Curcuma-sp_S668', 'Curcuma-sp_S665', 'Curcuma-sp_S667', 'Curcuma-aurantiaca_S675', 'Curcuma-reclinata_S678', 'Curcuma-sulcata_S606', 'Curcuma-inodora_S602', 'Curcuma-decipiens_S600', 'Curcuma-decipiens_S601', 'Curcuma-pseudomontana_S526', 'Curcuma-pseudomomtana_S605', 'Curcuma-caulina_S540', 'Curcuma-cannanorensis_S164', 'Curcuma-karnatakensis_S676', 'Curcuma-mutabilis_S604', 'Curcuma-bhatii_S165', 'Curcuma-scaposa_S541', 'Curcuma-angustifolia_S582', 'Curcuma-angustifolia_S575', 'Curcuma-angustifolia_S529', 'Curcuma-angustifolia_S674', 'Curcuma-angustifolia_S546', 'Curcuma-candida_S539', 'Curcuma-candida_S156', 'Curcuma-vamana_S167')

#for each tree
for(t in 1:length(all_trees)){
  #check if that tree is monophyletic, following our tips
  #the try here is to avoid any break in the loop, since the is.monophyletic can stop if any error appears
  a <- try(is.monophyletic(phy=all_trees[[t]], tips = tips_c_v))
  #if our tree is monophyletic
  if(a == TRUE){
    #get the tree name
    tree_name <- names(all_trees[t])
    #get the tree
    my_tree <- all_trees[[t]]
    #write a file with all tree names (i.e., the gene name)
    write(tree_name, file = 'gene_names_c_v.txt', append = TRUE)
    #write a file with all monophyletic trees, including their names
    write.tree(phy = my_tree, file = './parental_c_v/parental_c_v.nwk', append = TRUE, tree.names = tree_name)
    #print(t)
  }
}
#import the tree
parental_c_v <- read.tree('./parental_c_v/parental_c_v.nwk')

#previous function - it stops when some error shows up
#parental_c_v <- all_trees[unlist(lapply(X=all_trees, FUN=is.monophyletic, tips=tips_c_v))]

#Parental tree 2 - all trees where curcuma, vamana and roscoeana are monophyletic
tips_c_v_r <- c('Curcuma-mangga_S669', 'Curcuma-mangga_S670', 'Curcuma-amada_S673', 'Curcuma-amada_S580', 'Curcuma-sp_S659', 'Curcuma-sp_S660', 'Curcuma-sulphurea_S576', 'Curcuma-plicata_C343', 'Curcuma-mangga_S677', 'Curcuma-amada_S577', 'Curcuma-amada_S672', 'Curcuma-mangga_S581', 'Curcuma-amada_S671', 'Curcuma-aromatica_S547', 'Curcuma-aromatica_S561', 'Curcuma-aromatica_S528', 'Curcuma-picta_S514', 'Curcuma-rubescens_S661', 'Curcuma-repens_S516', 'Curcuma-rubescens_S578', 'Curcuma-montana_S544', 'Curcuma-montana_S530', 'Curcuma-montana_S531', 'Curcuma-viridiflora_S572', 'Curcuma-prakasha_S560', 'Curcuma-petiolata_S664', 'Curcuma-sp_S666', 'Curcuma-prakasha_S527', 'Curcuma-prakasha_S511', 'Curcuma-prakasha_S512', 'Curcuma-sp_S563', 'Curcuma-sp_C342', 'Curcuma-plicata_C295', 'Curcuma-plicata_S204', 'Curcuma-plicata_S657', 'Curcuma-rubrobracteata_S579', 'Curcuma-rubrobracteata_S559', 'Curcuma-sp_S668', 'Curcuma-sp_S665', 'Curcuma-sp_S667', 'Curcuma-aurantiaca_S675', 'Curcuma-reclinata_S678', 'Curcuma-sulcata_S606', 'Curcuma-inodora_S602', 'Curcuma-decipiens_S600', 'Curcuma-decipiens_S601', 'Curcuma-pseudomontana_S526', 'Curcuma-pseudomomtana_S605', 'Curcuma-caulina_S540', 'Curcuma-cannanorensis_S164', 'Curcuma-karnatakensis_S676', 'Curcuma-mutabilis_S604', 'Curcuma-bhatii_S165', 'Curcuma-scaposa_S541', 'Curcuma-angustifolia_S582', 'Curcuma-angustifolia_S575', 'Curcuma-angustifolia_S529', 'Curcuma-angustifolia_S674', 'Curcuma-angustifolia_S546', 'Curcuma-candida_S539', 'Curcuma-candida_S156', 'Curcuma-candida_C346', 'Curcuma-vamana_S167', 'Curcuma-roscoeana_S656', 'Curcuma-roscoeana_S593', 'Curcuma-roscoeana_S163', 'Curcuma-roscoeana_C400', 'Curcuma-myanmarensis_S168') 

for(t in 1:length(all_trees)){
  a <- try(is.monophyletic(phy=all_trees[[t]], tips = tips_c_v_r))
  if(a == TRUE){
    tree_name <- names(all_trees[t])
    my_tree <- all_trees[[t]]
    write(tree_name, file = 'gene_names_c_v_r.txt', append = TRUE)
    write.tree(phy = my_tree, file = './parental_c_v_r/parental_c_v_r.nwk', append = TRUE, tree.names = tree_name)
    #print(t)
  }
}

parental_c_v_r <- read.tree('./parental_c_v_r/parental_c_v_r.nwk')
#parental_c_v_r <- all_trees[unlist(lapply(X=all_trees, FUN=is.monophyletic, tips=tips_c_v_r))]

#Parental tree 3 - all tree where vamana, roscoeana and hitcheniopsis are monophyletic


tips_v_r_h <- c('Curcuma-thorelii_S590', 'Curcuma-thorellii_C250', 'Curcuma-parviflora_S592', 'Curcuma-pygmaea_Z48', 'Curcuma-gracillima_S1', 'Curcuma-prasina_Z1145', 'Curcuma-harmandii_C199', 'Curcuma-rhabdota_C334', 'Curcuma-rhabdota_S599', 'Curcuma-sparganiifolia_Z1140', 'Curcuma-sparganifolia_S594', 'Curcuma-alismatifolia_S591', 'Curcuma-alismatifolia_S158', 'Curcuma-involucrata_S205', 'Curcuma-campanulata_S598', 'Curcuma-leonidii_Z1146', 'Curcuma-parviflora_S159', 'Curcuma-parviflora_S597', 'Curcuma-parviflora_S595', 'Curcuma-vamana_S167', 'Curcuma-roscoeana_S656', 'Curcuma-roscoeana_S593', 'Curcuma-roscoeana_S163', 'Curcuma-roscoeana_C400', 'Curcuma-myanmarensis_S168')

for(t in 1:length(all_trees)){
  a <- try(is.monophyletic(phy=all_trees[[t]], tips = tips_v_r_h))
  if(a == TRUE){
    tree_name <- names(all_trees[t])
    my_tree <- all_trees[[t]]
    write(tree_name, file = 'gene_names_v_r_h.txt', append = TRUE)
    write.tree(phy = my_tree, file = './parental_v_r_h/parental_v_r_h.nwk', append = TRUE, tree.names = tree_name)
    #print(t)
  }
}

parental_v_r_h <- read.tree('./parental_v_r_h/parental_v_r_h.nwk')

############




############
#Copy the gene alignments for specific folders

gene_names_c_v <- read.table('./parental_c_v/gene_names_c_v.txt')
gene_names_c_v <- paste0('../0_alignments/all_genes/', gene_names_c_v$V1, '.fas')
gene_names_c_v_r <- read.table('./parental_c_v_r/gene_names_c_v_r.txt')
gene_names_c_v_r <- paste0('../0_alignments/all_genes/', gene_names_c_v_r$V1, '.fas')
gene_names_v_r_h <- read.table('./parental_v_r_h/gene_names_v_r_h.txt')
gene_names_v_r_h <- paste0('../0_alignments/all_genes/', gene_names_v_r_h$V1, '.fas')

#copy the files to a new directory
c_v_directory <- '../0_alignments/parental_alignments/parental_c_v/'
#file.copy(gene_names_c_v, c_v_directory)
c_v_files <- list.files(c_v_directory, full.names = TRUE)

c_v_r_directory <- '../0_alignments/parental_alignments/parental_c_v_r/'
#file.copy(gene_names_c_v_r, c_v_r_directory)
c_v_r_files <- list.files(c_v_r_directory, full.names = TRUE)

v_r_h_directory <- '../0_alignments/parental_alignments/parental_v_r_h/'
#file.copy(gene_names_v_r_h, v_r_h_directory)
v_r_h_files <- list.files(v_r_h_directory, full.names = TRUE)

#drop the undesired samples

undesired <- c('Curcuma-sp_Z54', 'Curcuma-sp_S662', 'Curcuma-manggaAndamans_C005')

for(g in 1:length(c_v_files)){
  rm.sequence.fasta(infile = c_v_files[g],
                    outfile = c_v_files[g],
                    to.rm = undesired)
}

for(g in 1:length(c_v_r_files)){
  rm.sequence.fasta(infile = c_v_r_files[g],
                    outfile = c_v_r_files[g],
                    to.rm = undesired)
}

for(g in 1:length(v_r_h_files)){
  rm.sequence.fasta(infile = v_r_h_files[g],
                    outfile = v_r_h_files[g],
                    to.rm = undesired)
}











#Alignment - i did not run this part... the genes are already aligned!

#function to export the msa alignment to a fasta file
alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

#msa uses clustalw as default

#c_v

for(g in 1:length(c_v_files)){
  gene <- readDNAStringSet(c_v_files[g])
  gene_name <- basename(c_v_files[g])
  gene_name <- str_replace(gene_name, '.fas', '')
  align_gene <- msa(gene)
  alignment2Fasta(alignment = align_gene, filename = paste0(dirname(c_v_files[g]), '/', gene_name, '_aligned.fas'))
}






############
#Concatenation

gene_names_c_v <- as.data.frame(list.files('../0_alignments/parental_alignments/parental_c_v/'))
colnames(gene_names_c_v) <- 'cv'
gene_names_c_v <- paste0('../0_alignments/parental_alignments/parental_c_v/', gene_names_c_v$cv)

gene_names_c_v_r <- as.data.frame(list.files('../0_alignments/parental_alignments/parental_c_v_r/'))
colnames(gene_names_c_v_r) <- 'cvr'
gene_names_c_v_r <- paste0('../0_alignments/parental_alignments/parental_c_v_r/', gene_names_c_v_r$cvr)

gene_names_v_r_h <- as.data.frame(list.files('../0_alignments/parental_alignments/parental_v_r_h/'))
colnames(gene_names_v_r_h) <- 'vrh'
gene_names_v_r_h <- paste0('../0_alignments/parental_alignments/parental_v_r_h/', gene_names_v_r_h$vrh)


#read the files as multiFASTA
gene_c_v <- read.multiFASTA(gene_names_c_v)
gene_c_v_r <- read.multiFASTA(gene_names_c_v_r)
gene_v_r_h <- read.multiFASTA(gene_names_v_r_h)

#concatenate them
gene_c_v <- concatenate(gene_c_v)
gene_c_v_r <- concatenate(gene_c_v_r)
gene_v_r_h <- concatenate(gene_v_r_h)

#write the concatenated files
write.FASTA(gene_c_v, file = '../0_alignments/parental_alignments/parental_c_v/parental_c_v_concat.fasta')
write.FASTA(gene_c_v_r, file = '../0_alignments/parental_alignments/parental_c_v_r/parental_c_v_r_concat.fasta')
write.FASTA(gene_v_r_h, file = '../0_alignments/parental_alignments/parental_v_r_h/parental_v_r_h_concat.fasta')

#Maximum Likelihood trees were created using MEGA, with the default parameters.

#I used the parental_x_y_concat.fasta files
#The Maximum Likelihood trees were exported in Newick format as parental_x_y_ML.nwk

#For the dated tree, I used the RelTime function in MEGA (see the Mello, 2018 paper from Molecular Biology and Evolution)
#Globba was used as an outgroup
#Camptandra was used as the constrain point, with a normal distribution of mean 15, and 3.0 of standard deviation

#The outputs are Newick trees - parental_x_y_ML_calibrated.nwk


############

#For biogeobears, we need to drop some tips (i.e., removing redundant assessions, keeping only one individual per species)
#We only kept the tips that are present in the Locality_table_for_BioGeoBears table
#We are going to use the tree generated in RelTime in MEGA




parental_c_v_calib <- read.tree(file = './parental_c_v/parental_c_v_calibrated.nwk')
parental_c_v_r_calib <- read.tree(file = './parental_c_v_r/parental_c_v_r_calibrated.nwk')
parental_v_r_h_calib <- read.tree(file = './parental_v_r_h/parental_v_r_h_calibrated.nwk')

#I can't extract the trees with Globba... so I'll not keep it.
tips_to_keep <- c('Curcuma-mangga_S669', 'Curcuma-sp_S659', 'Curcuma-sp_S660', 'Curcuma-sulphurea_S576', 'Curcuma-mangga_S677', 'Curcuma-aromatica_S547', 'Curcuma-picta_S514', 'Curcuma-rubescens_S661', 'Curcuma-montana_S544', 'Curcuma-viridiflora_S572', 'Curcuma-petiolata_S664', 'Curcuma-prakasha_S527', 'Curcuma-sp_S563', 'Curcuma-sp_C342', 'Curcuma-rubrobracteata_S579', 'Curcuma-sp_S668', 'Curcuma-sp_S667', 'Curcuma-aurantiaca_S675', 'Curcuma-reclinata_S678', 'Curcuma-sulcata_S606', 'Curcuma-inodora_S602', 'Curcuma-decipiens_S600', 'Curcuma-pseudomontana_S526', 'Curcuma-caulina_S540', 'Curcuma-cannanorensis_S164', 'Curcuma-karnatakensis_S676', 'Curcuma-mutabilis_S604', 'Curcuma-bhatii_S165', 'Curcuma-scaposa_S541', 'Curcuma-angustifolia_S582', 'Curcuma-candida_S539', 'Curcuma-thorelii_S590', 'Curcuma-pygmaea_Z48', 'Curcuma-gracillima_S1', 'Curcuma-prasina_Z1145', 'Curcuma-harmandii_C199', 'Curcuma-rhabdota_C334', 'Curcuma-sparganiifolia_Z1140', 'Curcuma-alismatifolia_S591', 'Curcuma-involucrata_S205', 'Curcuma-leonidii_Z1146', 'Curcuma-parviflora_S159', 'Curcuma-vamana_S167', 'Curcuma-roscoeana_S656', 'Curcuma-myanmarensis_S168', 'Curcuma-rhomba_S160', 'Curcuma-cotuana_Z1144', 'Curcuma-pambrosima_S319', 'Curcuma-vitellina_S3', 'Curcuma-singularis_C253', 'Curcuma-newmanii_S314', 'Curcuma-xanthella_S320', 'Curcuma-pierreana_S596', 'Curcuma-singularis_S2', 'Curcuma-sp_C285', 'Curcuma-arida_Z1143', 'Curcuma-arida_C263', 'Curcuma-kirirom_Z1142', 'Curcuma-rhomba_C256', 'Curcuma-cochinchinensis_C147', 'Curcuma-putii_S489', 'Smithatris-supraneanae_S542', 'Curcuma-flammea_S318', 'Curcuma-glans_Z176', 'Curcuma-corniculata_S655', 'Curcuma-ecomata_S658', 'Curcuma-cochinchinensis_C255', 'Curcuma-flaviflora_S315')
 
parental_c_v_kept <- keep.tip(parental_c_v_calib, tip = tips_to_keep)
parental_c_v_r_kept <- keep.tip(parental_c_v_r_calib, tip = tips_to_keep)
parental_v_r_h_kept <- keep.tip(parental_v_r_h_calib, tip = tips_to_keep)

#Rename the tips
#BioGeoBears doesn't like special characters, so we will replace every hyphen (-)
#with a underscore (_)

for(t in 1:length(parental_c_v_kept$tip.label)){
  #print(paste(t, parental_c_v_kept$tip.label[t], sep = ' '))
  old_tip <- parental_c_v_kept$tip.label[t]
  new_tip <- str_replace(old_tip, '-', '_')
  parental_c_v_kept$tip.label[parental_c_v_kept$tip.label==old_tip] <- new_tip
}

for(t in 1:length(parental_c_v_r_kept$tip.label)){
  #print(paste(t, parental_c_v_r_kept$tip.label[t], sep = ' '))
  old_tip <- parental_c_v_r_kept$tip.label[t]
  new_tip <- str_replace(old_tip, '-', '_')
  parental_c_v_r_kept$tip.label[parental_c_v_r_kept$tip.label==old_tip] <- new_tip
}

for(t in 1:length(parental_v_r_h_kept$tip.label)){
  #print(paste(t, parental_v_r_h_kept$tip.label[t], sep = ' '))
  old_tip <- parental_v_r_h_kept$tip.label[t]
  new_tip <- str_replace(old_tip, '-', '_')
  parental_v_r_h_kept$tip.label[parental_v_r_h_kept$tip.label==old_tip] <- new_tip
}

write.tree(parental_c_v_kept, file = './parental_c_v/parental_c_v_kept.nwk')
write.tree(parental_c_v_r_kept, file = './parental_c_v_r/parental_c_v_r_kept.nwk')
write.tree(parental_v_r_h_kept, file = './parental_v_r_h/parental_v_r_h_kept.nwk')

















#############Eliska's part >>>

snapp_trees <- read.tree(file="snapp_with_hybrids.trees.nwk")

snapp_trees <- read.tree(file="snapp.trees.nwk")

stromy.mono.s1 <- snapp_trees[unlist(lapply(X=snapp_trees, FUN=is.monophyletic, tips=c("bhatii_S165", "cannanorensis_S164", "myanmarensis_S168","roscoeana_S163", "vamana_S167")))]
stromy.mono.s1
write.tree(phy=stromy.mono.s1, file="stromy_mono_s1.nwk", tree.names=TRUE)
stromy.mono.s2 <- snapp_trees[unlist(lapply(X=snapp_trees, FUN=is.monophyletic, tips=c("alismatifolia_S158", "parviflora_S159", "myanmarensis_S168","roscoeana_S163", "vamana_S167")))]
stromy.mono.s2
write.tree(phy=stromy.mono.s2, file="stromy_mono_s2.nwk", tree.names=TRUE)











# Stromy
T1 <- read.tree(file="tree1_32percent.nwk")
T1
T7 <- read.tree(file="tree7_2percent.nwk")
T7
snapp_trees <- read.tree(file="snapp.trees_ed.nwk")
names(snapp_trees) <- paste("snapp_tree_", seq(from=1, to=length(snapp_trees), by=1), sep="")
snapp_trees

## Vzdálenosti mezi stromy

# T1
# Penny and Hendy (1985, originally from Robinson and Foulds 1981) nebo Robinson-Foulds distance (only depends on the topology)
dists.T1 <- rep(x=NA, times=length(snapp_trees))
for (T in 1:length(snapp_trees)) {
  dists.T1[T] <- dist.topo(x=T1, y=snapp_trees[[T]], method="PH85") # Penny and Hendy (1985, originally from Robinson and Foulds 1981)
  # 	dists.T1[T] <- RF.dist(tree1=T1, tree2=snapp_trees[[T]], normalize=FALSE, check.labels=TRUE, rooted=TRUE) # Robinson-Foulds distance (only depends on the topology)
}
names(dists.T1) <- names(snapp_trees)
dists.T1
png(filename="dists_t1.png", width=1000, height=3000)
barplot(height=sort(x=dists.T1), horiz=TRUE, main="T1", las=1)
dev.off()

# T7
# Penny and Hendy (1985, originally from Robinson and Foulds 1981) nebo Robinson-Foulds distance (only depends on the topology)
dists.T7 <- rep(x=NA, times=length(snapp_trees))
for (T in 1:length(snapp_trees)) {
  dists.T7[T] <- dist.topo(x=T7, y=snapp_trees[[T]], method="PH85") # Penny and Hendy (1985, originally from Robinson and Foulds 1981)
  # 	dists.T7[T] <- RF.dist(tree1=T7, tree2=snapp_trees[[T]], normalize=FALSE, check.labels=TRUE, rooted=TRUE) # Robinson-Foulds distance (only depends on the topology)
}
names(dists.T7) <- names(snapp_trees)
dists.T7
png(filename="dists_t7.png", width=1000, height=3000)
barplot(height=sort(x=dists.T7), horiz=TRUE, main="T7", las=1)
dev.off()

# Srovnání
dists.srovnani <- cbind(dists.T1, dists.T7, rep(x=NA, times=length(snapp_trees)), rep(x=NA, times=length(snapp_trees)), rep(x=NA, times=length(snapp_trees)))
colnames(dists.srovnani) <- c("dists.T1", "dists.T7", "T1==T7", "T1<T7", "T1>T7")
dists.srovnani[,3] <- dists.T1==dists.T7
dists.srovnani[,4] <- dists.T1<dists.T7
dists.srovnani[,5] <- dists.T1>dists.T7
head(dists.srovnani)
write.table(x=dists.srovnani, file="vzdalenosti.tsv", quote=FALSE, sep="\t")

# Geodesic distance

# Načtení
stromy <- read.tree(file="vse.nwk")
names(stromy) <- c("T1", "T7", paste("snapp_tree_", seq(from=1, to=length(snapp_trees), by=1), sep=""))
stromy

# Distance
dist.geo <- dist.multiPhylo(x=stromy, method="geodesic")
dist.geo <- as.matrix(dist.geo)
colnames(dist.geo) <- rownames(dist.geo) <- names(stromy)
dist.geo.srovnani <- cbind(t(dist.geo[1:2,3:812]), t(dist.geo[1:2,3:812])[,1]==t(dist.geo[1:2,3:812])[,2], t(dist.geo[1:2,3:812])[,1]<t(dist.geo[1:2,3:812])[,2], t(dist.geo[1:2,3:812])[,1]>t(dist.geo[1:2,3:812])[,2])
colnames(dist.geo.srovnani) <- c("dists.T1", "dists.T7", "T1==T7", "T1<T7", "T1>T7")
head(dist.geo.srovnani)
write.table(x=dist.geo.srovnani, file="vzdalenosti_geo.tsv", quote=FALSE, sep="\t")
png(filename="dists_geo_t1.png", width=1000, height=3000)
barplot(height=sort(x=dist.geo.srovnani[,1]), horiz=TRUE, main="T1", las=1)
dev.off()
png(filename="dists_geo_t7.png", width=1000, height=3000)
barplot(height=sort(x=dist.geo.srovnani[,2]), horiz=TRUE, main="T7", las=1)
dev.off()

