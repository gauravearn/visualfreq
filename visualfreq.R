#!/usr/bin/R
# Author Gaurav Sablok
# Universitat Potsdam
# Date 2024-3-21
# a gene alignment to the visual alignment priting 
# for the creation of the alignment maps
visualfreq <- function(inputfile){
suppressPackageStartupMessages(library(argparser, pos = "package:base"))
suppressPackageStartupMessages(library(methods, pos = "package:base"))
suppressPackageStartupMessages(library(Biostrings, pos = "package:base"))
suppressPackageStartupMessages(library(ggmsa, pos = "package:base"))
suppressPackageStartupMessages(library(ggplot2, pos = "package:base"))
suppressPackageStartupMessages(library(ape, pos = "package:base"))
suppressPackageStartupMessages(library(odseq, pos = "package:base"))
suppressPackageStartupMessages(library(venn, pos = "package:base"))
suppressPackageStartupMessages(library(msa, pos = "package:base"))
suppressPackageStartupMessages(library(DECIPHER, pos = "package:base"))
suppressPackageStartupMessages(library(ggtree, pos = "package:base"))
suppressPackageStartupMessages(library(phytools, pos = "package:base"))
suppressPackageStartupMessages(library(reticulate, pos = "package:base"))
  alignment <- DNAMultipleAlignment(msaMuscle(readDNAStringSet(file = inputfile)))
  alignmentwrite <- msaConvert(alignment, type = "ape::DNAbin")
  write.FASTA(alignmentwrite, file = "alignment.fasta")
  fasta <- "/home/gaurav/Desktop/arabidopsis/alignment.fasta"
  ggmsa(fasta, 164, 213, color="Chemistry_AA")
  complete_likelihood <- TreeLine(myXStringSet = 
          readDNAStringSet(file = "alignment.fasta"),myDistMatrix = 
         DistanceMatrix(RemoveGaps(StaggerAlignment(AlignSeqs(readDNAStringSet(
                     file = "alignment.fasta"))), 
                           removeGaps = "all", processors = 3), type = "dist"), 
       method = "complete", cutoff = 0.05, showPlot = TRUE, reconstruct = TRUE)
  WriteDendrogram(complete_likelihood, file = "maximum_likelihood.txt")
  ggmsa(DNAMultipleAlignment(msaMuscle(readDNAStringSet(
 file = "alignment.fasta"))), 
 color = "Shapely_NT", font = "DroidSansMono", 
 char_width = 0.5,seq_name = TRUE + geom_seqlogo(color = "Shapely_NT")) 
  +geom_msaBar() + geom_GC()
  ggsave(alignment_plot, width = 4, height = 4)
  sink(file = "distance_matrix_arabiopsis.txt")
  dist.dna(read.FASTA(file = "alignment.fasta", type = "DNA"))
  sink()
}
