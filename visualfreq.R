#!/usr/bin/R
# Author Gaurav Sablok
# Universitat Potsdam
# Date 2024-3-21
# a gene alignment to the visual alignment priting 
# for the creation of the alignment maps
visualfreq <- function(inputfile, alignmentstart, alignmentend, outputmatrixfilename){
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
  fasta <- paste(getwd(),"alignment.fasta", sep = "/")
  readfastanogaps <- RemoveGaps(readDNAStringSet(file = fasta))
  write.dna(readfastanogaps, "alignmentnogaps.fasta", format = "fasta")
  nogapsfasta <- paste(getwd(),"alignmentnogaps.fasta", sep = "/")
  ggmsa(nogapsfasta, alignmentstart, alignmentend,  color = "Shapely_NT", font = "DroidSansMono", 
      char_width = 0.5,seq_name = TRUE) + geom_seqlogo(color = "Shapely_NT") +geom_msaBar()
  ggsave("alignment_plot.pdf")
}
