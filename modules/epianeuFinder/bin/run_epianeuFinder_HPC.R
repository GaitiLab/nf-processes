#!/usr/bin/env Rscript

library(argparse)
library(epiAneufinder)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(SeuratObject)
library(tidyr)
library(ggplot2)

parser <- ArgumentParser(description='Process cellranger-arc fragments files with epiAneuFinder using CI and HPC.')

parser$add_argument("-i", "--input", type="character", default=NULL, 
                    help="Input fragment file generated from cellranger-arc",
                    metavar="fragments")
parser$add_argument("-o", "--outdir", type="character", default=NULL, 
                    help="Output directory for epiAneuFinder",
                    metavar="outdir")
parser$add_argument("-b", "--blacklist", type="character", default=NULL, 
                    help=".bed file for blacklist regions of organism genome",
                    metavar="blacklist")
parser$add_argument("-g", "--genome", type="character", default="BSgenome.Hsapiens.UCSC.hg38", 
                    help="Organism reference genome from library BSgenome.",
                    metavar="genome")
parser$add_argument("-w", "--window", type="integer", default=1e5, 
                    help="Window size for epiAneuFinder.",
                    metavar="window")
parser$add_argument("-m", "--minfrags", type="integer", default=20000, 
                    help="Min fragment count for inclusion.",
                    metavar="frags")
parser$add_argument("-n", "--numcores", type="integer", default=2, 
                    help="Number of cores to use for processing",
                    metavar="cores")
parser$add_argument("-t", "--title", type="character", default="epiAneuFinder", 
                    help="Title for the final karyotype plot.",
                    metavar="title")

args <- parser$parse_args()


epiAneufinder(input=args$input,
              outdir=args$outdir,
              blacklist=args$blacklist,
              windowSize=args$window, 
              genome=args$genome,
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo=args$title, 
              ncores=args$numcores,
              minFrags=args$minfrags)

