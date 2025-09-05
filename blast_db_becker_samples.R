library(metagMisc)
library(nationalparkcolors)
library(phyloseq)
library(microbiome)
library(CoDaSeq)
library(zCompositions)
library(cowplot)
library(googleway)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(vegan)
library(tidyverse)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(corncob)
library(DESeq2)
library(seqinr)
library(rBLAST)
library(tidyverse)

#Load Libraries
otu <- read.table("silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)
samples$Region[samples$Region == 'Curacao'] <- 'Curaçao'
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps

#Remove problematic samples
#Remove problematic samples from ps
ps<- remove_samples(c("17818", "17797"), ps)
summarize_phyloseq(ps)


#Remove samples that have less than 500 reads
ps <- prune_samples(sample_sums(ps) >= 500, ps)
ps #7559 taxa across 132 samples
#Prune taxa no longer in samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps) 
ps #7451 taxa in 132 samples

summarize_phyloseq(ps)
ps_ra<- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

#Read in Becker Bioindicator FASTA
bioindics<- read.fasta(file = "Bioindicator_ASVs_Becker_etal_2020.fasta", seqtype="DNA",
                       as.string = TRUE, forceDNAtolower = FALSE)
bioindics1<- unlist(bioindics)

#Take the rownames of the OTU Table and turn them into a FASTA file
sequences<- colnames(otu)
names(sequences)<- paste('ASV',1:7559, sep="")
write.fasta(as.list(sequences), names = names(sequences), file.out = "dendro_ASVs.fasta", open = "w")

#Read in the dendro FASTA file
dendro_seqs<- read.fasta(file = "dendro_ASVs.fasta", seqtype = "DNA", as.string = TRUE,
                         forceDNAtolower = FALSE)
dendro_seqs1<- unlist(dendro_seqs)



#Blast Becker Sequences against Dendrogyra Sequences
# Created the database on the command line
#makeblastdb -in "${DB_FILE}" -dbtype "nucl" -out "dendro_db" -logfile "makedb.log"

blast_out<- system2(command = "/usr/local/ncbi/blast/bin/blastn", 
        args = c("-db dendro_db -query Bioindicator_ASVs_Becker_etal_2020.fasta -outfmt 6 -evalue 10e-6 -ungapped"),
        wait = TRUE,
        stdout = TRUE)

colnames <- c("qseqid",
              "sseqid",
              "pident",
              "length",
              "mismatch",
              "gapopen",
              "qstart",
              "qend",
              "sstart",
              "send",
              "evalue",
              "bitscore")

tidy_blast <- blast_out %>%
  as_tibble %>% 
  separate(col = value, 
           into = colnames,
           sep = "\t",
           convert = TRUE)

#Get the ASVs that had 100% identity
tidy_blast_100<- tidy_blast[tidy_blast$pident==100,]
tidy_blast_100<- tidy_blast_100[tidy_blast_100$length==126,]

#Get sequences for ASVs with 100% identity
iden_asvs<- tidy_blast_100$sseqid
iden_seqs<- dendro_seqs1[names(dendro_seqs1) %in% iden_asvs]

#Subset phyloseq object by the sequences
ps_ra_subset<- prune_taxa(iden_seqs, ps_ra)

#Create dataframe 
otu.asvs<- as(otu_table(ps_ra_subset), "matrix")
taxon.asvs<- as(tax_table(ps_ra_subset), "matrix")
meta.asvs<- as(sample_data(ps_ra_subset), "matrix")
otu.asvs<- as.data.frame(otu.asvs)
otu.asvs<- rownames_to_column(otu.asvs, var="Sample")

meta.asvs<- as.data.frame(meta.asvs)
meta.asvs<- rownames_to_column(meta.asvs, var="Sample")
otu.asvs.meta<- merge(meta.asvs, otu.asvs, "Sample")

otu.asvs.long<- melt(otu.asvs.meta, id.vars=c("Sample", "Reef", "Region", "Collection.Date", "Season"), variable.name = "ASV", value.name="Proportion" )
otu.asvs.long.nonzero<- otu.asvs.long[otu.asvs.long$Proportion>0,]

taxon.asvs<- as.data.frame(taxon.asvs)
genera_labels<- c("Vibrio ASV1", "Vibrio ASV2", "Shimia", "Vibrio ASV3",
                  "Fusibacter", "Thalassobius", "Fusibacter", "Arcobacteraceae",
                  "Halarcobacter", "Algicola", "Vibrio ASV4", "Vibrio ASV5",
                  "Halodesulfovibrio")
genera_labeller <- function(variable,value){
  return(genera_labels[value])
}

region_palette<-c("Belize"="#003459", "Curaçao"="#75dddd")

pdf("Bioindicator_ASVs_Relative_Abundance.pdf",bg = "white",width=8.5, height=11)
p1<- ggplot(otu.asvs.long.nonzero, aes(x=Region, y=as.numeric(Proportion), fill=Region))+
  geom_boxplot()+
  scale_fill_manual(values = region_palette)+
  theme_classic() +
  facet_wrap(~ASV, nrow=5, ncol=3, labeller = genera_labeller)+
  theme(strip.text = element_text(face="bold", size = 16))+
  theme(strip.background = element_rect( fill="white"))+
  theme(legend.title=element_blank()) +
  theme(text=element_text(size=14))+
  theme(axis.text.y = element_text(color="black", size = 16))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_text(size = 20, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  ylab('Relative Abundance')
p1
dev.off()