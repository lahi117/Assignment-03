# this is the start of the second part of bioinformatics assignment-03
library("seqinr") 
library("R.utils") 
library("rBLAST") 
library("ape") 
library("ORFik") 
library("Biostrings") 
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=TRUE)
#library("rBLAST")
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl", "-parse_seqids")
download.file(url = "https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile = "sample.fa.gz")
R.utils::gunzip("sample.fa.gz",overwrite=TRUE)
#library("rBLAST")
#makeblastdb("sample.fa",dbtype="nucl","-parse_seqids")
SAMPLE <- read.fasta("sample.fa") 
sam <- SAMPLE[[70]]
str(sam)
seqinr::getLength(sam)
seqinr::GC(sam)
myblastn_tab
res <- myblastn_tab(myseq = SAMPLE, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(res)

head(res)
hits <- as.character(res$sseqid[1:3])
hits
db <- read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")

tophit <- db[which(names(db) %in% hits[1])] # extract the names of the top hit
tophit #[1:50]

seqinr::write.fasta(tophit,names=names(tophit),file.out = "tophit.fa")
makeblastdb("tophit.fa",dbtype="nucl", "-parse_seqids")
res <- myblastn(myseq = SAMPLE, db = "tophit.fa")
cat(res,fill=TRUE)
mutator

sam_mut <- mutator(myseq=sam,100)
sam_mut_ <- DNAString(c2s(sam_mut))
sam_ <- DNAString(c2s(sam))
aln <- Biostrings::pairwiseAlignment(sam_,sam_mut_)
pid(aln)
nmismatch(aln)

