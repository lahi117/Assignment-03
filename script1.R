# this is the start of the second part of bioinformatics assignment-03

#import libraries
library("seqinr") 
library("R.utils") 
library("rBLAST") 
library("ape")  
library("ORFik")
library("Biostrings")
# download mutblast_functions from marks link
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R") 

# Download the E.colai data sequences
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")

# unzip file
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=TRUE)

# make blast data base
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl", "-parse_seqids")

# download the sequence need to be compared
download.file(url = "https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile = "sample.fa.gz")

#unzip file
R.utils::gunzip("sample.fa.gz",overwrite=TRUE)

#read fasta file
SAMPLE <- read.fasta("sample.fa") 

# we chose 70 sequence
sam <- SAMPLE[[70]]

# check the structure of sam
str(sam)

# get sequence lenth ( lenth of 70 th one)
seqinr::getLength(sam) 

# get gc propotion of the sequence
seqinr::GC(sam)

# run myblastn function
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

write.fasta(sam,names="sam",file.out = "sam.fa")
makeblastdb(file="sam.fa",dbtype = "nucl")

sam_mut <- mutator(myseq=sam,100)
res <- myblastn_tab(myseq = sam_mut, db = "sam.fa")
res
