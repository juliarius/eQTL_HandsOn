## Author: Júlia Rius Bonet
## Date: 22/11/2018
## eQTL Hands-On
#--------------------important--------
#bin--> scripts needed for the pipeline
#input-->unprocessed and processed files
#tmp-->intermediate files
#result--> plots and other result files.

##-------------------- Task 1-----------------

echo "Running command1 to do X"

command1            # This command does X 
                    # Answer to Q1: we need to set the option --opt 
echo -e "\tdone!"

# [...]

#-------------------------- cis eQTL mapping---------------------
#######TASK 1#####
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g # Download the genotype
#######TASK 2#####
PATH=$PATH:$PWD/bin # especifiquem el path
cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt  #to get the geuvadis samples
bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz # we select the common samples, biallelic SNPs and indels, MAF >= 0.05, and exclude the duplicates.
filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz # to select variants with more than 10 individuals per genotype group and compress it so we need it to index it.

#Q1: What do the bcftools options employed mean? 

#view -v snps --> to only view biallelic snps
#indels -m2 -M2 -9 0.05:minor -S --> to exclude the samples (-s) with an MAF minor than 0.05
#-M --> Output sites where REF allele is N="2" 
#-f --> skips the columns that don't have any of the strings listed
#-Ob--> to extract the genotypes
#norm -d --> to normalize indels (-d) the noise observed  
#we stored everything in tmp/genotypes.chr22.vcf.gz
#filter.genotype.py -t 10 -g --> to select those variants with more than 10 individuals per genotype group
#tabix -p --> to index the VCF file
#Q2: How many variants do you get in input/processed/genotypes.chr22.vcf.gz? 
zcat input/processed/genotypes.chr22.vcf.gz | grep -v "#" | wc -l
#74656
#Q3: How many samples do you have before and after subsetting?
bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > tmp/stats.before
less tmp/stats.before
#2504
bcftools stats input/processed/genotypes.chr22.vcf.gz > tmp/stats.after
less tmp/stats.after
#445

##################TASK3#####################3
#Q1: Which version of GENCODE is GEUVADIS using? V12
#Q2: To which genome assembly does this annotation correspond? GRCh37
#Q3: How many protein coding genes are annotated in the last version (v29)?19940
#Q4: Which command do you use to do this?
PATH=$PATH:$PWD/bin
# to download the gencode anotation v12
release=12
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz
mv gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz # we moved it to gencode.annotattion.gtf.gz
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed # we select the 'protein coding' and 'lincRNA' genes and put them in a BED file
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_lenght\"exon_number\"" | gtf2bed.sh > result/task3_5.bed

#Q5:But how to get the TSS positions and the gene lengths from it?
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed #to get the lenght
awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed #to get the tss positions
sed -i "s/^chr//" tmp/gencode.annotation.bed
#Q6: to which BED coordinates would correspond the GTF coordinates chr1 10 20? 
#To chr1 9 20 coz BED is 0 based and GTF is 1 based
#Q7: Why do we need to use tmpfile below?
#it creates a temporary file to editin wich we can work on and once it's closed it gets delated.

##############TASK 4 ################
#We need to join the BED file with the expression file
join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv  # we row orded them by gene ID so we can join the two files
awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv # Subset chr22 (same as the VCF file)
paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed # to recover the column order and to sort rows by chromosome and start position
cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed # to recover the header
#NORMALITZATION
normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed
# Compress and index the final gene expression file
bgzip tmp/genes.chr22.norm.bed
tabix -p bed tmp/genes.chr22.norm.bed.gz
mv tmp/genes.chr22.norm.bed.gz* input/processed

#Q1: Of all genes considered, which have lower expression levels, protein-coding or lincRNA? lincRNA
#Q2: Why do we need gene expression to be normal? So that each gene has a distribution close to standard normal.
#Q3: How would you check that quantile normalization worked? Iwould plot a boxplot and i'd check if the results are more evened out. 
#Q4: and that gene expression of a gene follows a normal distribution?I would make a plot to check if it follows the right distribution, if it's close to the standard.
##############TASK 5 ###############
#before:
check.norm.R -i tmp/genes.chr22.rpkm.bed -o result/plots/check.normbefore.pdf
#After:
check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf

#Q1: What can you see? we can see that the box plots after normalized loked more even,the median is much more similar among the boxes and the quartiles as well, moreover me can see that in the normalized ones they follow a  linear distribution which is what we wanted.
##############################COVARIATES#################3
#############TASK 6###############
#Q1: Which ones would you select? gender,laboratory,population and super population.
head -1 input/unprocessed/geuvadis/geuvadis.metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}' #genotype
head -1 input/unprocessed/1000g/1000g.phase3_metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}' #gene expression
#############TASK 7################
#Q1: What do the parameters employed mean? 
QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression
QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes
#--bed coz it's a bed file --center and --scale can be used to enforce centering and scaling of the phenotype values prior to the PCA
# --maf 0.05 to only consider variant sites with a Minor Allele Frequency (MAF) above 5% --distance 50000 to only consider variant sites separated by at least 50kb
#Q2: Which information do the output files contain?
#result/expression --> The standard deviation,the variance explained, the cumulative variance explained
#result/genotypes --> The same but with the genotype data with ariant sites in linkage equilibrium.
pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf
pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf
#Q3: What can you observe in the plots? we can see that in expression pca the dots are more or less distributed symilarli in boths pca but in genotypes you can distinguish two clusters well separated
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca$
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color pop --out result/plots/genotypes.pca.pop.ge
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color gender --out result/plots/genotypes.pca.gender.pdf ACABAR AMB POP GENDER I SUPER POP
# to determine which fraction of the total variance in the expression data is explained by each potential covariate
join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt
sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf
#Q4. With this information, which covariates seem more relevant to explain the variability in the data? pop and super_pop and lab.
#Q5: Which are the factors that explain more variance? the population and the laboratory
#############TASK 8############### 
# to check how much variance do the first 5 PEER explain in comparison with the known factors
peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf
#Q1: How much variance do they explain? On average is it more or less than the explained by the known factors? araound 55% a bit more.
# generate a covariate file in the format required by QTLtools
join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'
gzip input/processed/covariates.tsv
#############TASK 9###############
#for each phenotype (i.e. gene expression), linear regressions between it and the genotypes of the variants in a window of 1Mb around the TSS.
QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 0.01 --out result/nominals.txt
#Which information contains each field in the ouptut file?(nominals.txt)
#Associations between genotype dosages and phenotype quantifications are measured with linear regressions:
#The columns are:
#-1. The phenotype ID
#-2. The chromosome ID of the phenotype
#-3. The start position of the phenotype
#-4. The end position of the phenotype
#-5. The strand orientation of the phenotype
#-6. The total number of variants tested in cis
#-7. The distance between the phenotype and the tested variant (accounting for strand orientation)
#-8. The ID of the tested variant
#-9. The chromosome ID of the variant
#-10. The start position of the variant
#-11. The end position of the variant
#-12. The nominal P-value of association between the variant and the phenotype
#-13. The corresponding regression slope
#-14. A binary flag equal to 1 is the variant is the top variant in cis
#Q1: Are there pairs genotype-phenotype with exactly the same p-value and effect size (β)?  How is this possible? coz there is snips that are inherited in clusters
pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf
#Have a look at the p-value distribution. Q2: What do you observe?
#You have on the surface a set of well-behaved p-values. The flat distribution along the bottom is all your null p-values, which are uniformly distributed between 0 and 1.
plink --ld rs36084991 rs5746938 --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ld2 
#Q3: Which SNPs did you select? What do you observe?SE ve dónde está el snip en que frequencia cada una de las variaciones i la  frequencia esperada
rs36084991 rs5746938 
#Haplotype     Frequency    Expectation under LE
#   ---------     ---------    --------------------
#         GTA      0.246067                0.060549
#          GA     -0                       0.185518
#         GTT     -0                       0.185518
#          GT      0.753933                0.568414

#   In phase alleles are GTA/GT

##########################################cis eQTL mapping (permutation pass)###
QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --out result/permutations.txt #performs 1000 permutations of the phenotypes and carries out linear regressions on the permuted data analogously to the nominal pass
#############TASK 10###############
for j in $(seq 1 16); do
	echo "cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --chunk $j 16 --out result/permutations_$j.txt"
done | xargs -P4 -n14 QTLtools
cat result/permutations_*.txt > result/permutations.txt; rm result/permutations_*.txt
R
p <- read.table("result/permutations.txt")                                                      # Read input file
pdf("result/plots/pv-correlation.pdf",  paper = 'a4r', width = 9, height = 6)                   # Open PDF device
plot(p[, 18], p[, 19], xlab = "pv (perm)", ylab = "pv (beta)")                                  # Plot p-values
abline(0, 1, col = "red")                                                                       # Add red line 1=1
plot(-log10(p[, 18]), -log10(p[, 19]), xlab = "-log10 pv (perm)", ylab = "-log10 pv (beta)")    # Repeat in -log10 space to check the behaviour of the small p-values.
abline(0, 1, col = "red")
dev.off()                                                                                       # Close device
quit("no")                                                                                      # Exit R
######################################Multiple testing correction###########
############# TASK 11 #############
#Q1: How many significant eQTLs do we find in each case in comparison with the nominal pass? 
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'bonferroni' --alpha 0.05 --out tmp/bonferroni.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'fdr' --alpha 0.05 --out tmp/fdr.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'perm-fdr' --alpha 0.05 --out result/eqtls.tsv 
wc -l result/nominals.txt #37941 
wc -l result/eqtls.tsv #12045
############TASK 12 ##############
eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose #  han salido 9 plots
#####################################eQTL functional analysis #####################
########### TASK 13 ##############
# Download from ftp server
rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl

# Get chr, start, end and feature name in BED format
zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed

# Merge overlapping features of the same type 
# e.g. chr1 100 200 feat1            chr1 100 300 feat1
#      chr1 150 300 feat1     =>     chr1 100 250 feat2
#      chr1 100 250 feat2
for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do 
  bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
done > input/processed/ERB.collapsed.bed

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i "s/^chr//" input/processed/ERB.collapsed.bed
#Perform the enrichment of top eQTLs:

for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do 
  QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat" 
done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt

plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf

#Q1: Which are the top enriched features? Which kind of factors are they? Algunos de los top enriched son:H3K36me3, H3K79me2, H3K4me1, H4K20me1, H3K9ac, H2AZ, H3K4me2,son marcadores de histonas 
#Q2: What does an odds ratio lower than one mean? OR<1 Exposure associated with lower odds of outcome than the expected ones.
########### TASK 14 ##############
sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv # da error utilizo el resultado de la pag.web de ARI

#Q1: Which kind of consequences have they, according to the VEP? In which proportion? 
#intron_variant: 47%
#upstream_gene_variant: 15%
#downstream_gene_variant: 14%
#non_coding_transcript_variant: 11%
#NMD_transcript_variant: 6%
#regulatory_region_variant: 2%
#intergenic_variant: 2%
#non_coding_transcript_exon_variant: 1%
#3_prime_UTR_variant: 1%
#Others

#Q2: How many eQTLs are high impact variants? 23  Which consequences are related to those high impact variants?  frameshift_variant, stop_gained, splice_acceptor_variant,non_coding_transcript_variant,NMD_transcript_variant,splice_donor_variant.
#Q3: Out of all high impact variants, how many of them are falling in acceptor splice sites of protein coding genes? 6
########## TASK 15 ##############  GO enrichment##########
cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt # to generate a list of genes
awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt #We will use as background all the genes (PC and lincRNA) in chr22

#Q1: In which biological processes are your eGenes enriched? response to lipopolysaccharide and response to molecule of bacterial origin  Which molecular functions and components correspond to those processes? Ras guanyl-nucleotide exchange factor activity 
########### TASK 16 ##################4. eQTL and GWAS co-localization#######
# Generate input files for QTLtools rtc
grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input
cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait

# Download the file 'hotspots_b37_hg19.bed' from QTLtools website
wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed

# Run RTC
QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt

#Q1: How many pairs of variants have a RTC value above 0.9? 39 
awk '{if($20>0.9) print $20}' result/rtc.txt | wc -l
#Q2: For each pair, we have a GWAS hit and an eQTL. Find one example so that the gene to which the eQTL is associated is relevant for the trait/disease to which the GWAS variant is associated. Explore the literature and the biological databases that you know to gather more information. 
awk '$20=="1"' result/rtc.txt  # rs909685 rs909685 ENSG00000100321.10 ENSG00000100321.10 # involved in short-term and long-term synaptic plasticity 

#Q3: Which consequences, according to the variant effect predictor, do these co-localized eQTL variants have? Intergenic region

########################################eQTL fine-mapping CAVIAR ###############  
# Generate the ID/Z-scores input file. Select your favourite gene (e.g. gene=ENS00000000000.0).
# Set k (number of variants) to 50
gene=ENSG00000099968.13
compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z

# Generate the LD matrix 
plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene
#############  TASK 17 ##########
CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene #Obtain a credible set of causal variants with probability ρ=0.95.
#Q1: How many variants are there in the credible (ρ=0.95) set? 4 For each of these variants, which is the probability to be causal? 
awk '$1=="rs8137591"' result/ENSG00000198951.6_post
awk '$1=="rs8141347"' result/ENSG00000198951.6_post  
awk '$1=="rs1978967"' result/ENSG00000198951.6_post  
awk '$1=="rs7290691"' result/ENSG00000198951.6_post    
#Q2: Which are the p-values and effect sizes of these variants? 
#rs8137591
awk '{print $1, $12, $13}' <(awk '$8=="rs8137591"' <(awk '$1=="ENSG00000099968.13"' result/nominals.txt))
#rs8141347
awk '{print $1, $12, $13}' <(awk '$8=="rs8141347"' <(awk '$1=="ENSG00000099968.13"' result/nominals.txt))
#rs1978967
awk '{print $1, $12, $13}' <(awk '$8=="rs1978967"' <(awk '$1=="ENSG00000099968.13"' result/nominals.txt))
#rs7290691
awk '{print $1, $12, $13}' <(awk '$8=="rs7290691"' <(awk '$1=="ENSG00000099968.13"' result/nominals.txt))
#Q3: How are they in comparison to the p-values and effect sizes of other variants tested for the same gene?  smaller
awk '{print $1, $12, $13}' <(awk '$1=="ENSG00000099968.13"' result/nominals.txt) 
#Q4: Which consequences, according to the variant effect predictor, do these variants have? intron_variant: 78% non_coding_transcript_variant: 13% NMD_transcript_variant: 9%
#http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=fO6rGseZF6zndlLq-4757169
############## TASK 18 ##########
# Generate a LocusZoom plot. Use as SNP any of the colocalized or fine-mapped variants. 
# Define the gene corresponding to the co-localized or fine-mapped variants of interest
gene=ENSG00000099968.13 
cat <(echo "MarkerName P.value") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$geneENSG00000099968.13

############## TASK 19 ##########
# Add, commit and push these changes to a GitHub repository with the same name.
