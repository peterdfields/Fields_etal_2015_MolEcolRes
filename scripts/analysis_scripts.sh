### Commands arising from Molecular Ecology Resources Note,Transcriptome resources for two highly divergent Silene latifolia populations, by Peter D. Fields, Laura A. Weingartner, and Lynda Delph.

### much of the following manipulation of files relies on a tutorial provided at http://khmer-protocols.readthedocs.org/en/v0.8.4/mrnaseq/1-quality.html, and relies on tools provided by Brown et al. (2013), though digital normalization was not used in the submitted assembly.

# Verification of GENEWIZ removal of adapter sequences with Trimmomatic
java -jar trimmomatic-0.32.jar PE *_R1.fastq.gz *_R2.fastq.gz s1_pe s1_se s2_pe s2_se ILLUMINACLIP:/usr/local/share/adapters/TruSeq3-PE.fa:2:30:10

# interleave paireded reads
khmer/scripts/interleave-reads.py s1_pe s2_pe | gzip -9c > sample_name.pe.fq.gz

# combine single end reads
cat s1_se s2_se | gzip -9c > sample_name.se.fq.gz

# use fastqc to verify removal of adapters
fastqc sample_name.pe.fq.gz
fastqc sample_name.se.fq.gz
## fastqc .html reports are assessed using Firefox browser

# quality trim reads w/fastx-toolkit
gunzip -c sample_name.pe.fq.gz | fastq_quality_filter -Q33 -q 30 -p 50 | gzip -9c > sample_name.pe.qc.fq.gz
gunzip -c sample_name.se.fq.gz | fastq_quality_filter -Q33 -q 30 -p 50 | gzip -9c > sample_name.se.qc.fq.gz

# extract still paired reads
for i in *.pe.qc.fq.gz
do
   khmer/scripts/extract-paired-reads.py $i
done

# rename files created by extract-paired-reads.py script
# paired
for i in *.pe.qc.fq.gz.pe
do
   newfile="$(basename $i .pe.qc.fq.gz.pe).pe.qc.fq"
   mv $i $newfile
   gzip $newfile
done

# single
for i in *.pe.qc.fq.gz.se
do
  otherfile="$(basename $i .pe.qc.fq.gz.se).se.qc.fq.gz"
  gunzip -c $otherfile > combine
  cat $i >> combine
  gzip -c combine > $otherfile
  rm $i combine
done

# split paired reads
khmer/scripts/split-paired-reads.py sample_name.pe.qc.fq.gz

# concantenate outputs from split-paired-reads and include single end reads into one of the sets of reads
cat *.1 > left.fq
cat *.2 > right.fq

gunzip -c *.se.qc.fq.gz >> left.fq

# reads from left.fq and right.fq were then used as inputs for blacklight. See blacklight_script.txt

#### Identify SNP polymorphisms
# index meta assembly reference with bowtie2
bowtie2-build Trinity.fasta meta_assembly

# align reads to assembled transcriptome, creating sam file
bowtie2 -x meta_assembly -1 sample_name_left.fq -2 sample_name_right.fq -S sample_name.sam

# compress .sam to save space
gzip sample_name.sam

# use samtools to generate .bam file and remove low quality reads
samtools view -q 25 -bS sample_name.sam.gz > sample_name.bam

# use samtools to sort .bam file
samtools sort sample_name.bam sample_name.sorted

# use samtools to index assembly reference
samtools index Trinity.fasta

# use samtools to create pileup file and bcftools to generate vcf file
samtools mpileup -uf Trinity.fasta *.sorted.bam | bcftools view -bvcg - > all_samples.raw.bcf
bcftools view all_samples.raw.bcf | vcfutils.pl varFilter -D100 > silene_transcriptomes.vcf 

# use vcftools to create .vcf with only high quality SNPs
vcftools --vcf silene_transcriptomes.vcf --recode --recode-INFO-all --out silene_transcriptome_highquality.vcf --minQ 20

### Annotations with Trinotate, updated to most recent version of Trinotate
# create pep file
TransDecoder.LongOrfs -t Trinity.fasta

# search transcripts Trinotate based databases
blastx -query Trinity.fasta -db uniprot_sprot.trinotate.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

# search against .pep file
blastp -query transdecoder.pep -db uniprot_sprot.trinotate.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

# run hmmer
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm transdecoder.pep > pfam.log

# run signalp
signalp -f short -n signalp.out transdecoder.pep

# run tmhmm
tmhmm --short < transdecoder.pep > tmhmm.out

# run rnnammer
RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer ~/bioinformatics/rnammer_v1.2/rnammer

# retrieve pre-generated trinotate sql database
wget "ftp://ftp.broadinstitute.org/pub/Trinity/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz" -O Trinotate.sqlite.gz

gunzip Trinotate.sqlite.gz

# load files into database
get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep transdecoder.pep
# load protein hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
# load transcript hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
# load pfam hits
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
# load transmembrane domain hits
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.ou
# load signal peptide predictions
Trinotate Trinotate.sqlite LOAD_signalp signalp.out

# output Trinotate report
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls


Literature cited:
Brown, C. Titus; Scott, Camille; Crusoe, Michael; Sheneman, Leigh;
Rosenthal, Josh; Adina Howe (2013): khmer-protocols documentation.
figshare. http://dx.doi.org/10.6084/m9.figshare.878460