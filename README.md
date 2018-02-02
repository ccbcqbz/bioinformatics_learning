# Goal

- Learn Bioinformatics basics and understand the computational/storage
requirement of the problems in the area. 
- Focus on *high performance computing, big data, and machine learning*.

# Background

Bioinformatics is one of the areas that have been rapidly advancing.
For example, the cost of sequencing a whole human genome is getting
cheaper and cheaper, and
[it is outpacing outpace Moore's law](https://www.genome.gov/sequencingcosts/).
[An article](https://techcrunch.com/2017/01/10/illumina-wants-to-sequence-your-whole-genome-for-100/)
said we will be able to sequence whole human genome for $100 in near
future. We expect that the rapid progress on cost reduction will allow
us to produce a large amount of genome data and can invent
technologies that will fundamentally change our life such as cancer
detection or drug discovery/personalization.

Storing and processing a high volume of genome data is very
challenging. For example, the human genome contains approximately 3
billions base pairs. Also, genome data processing (e.g., assembly,
sequence alignment) is a compute intensive job that can take easily
more than one day to complete. 

Cloud has been a natural place for such massive data processing as its
cost advantage and rich services offered (e.g.,
[Genomics in the Cloud](https://aws.amazon.com/health/genomics/),
[Google Genomics](https://cloud.google.com/genomics/)).
(Supercomputers and in-house datacenters are other options, but we
don't discuss here.) Companies like
[DNAstack](https://www.dnastack.com) or
[Seven Bridges Genomics](https://www.sevenbridges.com/) have been
building platforms for genomics data analysis and sharing on top of
cloud.
    
# Learning Material

This section summarizes material (e.g., course, papers, software) that
would be useful for learning Bioinformatics. (*Disclaimer*: I didn't go
thorough all of them mentioned here.)

## Overview

Coursera courses:

- [Bioinformatics Specialization](https://www.coursera.org/specializations/bioinformatics).
  This specialization consists of seven courses and covers various
  topics such as motif search, genome assembly with De Bruijn graphs,
  Burrows-Wheeler read mapping, profile hidden Hidden Markov Models
  for sequence alignment. The final course uses the BaseSpace cloud
  platform developed by Illumina.
- [Genomic Data Science Specialization](https://www.coursera.org/specializations/genomic-data-science#courses)

Text book, etc.:
- [Phillip Compeau, Pavel Pevzner, "Bioinformatics Algorithms - An Active Learning Approach"](http://bioinformaticsalgorithms.com/). This is a textbook used by Bioinformatics Specialization.
- [Rosalin: A platform for learning bioinformatics and programming through problem solving](http://rosalind.info/problems/locations/). Most (?) of the programming problems in this site also come from Bioinformatics Specialization.
- [String Algorithms and Algrorithms in Computational Biology - Gusfield](http://web.cs.ucdavis.edu/~gusfield/cs224f11/)

Software list:
- [Awesome Bioinformatics](https://github.com/danielecook/Awesome-Bioinformatics). A curated list of Bioinformatics software, resources, and libraries.

## Sequence Alignment

DNA/RNA sequence and alignment is one of the popular problems in
Bioinformatics. Burrows–Wheeler transform is one of the classic
algorithms in the area, and various efforts have been made to run a
sequence algorithm in parallel/distributed environments. In particular, the SIGMOD papers 
listed here are very good for me to understand the problem and
how recent technologies on parallel/distributed computing is used. 

Please see a later section to know the details of the sequence alignment problem.

Papers:

- [A. Roy et al., "Massively parallel processing of whole genome sequence data: an in-depth performance study". SIGMOD, 2017](https://dl.acm.org/citation.cfm?id=3064048). This paper
  describes a system called [GESALL](http://gesall.cs.umass.edu/). It
  aims to run existing genomic data analysis programs (e.g.,
  [GATK](https://software.broadinstitute.org/gatk/)) without rewriting.
- [C. Miller et al., "bam.iobio: a web-based, real-time, sequence alignment file inspector". Nature Methods, 2014](http://www.nature.com/nmeth/journal/v11/n12/full/nmeth.3174.html)
- [F. Nothaft, et al., "Rethinking data-intensive science using scalable analytics systems". SIGMOD, 2015](https://amplab.cs.berkeley.edu/wp-content/uploads/2015/03/adam.pdf).
  This paper describes a system called
  [ADAM](https://github.com/bigdatagenomics/adam). It uses Apache
  Spark and Parquet to run (part of) sequencing and alignment in parallel.
- [A. Dobin, et al., "STAR: ultrafast universal RNA-seq aligner". Bioinformatics, 2013](https://academic.oup.com/bioinformatics/article/29/1/15/272537/STAR-ultrafast-universal-RNA-seq-aligner)
- [H. Li, A Statistical Framework for SNP Calling, Mutation Discovery, Association Mapping and Population Genetical Parameter Estimation from Sequencing Data". Bioinformatics 27.21 (2011)](https://www.ncbi.nlm.nih.gov/pubmed/21903627)
- [Y. Li, et al., "WHAM: a high-throughput sequence alignment method". SIGMOD, 2011](http://pages.cs.wisc.edu/~jignesh/publ/wham.pdf)
- [M. Zaharia, et al., "Faster and more accurate sequence alignment with SNAP". arXiv preprint, 2011](http://arxiv.org/abs/1111.5572)
- [H. Li and R. Durbin, "Fast and accurate long-read alignment with Burrows–Wheeler transform". Bioinformatics, 2010](https://academic.oup.com/bioinformatics/article/26/5/589/211735/Fast-and-accurate-long-read-alignment-with-Burrows)
- [H. Licorresponding and N. Homer, "A survey of sequence alignment algorithms for next-generation sequencing", Briefings in Bioinformatics, 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2943993/)
  
Videos:
- [J. Bloom and T. Poterba, "Scaling Genetic Data Analysis with Apache Spark". Spark Summit, 2017](https://www.youtube.com/watch?v=pyeQusIN5Ao)

Software:
- [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/): A
  software package for mapping low-divergent sequences against a large
  reference genome.
- [SAMtools](http://www.htslib.org/): A suite of programs for
  interacting with high-throughput sequencing data.
- [Genome Analytic Toolkit (GATK)](https://software.broadinstitute.org/gatk/):
  A toolkit offering a wide variety of tools with a primary focus on
  variant discovery and genotyping. [This article]( https://software.broadinstitute.org/gatk/blog?id=9645) and [this article](https://blog.cloudera.com/blog/2016/04/genome-analysis-toolkit-now-using-apache-spark-for-data-processing/) explain the vision of GATK4 (e.g., Spark-based like ADAM).
- [GESALL](http://gesall.cs.umass.edu/): A system for end-to-end processing of the genomic data.
- [ADAM](https://github.com/bigdatagenomics/adam): A genomics analysis
  platform with specialized file formats built using Apache Avro,
  Apache Spark and Parquet.
- [SNAP](http://snap.cs.berkeley.edu/): A new sequence aligner that is
  3-20x faster and just as accurate as existing tools.
- [STAR](https://giathub.com/alexdobin/STAR/): An ultrafast universal RNA-seq aligner.
- [Bam.iobio](http://bam.iobio.io/)
- [Hail](https://hail.is/): An open-source, scalable framework for exploring and analyzing genomic data. Starting from genetic data in VCF, BGEN or PLINK format.
- [GenomicsDB](https://github.com/Intel-HLS/GenomicsDB/wiki): A database built on top of the [TileDB](http://istc-bigdata.org/tiledb/index.html), which is a system for efficiently storing, querying and accessing sparse matrix/array data. 

### Data Science and Machine Learning

Machine learning is one of the areas that we see very rapid progress,
and people are trying to apply the technology to solve Bioinformatics
problems (e.g., why not try LSTM instead of HMM as people do for
speech recognition). For example, one of the goals would be to predict phenotypes, such as disease risks, from a genotype. 

Papers:
- [T. Ching, et al., "Opportunities and obstacles for deep learning in biology and medicine". 2017](https://greenelab.github.io/deep-review/) This paper summarizes the current status of deep learning research and its challenges. The following three research areas are discussed with various example research work: (1) disease and patient categorization, (2) fundamental biological study, (3) treatment of patients. 
- [M. Leung, et al., "Machine learning in genomic medicine: a review of computational problems and data sets". Proceedings of the IEEE, 2016](http://ieeexplore.ieee.org/document/7347331/). This paper summarizes how machine learning can solve problems in genomic medicine and proposes a cell variables approach where models are trained to predict how genotype influences “cell variables” such as concentrations of proteins.
- [R. Poplin, et al., "Creating a Universal SNP and Small Indel Variant Caller with Deep Neural Networks". 2016](https://doi.org/10.1101/092890). DeepVariant paper
- [R. Irizarry, et al., "Biomedical data science". Course notes, EdX PH525.1x 2015](http://genomicsclass.github.io/book/)
- [B. Yoon, "Hidden Markov Models and their applications in biological sequence analysis". Curr Genomics, 2009](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2766791/)
- [P. Larrañaga, et al., "Machine learning in bioinformatics". Brief Bioinform, 2006](https://academic.oup.com/bib/article/doi/10.1093/bib/bbk007/264025/Machine-learning-in-bioinformatics)
- https://github.com/hussius/deeplearning-biology: A list of implementations of deep learning methods to biology

Programming contests:
- [PrecisionFDA Truth Challenge](https://precision.fda.gov/challenges/truth/results)
- Merck Molecular Activity Challenge
- [NF2 Tumor Genomics Hackathon](https://sv.ai/hackathon/)

Software:
- [MEME Suite](http://meme-suite.org/): Motif-based sequence analysis tool
- [Bioconductor](https://www.bioconductor.org/): Tools for the analysis and comprehension of high-throughput genomic data.
- [DragoNN](http://kundajelab.github.io/dragonn/): A toolkit to teach and learn about deep learning for genomics.

## Data Compression

Papers:
- [SAM/BAM format specification](http://samtools.github.io/hts-specs/SAMv1.pdf)
- [CRAM format specification version 3.0](https://samtools.github.io/hts-specs/CRAMv3.pdf)
- [Asymmetric Numeral Systems](https://en.wikipedia.org/wiki/Asymmetric_Numeral_Systems)

Software:
- [CRAM](http://www.ebi.ac.uk/ena/software/cram-toolkit): A file format
  and toolkit for highly efficient and tunable reference-based
  compression of sequence data.

## Cloud Computing

Papers:
- [Muir P et al., "The real cost of sequencing: scaling computation to keep pace with data generation". Genome Biology 17, 2016](https://doi.org/10.1186/s13059-016-0917-0)
- [David H et al., "A Million Cancer Genome Warehouse", UCB/EECS-2012-211, 2012](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2012/EECS-2012-211.html). This technical paper summarizes the technical requirements and challenges of genome processing very well.

## Genome Databases 

- [NCBI BLAST Database](https://blast.ncbi.nlm.nih.gov/smartblast/?LINK_LOC=BlastHomeLink)
- [Pfam](http://pfam.xfam.org/): A database of HMMs and multiple alignments representing protein families.
- UniProt: A comprehensive resource for protein sequence and annotation data.
- WormBase: An authoritative nonredundant database of nematode predicted protein sequences.

## Workflow Management

- [Illumina BaseSpace](https://basespace.illumina.com/apps/): A
  cloud-based genomics computing environment for next-generation
  sequencing (NGS) data management and analysis.
- [Galaxy](https://usegalaxy.org/): A web-based platform for data intensive biomedical research.

## GitHub

- https://github.com/lh3/bwa (Burrows-Wheeler Aligner)
- https://github.com/samtools (samtools, htslib, bcftools)
- https://github.com/broadinstitute (gatk, picard)
- https://github.com/enasequence/cramtools
- https://github.com/bigdatagenomics (adam, snappea)
- https://github.com/BD2KGenomics
- https://github.com/opencb
- https://github.com/kundajelab/dragonn
- https://github.com/googlegenomics

## Companies, Research Institutes, and Projects
   
- [BROAD Institute](https://www.broadinstitute.org/)
- [Illumina](https://www.illumina.com/)
- [Deep Genomics](https://www.deepgenomics.com/)
- [Grail](https://grail.com/)
- [Seven Bridges Genomics](https://www.sevenbridges.com/)
- [DNAstack](https://www.dnastack.com/#/)
- [DNAnexus](https://www.dnanexus.com/)
- [Color Genomics](https://www.color.com/)
- [Sentieon](https://www.sentieon.com/)
- [GENOOX](https://genoox.com/)
- [OpenCB](http://www.opencb.org/) OpenCB provides advanced open-source software for the analysis of high-throughput genomic data.
- [CANDLE - Exascale CANcer Distributed Learning Environment](https://cbiit.nci.nih.gov/ncip/hpc/candle)
- [Computational Genomics Research Group](http://compbio.berkeley.edu/)
- [National Cancer Institue Genomic Data Commons](https://gdc.cancer.gov/) 
- [iobio](http://iobio.io/): A team of developers in the Marth lab at the Center for Genetic Discovery at the University of Utah.

## Conferences, Workshops, and Journals

- Bioinformatics
- SIGMOD 
- [IEEE/ACM Transactions on Computational Biology and Bioinformatics](http://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=8857)

# Potential Areas for Improvement

This section discusses potential areas where we might be able to contribute to Bioinformatics and its practical applications. 

- *Modern engineering to existing research tools*: One problem we've been
  seeing is that some of the existing tools are written by researchers,
  not professional fulltime software engineers. The quality of overall code
  can be improved and that will allow people to use or extend the tools
  more reliably. 
- *Test dataset, benchmarking, quality evaluation*: One obstacle
  preventing people from making changes to existing tools or trying
  new tools is that there is no easy way to assess the quality of
  tools and compare them. Suppose that we make performance
  optimization to to BWE and that changed an output sequence slightly.
  How can we tell how it affects the final output of a workflow
  consisting of many pipelines?
- *Cloud*: TBD. How can we run existing tools efficiently in Cloud environment? Is there
  any problem specific to Bioinformatic? For example, can we simply upload large genome
  data to Cloud? Is it too slow? Is there any preprocessing needed?
- *Machine learning*: TBD.
- A way to present problems to computer scientists and engineers without
  requiring too much biology knowledge.
- Sophisticated web service for interactive workflow. 

The following areas are interesting, but there might not be right problems to be solved.
  
- *Parallel processing framework and batch task management*: MapReduce (Hadoop)
  and Spark might be just sufficient here. The scheduling problem is
  also a classic one in high performance computing (e.g., gang
  scheduling) if we consider only batch jobs. 
- *Storage systems*: Most data we manage is immutable. Might need to
  be compatible with file formats used by existing systems. We already
  have choices like HDFS, Bigtable, and Parquet. There is an issue
  with how to store genome data efficiently on these storage services,
  but the problem would probably be better solved by changing tools
  (not by inventing a new storage system).
- Common API and tool integration

## Project Ideas

- Cloud cost analysis: There are various solutions we can choose to
  satisfy our requirements, and we need to evaluate the solutions
  based on certain criteria. Cost in one of the important ones to consider.
- Code quality improvement to [BWA](https://github.com/openpcb/bwa).
  BWA is written in C with no unit test code. Related work:
  [this code](https://github.com/ytchen0323/cloud-scale-bwamem) and
  [this discussion](https://github.com/bigdatagenomics/adam/issues/1311).
- Sequence algorithm improvement. My naive question is here is whether we can 
  use more than one reference genome to improve the sequencing accuracy. For example,
  if the sequence speed is high enough, we should be able to sequence a given genome against 
  multiple reference genomes and can get a better result? One question is how important for us
  to improve sequencing accuracy. Will that affect the accuracy of the final result like variant calling?

# Miscellaneous Notes

## DNA, RNA, and Protein

- The Central Dogma of Molecular Biology: "DNA makes RNA makes protein.”
- DNA: Sequence of {A, C, G, T}
- RNA: Sequence of {A, C, G, U}
- Protein: Sequence of amino acids 
- DNA-to-RNA transcription: Transform a DNA string into an RNA string by
  replacing all occurrences of T with U.
- RNA-to-protein translation: Partition an RNA strand into
  non-overlapping 3-mers called codons. Then, convert each codon into
  one of 20 amino acids via the genetic code.
- Chromesome: A set of genes(?)
- Gene: A set of DNA sequences(?)

### Paired Read 

A paired-end (PE) sequence is a common practice in next-generation
sequencing. Given a test genome, lab processes are used to break the
DNA into small fragments. Then, a PE sequencer "reads" each double-stranded 
DNA fragment from both sides, and aligns the forward and reverse reads as 
*paired reads, with the distance between them known for a specific sequencer
(cited from GESALL paper).

```
Forward         Reverse
Read               Read 
o------>       <------o
=======================
Reference
```

Numbers:
- Short read: 50-250 bases in length
- A single human genome has approximately 3 billion bases.
- A single human genome sequenced at 65x coverage will produce
  approximately 1.4 billions reads of 150 base length (cited from the ADAM
  paper)
- For a human genome sample a sequencer can produce one billion short
  reads of 200-1000 bases each, totaling 0.5–1 TB of data, within three
  days and at a cost of a few thousand dollars (cited from the GESALL
  paper).
- About 20,000 genes in the human genome 

### Transcription

- The DNA sequence is first transcribed into precursor mRNA (pre-mRNA).
- The pre-mRNA is further processed to make a messenger RNA (mRNA).
  Various modifications take place during the RNA processing (e.g., splicing, polyadenylation).

Splicing:
- The DNA sequence contains both exons and introns.
- Splicing removes the introns from the pre-mRNA and connects 
  the exons together.
- In the standard model, splicing removes introns and retains all
  exons, but most genes may be spliced in different ways, so that exons
  are sometimes removed and/or introns are retained, which increases the
  variety of proteins.
- Splicing may occur after transcription is complete, but it frequently
  occurs concurrently with transcription.
  
Polyadenylation:
- A pervasive mechanism responsible for regulating mRNA function, stability, localization, and translation efficiency. 
- Append a sequence of adenine bases to the end of the mRNA.

## DNA Sequence / Alignment

For each read, we find the position in the genome that the read is
most likely to have come from. As an exact search is too expensive,
there has been an extensive amount of research that has focused on
indexing strategies for improving alignment performance. This process
is parallel per sequenced read (cited from the ADAM paper).

We can also define the problem as follow: Given an error model
specifying *s* substitutions, *i* insertions, and *d* deletions, a
valid alignment for a given query sequence *Q* is all the subsequences
in the reference sequence *R* that can be aligned by applying up to
*s* substitutions, *i* insertions, and *d* deletions between *Q* and
each matched subsequence in *R* (cited from the WHAM paper).

There are two approaches to the problem. One is hash-based
seed-and-extend methods. The other is prefix trie methods based on the
Burrows-Wheeler transform (BWT). An optimal algorithm would be
different depending on input types (e.g., short reads v.s. long reads)
and alignment goal (e.g., local v.s. global).

The above description is for the primary analysis, and there are steps
for second analysis and tertiary analysis such as preprocessing,
variant calling, and filtering. A variant call is a conclusion that
there is a nucleotide difference vs. some reference at a given
position in an individual genome or transcriptome. It is usually
accompanied by an estimate of variant frequency and some measure of
confidence
([reference](https://www.bioconductor.org/help/course-materials/2014/CSAMA2014/3_Wednesday/lectures/VariantCallingLecture.pdf)).

Here is a citation from DeepVariant paper: For example, the widely-used GATK uses logistic regression to model base errors, hidden Markov models to compute read likelihoods, and naive Bayes classification to identify variants, which are then filtered to remove likely false positives using a Gaussian mixture model with hand-crafted features capturing common error modes. These techniques allow the GATK to achieve high but still imperfect accuracy on the Illumina sequencing platform3,4. Generalizing these models to other sequencing technologies has proven difficult due to the need for manual retuning or extending these statistical models (see e.g. Ion Torrent8,9), a major problem in an area with such rapid technological progress.

## RNA Sequence (RNA-seq) / Alignment

Two key tasks make these analyses computationally intensive. The
first task is an accurate alignment of reads that contain mismatches,
insertions and deletions caused by genomic variations and sequencing
errors. The second task involves mapping sequences derived from
non-contiguous genomic regions comprising spliced sequence modules
that are joined together to form spliced RNAs. Although the first task
is shared with DNA resequencing efforts, the second task is specific
and crucial to the RNA-seq, as it provides the connectivity
information needed to reconstruct the full extent of spliced RNA
molecules. These alignment challenges are further compounded by the
presence of multiple copies of identical or related genomic sequences
that are themselves transcribed, making precise mapping difficult (cited from the STAR paper).

## Finding Protein-binding DNA/RNA

- ChIP-seq: common experimental protocol for DNA-binding proteins
- RIP- and CLIP-seq: RNA-binding proteins
- Motif search
- MEME
- Position Frequency Matrix (PFM) and beyond

### Regulatory Proteins and Transcription Factors

(mostly cited from Bioinformatics Specialization course)

It turns out that every plant cell keeps track of day and night
independently of other cells, and that just three plant genes, called
LCY, CCA1, and TOC1, are the clock’s master timekeepers. Such
regulatory genes, and the *regulatory proteins* that they encode, are
often controlled by external factors (e.g., nutrient availability or
sunlight) in order to allow organisms to adjust their gene expression.

For example, regulatory proteins controlling the circadian clock in
plants coordinate circadian activity as follows. TOC1 promotes the
expression of LCY and CCA1, whereas LCY and CCA1 repress the
expression of TOC1, resulting in a *negative feedback loop*. In the
morning, sunlight activates the transcription of LCY and CCA1,
triggering the repression of TOC1 transcription. As light diminishes,
so does the production of LCY and CCA1, which in turn do not repress
TOC1 any more. Transcription of TOC1 peaks at night and starts
promoting the transcription of LCY and CCA1, which in turn repress the
transcription of TOC1, and the cycle begins again.

LCY, CCA1, and TOC1 are able to control the transcription of other
genes because the regulatory proteins that they encode are
*transcription factors*, or master regulatory proteins that turn other
genes on and off. A transcription factor regulates a gene by binding
to a specific short DNA interval called a *regulatory motif*, or
*transcription factor binding site*, in the gene's *upstream region*, a
600-1000 nucleotide-long region preceding the start of the gene. For
example, CCA1 binds to AAAAAATCT in the upstream region of many genes
regulated by CCA1.

### Motif Finding

(mostly cited from Bioinformatics Specialization course)

The life of a bioinformatician would be easy if regulatory motifs were
completely conserved, but the reality is more complex, as regulatory
motifs may vary at some positions, e.g., CCA1 may instead bind to
AAGAACTCT. But how can we locate these regulatory motifs without
knowing what they look like in advance? We need to develop algorithms
for *motif finding*, the problem of discovering a “hidden message”
shared by a collection of strings.

```
Given a collection of strings, find a set of k-mers, one from each
string, that minimizes the score of the resulting motif.
```

To define scoring, consider a list of `t` DNA strings, where each
string has length `n`, and select a k-mer from each string to form a
collection Motifs, which we represent as a `t x k` motif matrix.

Our goal is to select k-mers resulting in the most “conserved” motif
matrix. One simple approach is to define Score(Motifs) as the number of
unpopular letters in the motif matrix. Another approach is to
construct a profile matrix as follow:

- Construct the 4 × k count matrix `Count(Motifs)` counting the
number of occurrences of each nucleotide in each column of the motif
matrix; the `(i, j)`-th element of `Count(Motifs)` stores the number of
times that nucleotide `i` appears in column `j` of Motifs.
- Divide all of the elements in the count matrix by `t`, the number of
rows in Motifs. This results in a profile matrix `P = Profile(Motifs)`
for which `P[i][j]` is the frequency of the i-th nucleotide in the j-th
column of the motif matrix. 

Every column of `Profile(Motifs)` corresponds to a probability
distribution. We can then score the motif based on the entropy of a
motif matrix, which is defined as the sum of the entropies of its
columns.

## Machine Learning Approach to Model the Genetic Basis of Disease Risks

Genome-Wide association studies:
- Predict disease risks
- Input: single-nuculeotide polymorphism (SNP) profiles of individuals

Evolutionary conservation:
- Detect conserved sequences
- Deleterious mutation (lower reproductive fitness) v.s. pathogentic (cause a disease)

Cell variables approach:
- One major goal of genomic medicine is to predict phenotypes, such as
  disease risks, from a genotype.
- By training models that predict how genotype defined by as a stretch
  of DNA sequence influences “cell variables,” such as concentrations of
  proteins, it hugely simplifies and modularizes the machine learning
  problem, and enables the exploration of therapies that target these
  crucial variables (cited from the paper written by M. Leung, et al).

Other way to classify ML:
- Genomics
- Proteomics
- Evolution
- Text Mining
- System Biology
- Microarray

## Note on the Hail Tech Talk

[Hail](https://hail.is/) is an open-source, scalable framework for exploring
and analyzing genomic data. [Its tech talk](https://www.youtube.com/watch?v=pyeQusIN5Ao)) describes a table 
for used for genetic data analysis `M`. `M[i,i]` represents a genotype of individual `i` at genome position `j`. 
 
```
M[i,j] = Genotype {
  “call”: “A/T”,
  “reads”: [10, 8],
  “quality”: 43,
  “p”: [43, 0, 52]
}
```
 
The metadata associated with individuals is ID (e.g.,  “NA12878”). The metadata associated with genome positions is following:
 
```
Genomic Locus { 
  “chromosome”: 1,
  “position”: 16123092,
  “reference: “A”,
  “alternate: “T”,
}
```
 
Additionally, there are two separate tables for metadata:
 
```
Locus-indexed table {
  “Gene”: “SHH”,
  “pred_impact” : “high”,
  “pop_frequency”: 0.102
}
 
ID-indexed table {
  “LDL”: 75.123,
  “ancestry”: “SAS”,
  “cohort”: “1KG”
}
```

TODO(kaneda): Describe operations to the tables.

## Cancer Related Info

- [Color Risk Model](https://github.com/color/risk-models): An implementation of the Claus risk model, which is a model of breast cancer risk from familial data.