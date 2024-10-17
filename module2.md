### **RNA-Seq Data Analysis**

---

### **7.3.1 Read Mapping**
- **Initial QC**: After RNA sequencing, review metrics such as **total read count**, **quality score distribution**, and **GC content**. Tools like **subSeq** can help assess whether the sequencing depth provides adequate detection power.
  
- **Eukaryotic vs. Prokaryotic Mapping**: 
  - **Prokaryotes**: Use aligners like **Bowtie** or **BWA** for contiguous mapping, as no splicing occurs.
  - **Eukaryotes**: Mapping is more complex due to **intron splicing**, requiring specialized aligners for spliced reads.
  
- **Mapping Approaches**:
  1. **Annotation-guided mappers**: Align reads to pre-built transcript databases (e.g., **RUM**, **SpliceSeq**). While accurate, they can't detect novel transcripts and may result in higher multi-mapping rates.
  2. **Ab initio spliced mappers**: Identify splicing junctions without relying on annotations.
     - **Exon-first methods**: Map exonic reads first, then predict splice junctions (e.g., **TopHat**, **MapSplice**).
     - **Seed-and-extend methods**: Begin with small k-mer seeds, extending the alignment to detect splice junctions (e.g., **STAR**, **HISAT2**).

- **Long-Read RNA-Seq**: Tools like **minimap2**, **deSALT**, and **uLTRA** handle long RNA reads (e.g., from PacBio or ONT platforms).

- **Key QC Metrics**:
  - **Mapping efficiency**: 70-90% of reads should map to the genome.
  - **rRNA reads**: A successful rRNA depletion step results in low rRNA reads (<10%).
  - **Duplicate reads**: High duplicates (up to 40-60%) are often due to **PCR overamplification** or highly expressed genes but may not always require removal.

- **De Novo Transcriptome Assembly**: For species without a reference genome, use de novo assemblers like **Trinity**, **rnaSPAdes**, or **StringTie2**. This is essential for identifying **novel transcripts** and for heavily fragmented genomes, such as tumor cells.

---

### **7.3.2 Quantification of Reads**
- **Read Counting Tools**:
  - Tools like **featureCounts** and **htseq-count** require aligned reads (SAM/BAM) and use gene feature annotations (GFF/GTF) for quantification.
  - **RSEM**, **Cufflinks**, and **eXpress** use model-based approaches, assigning reads probabilistically to transcripts.

- **Mapping-Independent Quantification**:
  - Tools like **kallisto** and **Salmon** use **pseudo-alignment** of k-mers to transcripts, significantly speeding up the process without full alignment.

---

### **7.3.3 Normalization**
- **Factors Influencing Read Count**:
  1. **Sequencing depth**: Deeper sequencing leads to more reads, which needs to be corrected.
  2. **Gene length**: Longer transcripts generate more reads.

- **Normalization Methods**:
  - **RPKM/FPKM**: Corrects for sequencing depth and gene length.
  - **Total count normalization**: Adjusts by dividing by the total number of mapped reads.
  - **Quantile normalization**: Ensures uniform read count distributions across samples.
  - **Advanced Methods**:
    - **DESeq2 (Relative Log Expression, RLE)**: Computes scaling factors based on the median of read count ratios for each gene across samples.
    - **edgeR (Trimmed Mean of M-values, TMM)**: Excludes highly expressed genes and computes scaling factors to adjust library size.

- **Iterative Normalization**: Methods like **DEGES** or **PoissonSeq** repeatedly remove differentially expressed genes (DEGs) to improve normalization accuracy.

---

### **7.3.5 Identification of Differentially Expressed Genes (DEGs)**
- **Data Distribution**: RNA-Seq data follows a **Poisson** or **negative binomial distribution** due to variance scaling with expression level.

- **DE Tools**:
  - **DESeq2**: Uses negative binomial generalized linear models and a Wald test to detect DEGs.
  - **edgeR**: Also based on negative binomial models, with exact tests and GLM likelihood ratio tests.
  - **limma**: A linear modeling approach adapted for RNA-Seq (voom function).
  - **Others**: **Cuffdiff2**, **EBSeq**, and **NOISeq** for Bayesian or nonparametric methods.

---

### **7.3.6 Multiple Testing Correction**
- **Multiple Comparisons**: RNA-Seq involves testing thousands of genes, increasing the risk of false positives.
  - **Bonferroni correction**: Controls the **family-wise error rate (FWER)** but is conservative.
  - **FDR (False Discovery Rate)**: **Benjamini-Hochberg** correction is more practical, allowing a small proportion of false positives while boosting statistical power. Genes are assigned **q-values**, where q < 0.05 suggests a 5% chance of false positives.

---

### **7.3.7 Gene Clustering**
- **Clustering Algorithms**:
  1. **Hierarchical clustering**: Builds a dendrogram based on gene expression similarities.
  2. **k-means clustering**: Groups genes into predefined clusters (k).
  - Both use **Pearson correlation coefficient** as a similarity measure. Tools like R can perform clustering for gene expression patterns and sample comparisons.

---

### **7.3.8 Functional Analysis of Identified Genes**
- **Functional Enrichment**: After identifying DEGs, tools like **Enrichr**, **GOseq**, or **DAVID** connect genes to **biological pathways** (e.g., KEGG, Reactome) and **Gene Ontology (GO)** terms.
  - **Gene Set Enrichment Analysis (GSEA)** uses the entire gene set instead of a filtered list, increasing sensitivity in detecting weaker signals.

- **Gene Networks**: Tools like **Cytoscape** reconstruct **gene regulatory networks**, visualizing interactions and co-expression patterns.

---

### **7.3.9 Differential Splicing Analysis**
- **Splicing Event Detection**: 
  - Methods like **DEXSeq** and **JunctionSeq** analyze exon usage.
  - Event-based tools like **MISO** and **rMATS** quantify splicing events (e.g., exon skipping, retained introns) using **percent-spliced-in (Psi)**.
  - Isoform-level tools like **Cuffdiff2** reconstruct transcript isoforms and test for differential expression.

---

### **7.4 Visualization of RNA-Seq Data**
- **Visualization Tools**:
  - Basic plots: **Histograms, boxplots, density plots**, etc., are useful for exploratory data analysis.
  - **IGV** or **UCSC Genome Browser**: Provide visual inspection of mapped reads, especially for splicing events.
  - **SpliceSeq**, **JunctionSeq**, and **DiffSplice** offer splicing-specific visualization.

---

### **7.5 RNA-Seq as a Discovery Tool**
- **Hypothesis Generation**: RNA-Seq identifies DEGs, novel transcripts, and splicing events, leading to new biological hypotheses.
  - **Validation**: RNA-Seq results are validated using **qPCR** or **Western blotting**.
  - **Clinical Applications**: RNA-Seq assists in disease classification (e.g., breast cancer subtypes) and resolving genetic uncertainties (e.g., hereditary cancer diagnosis).

