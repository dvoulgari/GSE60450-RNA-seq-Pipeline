# 🧬 GSE60450 RNA-seq Pipeline

A concise RNA-seq workflow for the **GSE60450** mouse mammary gland dataset, including SRA processing, QC, trimming, alignment, and differential expression.

---

## 🔧 Pipeline Steps

| Step | Description | Tools |
|------|-------------|--------|
| **1. SRA → FASTQ** | Convert `.sra` files to `.fastq.gz` | SRA Toolkit |
| **2. QC** | Quality check on raw reads | FastQC |
| **3. Trimming** | Adapter removal | Cutadapt |
| **4. Alignment** | STAR genome alignment + gene counts | STAR |
| **5. DE Analysis** | Differential expression & plots | R, DESeq2 |

---

## 📂 Repository Files

- sra_processing.sh # SRA → FASTQ conversion
- rnaseq.sh # QC, trimming, alignment pipeline
- deseq_analysis.R # DESeq2 differential expression

---

## ▶️ Running the Pipeline

**Bash:**
```bash
bash sra_processing.sh
bash rnaseq.sh
source("deseq_analysis.R")
