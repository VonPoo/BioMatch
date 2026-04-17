# BioMatch

BioMatch is a data-driven framework for comprehensive sample identification.

## Features

- Generate k-mer panels from reference genomes and population VCFs
- Count k-mers in sample FASTA files
- Evaluate sample identity
- Support for FASTQ, VCF and PLINK format files

## Installation

```bash
# Install from conda (推荐，包含所有R包依赖)
conda install -c bioconda biomatch

# 或从本地包安装 (注意：需要手动安装R依赖)
cd /path/to/biomatch
pip install .
```

### 依赖说明

BioMatch依赖以下软件和包：

- Python依赖: pyfaidx, pandas, numpy, biopython
- R依赖: tidyverse, caret, e1071, pls, deepKin
- 其他工具: bcftools, tabix, samtools, plink, plink2, parallel

## Usage

BioMatch supports several processing modes:

### 1. Generate Panel Only

```bash
biomatch --gen-panel --ref reference.fa --pop-vcf population.vcf --chr-set 33 --panel-name panel_name.fa
```

### 2. Generate Panel, Count & Eval

```bash
biomatch --gen-panel --ref reference.fa --pop-vcf population.vcf --chr-set 33 --panel-name panel_name.fa --count samples_dir --count-db count_results --eval-result eval_results
```

### 3. Count & Eval using Existing Panel (并行处理)

```bash
# 使用-t参数指定并行处理的线程数
biomatch --panel-name panel_name.fa --count samples_dir --count-db count_results -t 20 --eval-result eval_results
```

> **注意**: `-t` 参数控制并行处理的线程数，建议根据系统资源设置合理的值

### 4. Deepkin Eval on VCF/PLINK

```bash
biomatch --match samples.vcf --species human --eval-result eval_results
```

### 5. Default Eval on Count Results

```bash
biomatch --count-db count_results --eval-result eval_results
```

## License

MIT