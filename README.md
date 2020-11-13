# Cross Methylome Omnibus (CMO) test

CMO is a gene-level association test that can identify many significant and novel genes ignored by many benchmark methods.  Specifically,  CMO integrates genetically regulated DNAm in enhancers, promoters, and the gene body to identify additional disease-associated genes. Please cite the following manuscript for CMO analysis:

>  Wu et al. A gene-level methylome-wide association analysis identifies novel Alzheimer's disease genes. bioRxiv. doi:https://doi.org/10.1101/2020.07.13.201376. 

In this repo, we provide the following sources.

* CMO: the software for running CMO test

* Codes: all source codes for replicating the results present in the manuscript.

  

## Outline

1. [Installation](#Installation)
2. [Typical analysis and output](#Analysis)
3. [Command-line parameters](#Command)
4. [FAQ](#FAQ)



## <a name="Installation"></a>Installation

- Download and unpackage the CMO package from GitHub. Download through [this link](https://github.com/ChongWuLab/CMO/releases/tag/v1.0) or by the following commands:

  ~~~
  wget https://github.com/ChongWuLab/CMO/archive/v1.0.zip
  unzip v1.0.zip
  cd CMO-1.0
  ~~~

- Download and unzip the [constructed weights](). (The link will be available after the manuscript has been accepted.)


## <a name="Analysis"></a>Typical analysis and output

The CMO analysis takes pre-computed DNA methylation prediction models (included in the CMO package), enhancer-promoter interactions (included in the CMO package), and GWAS summary data to estimate the association between a gene and the trait of interest. We will use the [IGAP Alzheimer's summary data](http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php) (Lambert et al. 2013) as an example to illustrate how to use our methods. This example assumes you have set up the required environment and data, as illustrated in the previous section.

### Input: GWAS summary statistics

We write a wrapping code and try to support all publically available GWAS summary data. Please revise the column names to the commonly used ones  (such as BETA, SE, Z, A1, A2, SNP) if you meet some problems. We also calculate Z score vectors automatically if the raw GWAS summary data failed to include it. 

**Note:** The performance of CMO depends on the density of summary-level data. We highly recommend running CMO with raw summary-level data. Pre-process steps such as pruning and restricting to top SNPs may harm the performance.

### Performing the CMO

After we prepared the data, we can run CMO via the following single line.

```
Rscript CMO.R \
--sumstats ./Example/IGAP_chr22.txt \
--out ./Example/ \
--weights_dir /gpfs/research/chongwu/Chong/CMO/
--chr_id 22
```

This should take around one or two minutes, and we will see some intermediate steps on the screen. If everything works, the results will be saved into ./Example/. Briefly, we save the results of CMO (res_CMO_CHR22.txt) and individual MWAS results (res_MWAS_CHR22.txt).



### Output: Gene-disease association

The results are stored in a user-defined output file. For illustration, we explain the meaning of each entry in the first two lines of the output.

| Col. num. | Column name | Value           | Explanations                                          |
| --------- | ----------- | --------------- | ----------------------------------------------------- |
| 2         | CHR         | 22              | Chromosome ID                                         |
| 1         | geneID      | PNPLA3          | Feature/gene identifier, taken from gene_list file    |
| 2         | ensembl     | ENSG00000100344 | Ensemble ID                                           |
| 3         | P0          | 44319619        | Gene start                                            |
| 4         | P1          | 44360368        | Gene end                                              |
| 5         | n_enhancer  | 5               | Number of enhancers that linked to the gene           |
| 6- 11     | CpG sites   |                 | Information for how many CpG sites and SNPs are used. |
| 12        | CMO         | 0.000122        | P value for CMO test                                  |
| 13        | Runtime     | 1.3             | Running time for this gene (in second)                |

**Note:** We only store the results for genes with external weights. The genes without external weights will be ignored in the output file. We also save the results for standard MWAS. However, we only run MWAS for the CpG sites that are linked to genes. 

## <a name="Command"></a>Command-line parameters

### CMO.R

| Flag          | Usage                                                        | Default  |
| :------------ | ------------------------------------------------------------ | -------- |
| --sumstats    | Summary statistics (rds file and must have SNP and Z column headers) | Required |
| --out         | Path to output directory                                     | Required |
| --weights_dir | The weights directory                                        | Required |
| --chr_id      | The chromosome ID. We recommend parallelling the computations by chromosomes | Optional |

## <a name="Analysis"></a>FAQ

If you have questions, please submit an issue. We will summarize commonly asked questions here. 



## License

Maintainer: [Chong Wu](http://wuchong.org/index.html) (cwu3@fsu.edu)

[MIT](http://opensource.org/licenses/MIT)

Copyright (c) 2013-present, Chong Wu (cwu3@fsu.edu)