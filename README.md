This repository contain code used to produce results shown in article [Methylation-independent gene regulation by methyl-CpG-binding protein 2]()

Following published datasets has been used in this study 


| **SI. No** | **Sequencing data type**   |**Source** | Datasets |
|-----------------|--------------------------------------|-----------------|----- |
| 1\.        | CUT&RUN (Adult Brain Cortex of *Mecp2*, *Dnmt3a* WT and KO mice)         | [GSE150538](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150538)       | Raw data was analysed                          |
| 2\.        | MECP2 CUT&RUN (Adult hippocampus of *Mecp2*-TG1, *Mecp2* KO, and WT) | [GSE213752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213752)       | Processed Bigwig was used                          |
| 3.         | MECP2 ChIP-seq     (Frontal cortex)                                                  | [GSE67293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67293)     |  Raw data was analysed                        |
| 4.         | MECP2 ChIP-seq  (Forebrain of WT and *Mecp2* KO)                                                | [GSE139509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139509)   |  processed file was used                        |
| 5\.        | BS-seq (Brain Cortex of 10 week old *Dnmt3a* cKO and WT mice)          | [GSE103214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103214)         | processed file was used                        |
| 6\.        | MECP2 CUT&RUN ( WT and G118E )                                       | [GSE243009](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243009)  | processed file was used                        |
| 7\.        | RNA-seq (Supplementary Table S1)                                     | [Boxer et. al, 2020](https://www.sciencedirect.com/science/article/pii/S109727651930810X?via%3Dihub#app2) |  processed file was used                        |
| 8\.        | snRNA-seq (Supplementary Table S7)                                     | [Renthal et. al, 2018](https://www.nature.com/articles/s41593-018-0270-6#citeas) |  processed file was used                        |

*Note : Reserve node before running a shell script*
```
srun --partition=super --nodes=8 --pty --time=5:00:00 /bin/bash
```
**Step1 :**
```bash
$ cd MECP2_Project/scripts/shell_scripts/
[shell_script]$ bash getRequiredFiles.sh
```

**Step2 :**

```bash
[shell_script]$ perl Process_CUTandRUN.pl ../../data/CUTandRUN_seq_Data.txt
[shell_script]$ bash MECP2_MECP2WT_vs_IgG_and_MECP2KO_macs2.sh
[shell_script]$ bash process_peak_files.sh
```

**Step3 :** open `R/MECP2_project.qmd` in Rstudio

Run all the code in the data processing step.

**Step4 :**
```bash
[shell_script]$ perl callPeaks_Macs2.pl
[shell_script]$ bash bamCoverage.sh
[shell_script]$ perl Process_MECP2_ChIPseq.pl
```
**Step5 :**

```bash
[shell_script]$ bash MECP2_G118E_CNR.sh
[shell_script]$ bash Process_GEOdatasets.sh
[shell_script]$ bash Process_GSE139509.sh
[shell_script]$ bash Process_ENCODE.sh
[shell_script]$ bash Figure1.sh
[shell_script]$ bash Figure2.sh
[shell_script]$ bash Figure3.sh
[shell_script]$ bash Supplementary_figure_9.sh
```