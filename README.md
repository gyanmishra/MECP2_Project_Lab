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
```
$ git clone https://github.com/gyanmishra/MECP2_Project_Lab.git
$ cd MECP2_Project_Lab/scripts/shell_scripts/
```

**Step2 :** prerequisite
```bash

[shell_script]$ conda create deeptools3.1.2 
[shell_script]$ conda activate deeptools3.1.2
[shell_script]$ conda install deeptools==3.1.2
[shell_script]$ bash getRequiredFiles.sh
```


**Step3 :**

```bash
< Input directory > /archive/OBI/Neuroinformatics_Core/Stroud_lab/shared/MECP2_project_Gyan/

[shell_script]$ bash Figure1.sh < Input directory > 
[shell_script]$ bash Figure2.sh < Input directory >
[shell_script]$ bash Figure3.sh < Input directory >
[shell_script]$ bash Supplementary_figure_9.sh 
```

**Step4 :** open `R/MECP2_project.qmd` in Rstudio

Run all the code one by one for each figure including data processing step. 
