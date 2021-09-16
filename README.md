# mediation_GWAeQTL

A Nextflow pipeline used for mediation analysis for [cegwas2-nf](https://github.com/AndersenLab/cegwas2-nf) output.

## Execution of pipeline using Nextflow
```
git clone https://github.com/AndersenLab/mediation_GWAeQTL.git

cd mediation_GWAeQTL

nextflow mediation_GWAeQTL.nf --cegwas2dir = <path to cegwas2-nf result folder> --traitfile = <path to cegwas2-nf traitfile> 

```


## Required software packages that should be in users PATH

1. [nextflow-v19.07.0](https://www.nextflow.io/docs/latest/getstarted.html)
2. [R-tidyverse-v1.2.1](https://www.tidyverse.org/)
3. [R-mediation-v4.5.0](https://github.com/kosukeimai/mediation)
4. [R-MultiMed-v2.4.0](https://github.com/SiminaB/MultiMed)



## Pipeliine parameters

* --cegwas2dir

cegwas2-nf output folder

* --traitfile

cegwas2-nf input traitfile

* --transcripteQTL

Transcript-level eQTL detected in [WI-Ce-eQTL](https://github.com/AndersenLab/WI-Ce-eQTL) 

* --transcript_exp

Transcript-level expression data used for mapping eQTL in [WI-Ce-eQTL](https://github.com/AndersenLab/WI-Ce-eQTL) 

* --out

Add result folder name and path. The default is "mediation-*date*" in the input cegwas2-nf result folder

 
 

