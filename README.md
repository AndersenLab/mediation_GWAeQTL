# mediation_GWAeQTL

A Nextflow pipeline used for mediation analysis for [cegwas2-nf](https://github.com/AndersenLab/cegwas2-nf) or [NemaScan](https://github.com/AndersenLab/NemaScan) output.

Citations: https://doi.org/10.5281/zenodo.6595334,  https://doi.org/10.1038/s41467-022-31208-4

## Execution of pipeline using Nextflow
```
git clone https://github.com/AndersenLab/mediation_GWAeQTL.git

cd mediation_GWAeQTL

nextflow mediation_GWAeQTL.nf --gwa="<cegwas2nf or nemascan>" --gwa_dir="<path to GWA result folder>"

```


## Required software packages that should be in users PATH

1. [nextflow-v19.07.0](https://www.nextflow.io/docs/latest/getstarted.html)
2. [R-tidyverse-v1.2.1](https://www.tidyverse.org/)
3. [R-mediation-v4.5.0](https://github.com/kosukeimai/mediation)
4. [R-MultiMed-v2.4.0](https://github.com/SiminaB/MultiMed)
5. [R-data.table-v1.12.8](https://cran.r-project.org/web/packages/data.table/index.html)



## Pipeliine parameters

* --gwa

Choose between "cegwas2nf" or "nemascan" according to the GWA pipeline used for QTL mapping 

* --gwa_dir

GWA mapping output folder

* --transcripteQTL

Transcript-level eQTL detected in [WI-Ce-eQTL](https://github.com/AndersenLab/WI-Ce-eQTL) 

* --transcript_exp

Transcript-level expression data used for mapping eQTL in [WI-Ce-eQTL](https://github.com/AndersenLab/WI-Ce-eQTL) 

* --out

Add result folder name and path. The default is "mediation-*date*" in the input GWA mapping result folder

 
 

