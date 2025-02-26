# A-stress-response-axis-in-TAMs-dictates-phenotype-heterogeneity-and-immunotherapy-evasion

This repository provides code to reproduce the figures of (TODO: insert doi here)

Raw data can be found on GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290434

Docker images to run the jupyter notebooks can be found on dockerhub:
- Single-cell RNA sequencing: primelabzurich/single_cell:1.0
- Spatial RNA sequencing: primelabzurich/spatial:1.0


How to reproduce figures:

1. clone git repository:
   ```
   git clone git@github.com:primelab-zurich/A-stress-response-axis-in-TAMs-dictates-phenotype-heterogeneity-and-immunotherapy-evasion.git
   ```
2. pull docker images:
   ```
   docker pull primelabzurich/single_cell:1.0
   ```
   ```
   docker pull primelabzurich/spatial:1.0
   ```
3. run docker container with repo mounted (insert the path to cloned repo):
   ```
   docker run --rm -it -p 7777:8888 -v <path to cloned repo>:/home/sc primelabzurich/single_cell:1.0
   ```
   ```
   docker run --rm -it -p 7777:8888 -v <path to cloned repo>:/home/sc primelabzurich/single_cell:1.0
   ```
