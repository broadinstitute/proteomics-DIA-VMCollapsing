# proteomics-DIA-VMCollapsing

REQUIRES PYTHON 3.8 TO RUN 
Performs VM collapsing on Spectronaut and DIA-NN data

In yaml file:
  report: File name ('report.tsv' for DIA-NN, report name for Spectronaut)
  fasta path: Point to fasta file 
  engine: 'SN' for Spectronaut, 'DIANN' for DIA-NN
  SN/DIANN:
    heavies: List of heavy modifications 
    target_ptm: PTM to be collapsed on
    experimental: List of experimental modifications
    ptm_confidence_threshold: PTM site localization threshold that should be applied

In the script: Paste path to directory your files are for variable wd at the top of the script and hit run.
    
