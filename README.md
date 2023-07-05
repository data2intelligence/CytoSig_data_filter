Example code of filter cytokine response by correlation between ligand or receptor expression and cytokine target expression values.  
  
Please CD into the src folder and run "./download.py" first to download all data files.  

There are 61 datasets from the TCGA and GTEx cohorts.  
So, you may run "./correlation.py inx 61" where inx is a number between 0 and 60 to compute the correlation results for each dataset.  
We used the NIH high performance cluster (HPC) for parallel computation as "./hpc_submit.py 61". You may re-write this file for your local HPC.  

After computing the correlation results for each dataset, you can run "./correlation.py" to generate the filtered cytokine response data.  