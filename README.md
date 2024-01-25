# BAHD_SCPL_synteny

Scripts to create the synteny network of BAHD and SCPL genes. 

- `01_network_construction_BAHD_SCPL.R`: focuses on integrating and 
  visualizing synteny information from multiple sources and conducts quality 
  control to assess the reliability of the identified relationships.
  
  In detail the following steps are taken:
  
  

    1. Loading and extracting genes:
        - Read a table containing orthogroups information from Orthofinder and MCL outputs.
        - Define specific orthogroups related to BAHD and SCPL genes.
        - Create a list of genes from Orthofinder and MCL results and saves it.

    2. Synteny analysis using i-ADHore:
        - Define functions (`find_synteny` and `find_tandem`) to extract 
          synteny and tandem gene pairs from i-ADHore results.
        - Apply these functions to i-ADHore outputs generated from Orthofinder 
          and MCL inputs.
        - Generate binary matrices (`bin_mat` and `bin_mat_tandem`) to 
          represent the presence of synteny or tandem relationships.

    3. Synteny analysis using MCScanX:
        - Define functions (`find_synteny_mcscanx` and `find_tandem_mcscanx`) 
          to extract synteny and tandem gene pairs from MCScanX results.
        - Apply these functions to MCScanX outputs generated from Orthofinder 
          and MCL inputs.
        - Generate binary matrices (`bin_mat_mcscanx` and 
          `bin_mat_tandem_mcscanx`) to represent the presence of synteny or 
          tandem relationships.
          
    4. Integration and visualization:
        - Combine synteny information from i-ADHore and MCScanX from Orthofinder 
          and MCL inputs.
        - Generate a network plot (`graph_from_adjacency_matrix`) for 
           both tandem and synteny relationships.
        - Determine components in the complete graph.

    5. Quality Control and Visualization:
        - Check reliability of synteny information based on number of 
          neighboring genes in scaffolds.
        - Evaluate and visualize the reliability of the connections based on 
          the number of neighboring genes.
        - Remove genomic regions that do not link to others and plot separate 
          networks for linking and non-linking features.
        - Quantify and visualize the number of genes on chromosomes/scaffolds
          for linking and non-linking features.