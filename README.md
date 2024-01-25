# Analysis of syntenic relationships of BAHD and SCPL-containing genomic regions using a network approach

Scripts to create the synteny network of BAHD- and SCPL-containing genomic regions. 

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
          
- `01_network.qmd`: focuses on detecting communities, calculating topology
  parameters, and clique determination.
  
  In detail the following steps are taken:
  
  
    1. Preparing the network:
        - Library Loading: The code loads essential R libraries 
            (`dplyr`, `ggplot2`, `ggpubr`, `igraph`, and `scales`) for data 
            manipulation, visualization, and network analysis.
        - Network Construction: A binary matrix representing genomic relationships 
            is used to construct a network. The code focuses on the five 
            largest components with syntenic information from the 
            `bin_mat_complete_cut` network.
        - Manual Edge Removal: Certain edges in the network are manually 
            removed due to misclassification. 
            
    2. Community Detection:
        - Various community detection algorithms 
            (fastgreedy, walktrap, leading.eigenvector, label.propagation, 
            multilevel, infomap) are applied to identify clusters or communities 
            within the modified network (`net_cut`).
        - Biological Information Integration: Additional information on 
            orthogroups, MCL memberships, and family (BAHD/SCPL) assignments is 
            added to the network nodes. The results, including community 
            memberships and biological annotations, are stored in a data frame
            (`membership_community_gene_cluster_added.txt`) for further analysis.
            
    3. Topology parameter calculation:
        - Function definitions: Define a set of helper functions 
            for calculating topological parameters of networks.
            `calculate_topological_parameters` function computes various 
            network metrics for a given component.
            `get_topology_for_components` function iterates through 
            orthogroups, decomposes the network into components, 
            and calculates topological parameters for each component.
        - Function application: The functions are applied to orthogroups of 
            BAHD and SCPL genes using the `get_topology_for_components` function.
            The results are stored in the `topology_bahd` and `topology_scpl` 
            variables.
        - Plotting: Include a generic plotting function, 
            `plot_topology_parameter`, that generates scatter plots for 
            specified topological parameters.
            The various topological parameters are visualized, 
            including average local efficiency, average path length, 
            betweenness, closeness, constraint, diameter, degree, 
            eccentricity, edge betweenness, eigencentrality, global 
            efficiency, harmonic centrality, neighborhood size, radius, 
            strength, and shortest paths.
        - Linear models: Bootstrap analysis and t-tests are performed for 
            topological parameters to compare differences between 
            BAHD and SCPL gene families.

    4. Clique determination:
        - `find_largest_non_overlapping_cliques` function:
            Find the largest cliques in a given graph (`g`) without overlap.
            Repeatedly identify the largest clique, remove its nodes, 
            and continue until all nodes are assigned to cliques.
        - `greedy_clique Function`: Applies the 
            `find_largest_non_overlapping_cliques` function on each network 
            component of a given graph (`g`).
            Run the process multiple times (`n_rep` repetitions) to randomly 
            select the largest cliques.
            Select the optimal result based on the number of cliques and a 
            high median membership of genomic regions.
        - Application on Two Different Inputs:
            1. MCL groups and orthofinder orthogroups input:
                Convert the weighted adjacency matrix derived from 
                MCL groups and Orthofinder orthogroups into a binary matrix/graph.
                Apply the greedy_clique function on the graph.
                Assign components and cliques to genomic regions based on the 
                determined cliques.
            2. Only Orthofinder orthogroups input:
                Similar to the first input, convert a weighted matrix into a 
                binary matrix but with specific considerations (only edges
                that have Orthofinder information).
                Apply the `greedy_clique` function on the graph derived from 
                Orthofinder orthogroups only.
                Assign components and cliques to genomic regions based on the 
                identified cliques.

    5. Export and Visualization:
        - Export the results to a CSV file (`membership_community.csv`) 
            containing information about components and cliques for each 
            genomic region.
        - Export the network file derived from the adjacency matrix 
            (`bin_mat_cut_graphml_manually_removed_edges.xml`).
        - Plot the network using various community detection algorithms 
            (fastgreedy.community, walktrap.community, 
            leading.eigenvector.community, label.propagation.community, 
            multilevel.community, infomap.community, og, mcl, component_og_mcl).
            Each plot represents a network component, and nodes are colored 
            based on the assigned community/component.

