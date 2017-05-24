# MasterFinalProject
Scripts used in the Final Project: ACO for predict gene interactions from expression data (Master's Degree in Biostatistics and Bioinformatics at UOC)

In this study ACOGeneInteractions was implemented based on the previous ACOPants already implemented in Python. 
In order to use ACOPants for predicting gene interactions some modifications were done:

  1.  The input for ACOPants is a list of coordinates (x, y) of the nodes, and providing a length function to the algorithm this is able to   calculate the distances from node i to j. 
  In ACOGeneInteraction the input is the correlation matrix from raw data (DNAmicroarrays), where the column’s names represent the names of   the genes (nodes), and the distances between two nodes are given by the pair-wise correlation value in the matrix. In this case, nodes     have no coordinates so we have to modify the algorithm to accept this data format.

  2. ACOPants is implemented using the TSP approach and each agent (ant) visits each node (gene) once. 
  In the modified version, instead of going once through each node (gene) the algorithm will be forced to visit each node a specific number   of times which is set by a parameter called ‘times’. The purpose of this modification is to overcome the limitation of the algorithm        found in other studies: ACO is able to find as many interactions as number of nodes are present in the graph. The idea is to force the      algorithm to visit alternatives routes in order to find more gene interactions. 

  3. The outputs when using ACOPants is the list of nodes visited sequentially during the tour and the distance of the complete optimal       tour. 
  ACOGeneInteraction has been modified so also the output shows a list of interactions, pairs of genes that interact between them.
