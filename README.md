# Gene Set Enrichment Analysis (GSEA)

Implementation of the GSEA algorithm as described in the original paper: Subramanian, Aravind, et al. "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Available at: <http://www.pnas.org/content/102/43/15545.abstract>

#### Files in the repository:

* gsea.py: The main file with the gsea function. Also includes some supplementary functions.
* leukemia.txt: Example gene expression data.
* pathways.txt: Example gene sets data.
* results.txt: Example output data, obtained by running gsea function with 'leukemia.txt' and 'pathways.txt' files.

#### Instructions:

* Run the gsea.py file.
* Call the gsea(expressionProfiles, geneSets) function. (Example: gsea('leukemia.txt', 'pathways.txt')).
* Results.txt file is created. 
