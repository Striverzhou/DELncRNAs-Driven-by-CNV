# ceRNA-network-regulated-by-CNV-driven-DELncRNAs
In this study, CNV, lncRNAs, miRNAs, and mRNAs in The Cancer Genome Atlas (TCGA) data set were analyzed to explore a CNV-associated competing endogenous RNA (ceRNA) network related to cancer development.

# Start
# 1. Differential Expression Analysis
Samples derived from TCGA database were divided into the tumor or normal group. R package **edgeR** was applied to perform a differential gene expression between samples with the different groups to filter differentially expressed lncRNAs, miRNAs and mRNAs, |log2FC| > 1.5 and FDR < 0.05 were used as thresholds.

# 2. Identification of CNV-Driven DELncRNAs
#Merge all the segmentation files into a text.

**perl cnvMerge.pl** # Running the script in directory with all the segmentation files, it will produce a file named 'cnvMatrix.txt'.

**perl prepareChi.pl** # Preparing the input file for Chi-Square Test, the output file is 'chiInput.txt'.

**Rscript chi.R** #The chi-square test was used to analyze whether there is a significant difference in the CNV of lncRNA in normal tissue samples and tumor samples. The Bonferroni method was used to correct the p-value and lncRNAs with CNVs with adjusted. The output file is 'chiResult.txt'.

#Plot the circle map of lncRNA CNVs

**perl prepareRcircos.pl** # Preparing the input file for RCircos, the output file is 'lncRNA_Label.txt'.

**Rscript RCircos.R** # Plot the circle map

#the upregulated expression genes with copy number amplification and the downregulated expression genes with copy number deletion were selected as candidate CNV-driven DElncRNAs.

**Rscript cnv-expr-kruskal.R** # Kruskal-Wallis test was used to analyze correlation between CNV and lncRNA expression. The output file is 'cnv-expr-ks.txt'. Finally, the genes with pvalue < 0.05 were selected as the CNV-Driven DELncRNAs.

# 3. Construction of CNV-Driven ceRNA Network
#miRNAs that interacted with CNV-driven DElncRNAs were predicted through LncBase database and starBase database, and the predicted miRNAs were then intersected with DEmiRNAs to select the DEmiRNAs that interacted with CNV-driven DElncRNAs. Pearson correlation analysis was conducted to calculate the correlation between the expression of CNV-driven DElncRNAs and the DEmiRNAs mentioned earlier. The DEmiRNAs with r < −0.2 and p < 0.05 were selected as downstream CNV-associated DEmiRNAs regulated by lncRNAs.

#Similarly, mRNAs that had an interactive relationship with CNV-associated DEmiRNAs were predicted through starBase, miRDB, mirDIP, and TargetScan databases. The predicted mRNAs mentioned earlier were intersected with DEmRNAs, and the DEmRNAs with r < −0.2 and p < 0.05 were regarded as downstream CNV-associated DEmRNAs regulated by miRNAs.

#Based on the results mentioned earlier, lncRNAs, miRNAs, and mRNAs with complete regulatory relationships were selected to establish a corresponding ceRNA network for subsequent analysis. The ceRNA network was the CNV-associated ceRNA network visualized by the Cytoscape. 
