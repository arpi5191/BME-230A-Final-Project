# Space Gene Expression Analysis

Implements Principal Component Analysis (PCA) on a NASA dataset to identify the features contributing most significantly to the mice's gene expression. Performs differential gene expression analysis based on condition, library preparation, and strain features to uncover insights about mice in space compared to those on the ground.

# Instructions:
  1) Download metadata.csv and data.csv to a platform like Jupyter Notebook or Google Colab.
  2) Download the PCA.ipynb notebook and upload it to the same platform.
  3) Open the PCA.ipynb notebook and run it to perform Principal Component Analysis (PCA) on the dataset. The analysis will identify the most significant factors—condition, library preparation (libprep), and strain—that contribute to gene expression variation between mice in space and on the ground.
  4) Download the following R scripts: DeseqRNA_condition.r, DeseqRNA_libprep.r, DeseqRNA_strain.r. Open RStudio and run each script to analyze the differential gene expression based on the condition, libprep, and strain features from the dataset.
