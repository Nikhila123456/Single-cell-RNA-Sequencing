

Seurat allows to easily explore QC metrics and filter cells based on various quality criteria. A few commonly used QC metrics are given below:

1. Transcript (nCount) and gene (nFeature) abundance
	- Low-quality cells or empty droplets will often have very few genes/transcripts
	- Cell doublets or multiplets may exhibit an aberrantly high gene count
2. The percentage of reads that map to the mitochondrial genome
	- Low-quality / dying cells often exhibit extensive mitochondrial contamination
	- Mitochondrial QC metrics is calculated as percentage of counts originating from a set of mitochondrial genes (i.e. all genes starting with MT)
3. Ribosomal contents
	- Less than 50% ribosomal content is often preferred
