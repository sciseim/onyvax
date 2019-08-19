# onyvax_MS
Code and data for the MS: 
*Characterisation of two cancerous and two benign cells lines derived from the primary tumour site in prostate cancer patients*.

**Summary**:
featureCounts can be found in 'featureCounts'

*`Scripts used for various analyses downstream of TopHat`*:

- `Differential gene expression: 1-prepare-BAM-files-for-GFOLD.sh`
Copy BAM files from TopHat from each sample into a GFOLD directory.

- `Differential gene expression: 2-gfold.sh`
Run GFOLD.

- `featureCounts normalisation and heat map generation: featureCounts-to-heatmap.R`
Normalise and draw a heat map.
