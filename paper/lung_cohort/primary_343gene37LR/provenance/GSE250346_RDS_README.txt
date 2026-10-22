Contents have been revised since the original (pre-publication) version of this work. Please see important notes at the end of this README for details.

---XENIUM OUTPUT---
- Full output from Xenium have been included as .tar.gz files in order to allow the data to be viewed in the Xenium Explorer.
- Directory names from Xenium are in the format: [instrument serial number]__[slide serial number]__[sample]__[date]__[time]. 
- Samples on TMAs 1-4 were run separately from TMA5, and have slightly different output formats. See paper for details.
- Information about output files can be found here: https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/understanding-xenium-outputs
NOTE: GSE250346_output-XETG00048__0003789__VUILD104MFVUILD48LFVUILD105LF__20230308__003731.tar.gz contains Xenium output for three samples (VUILD104MF/VUILD105LF/VUILD48LF) which were merged due to FOV overlap. These samples have been split out in the transcripts.csv.gz files.


---SEURAT OBJECT (RDS FILE)---
- Seurat object containing expression data for nuclei and cell-level metadata.


---EMBEDDINGS FILES (GRAPHSAGE ANALYSIS)---
- The npy file for each sample contains the N by 50 embedding matrix where N is the the number of transcripts. The rows in the npy file corresponds to rows in the node-meta csv file where transcript locations, and gene labels can be found.


---IMAGES---
- Files ending in .svs.tif.gz are high res whole slide TIFFs of TMA slices, please see accompanying powerpoint slides for sample assignments- 
- Files ending in .tif.gz are single sample high res TIFFs that have been registered to be in a common coordinate system with the associated Xenium data


---PATHOLOGY ANNOTATIONS---
- GSE250346_HE_annotations.tar.gz is a tar containing pathology annotations in geojson format and a csv indicating which cells fall within which annotated regions. There is one geojson file per annotated region and they are organized into directories by sample. The geojsons are intended to accompany the per sample H&E images. Any annotation which was a poor example of the feature, did not contain any cells, or covered more than 95% of the cells in a sample was removed from both the csv and the geojsons. 


---DIRECTIONS FOR SAMPLE-LEVEL TRANSCRIPTS.CSV FILES---
- Each "transcript_id" and "cell_id" has the sample name added to it as a prefix, which makes these files different from the original Xenium outputs. Old sample names are included for posterity, for people who previously accessed these data before publication and want to match up old sample names to the new, corrected ones.
- Transcripts that were not assigned to a segmented cell or nucleus have "UNASSIGNED" in the "cell_id" slot.
- "overlaps_nucleus" is 1 if the transcript was in a segmented nucleus, and 0 if not.
- The column "qv" contains quality metrics for each transcript. For most purposes, transcripts with qv >= 20 should be retained, and all others should be filtered out.
- In the column "feature_name", some features start with "NegControlProbe", "NegControlCodeword", or "UnassignedCodeword". Again, for most purposes these should be filtered out, as theyare not considered real transcripts. 343 unique gene IDs should remain after filtering.
- Samples on TMAs 1-4 were run separately from TMA5, and have slightly different output formats. See paper for details. For this reason, samples on TMA5 have an additional column "codeword_index".

---IMPORTANT! CHANGES FROM PRIOR VERSIONS----
Prior to publication in Nature Genetics, the biorxiv preprint version of this data had some errors. We discovered that in the decoding step of Xenium processing, the file that decodes "codewords" from the Xenium platform was incorrectly assigning a number of codewords to an incorrect gene. 
For example, the gene MRC1 was originally (incorrectly) labeled as ACTA2, while ACTA2 was labeled as COL1A2, resulting in smooth muscle cells (SMCs) that appeared to express COL1A2 but not ACTA2, and myeloid populations that appeared to express ACTA2 but not MRC1.

This has now been corrected. See notes attached to the new Nature Genetics publication for info. If you are using a prior version of this data, please switch to the new upload for accuracy. However, users who downloaded a prior version of this data may be interested to know which gene names have been corrected to match their true identity. Thank you for your understanding and interest in our work!

Gene "swaps" were as follows; Any genes not listed are correct between all versions of the Seurat object/transcript files/etc.:
Original --> Corrected
ACTA2 --> MRC1
APLNR --> TGFB3
CCL18 --> GDF15
CCL21 --> HLA-DRA
CDKN2A --> FGF2
CEACAM6 --> ELN
COL1A2 --> ACTA2
COL4A3 --> LAMP3
CSPG4 --> APLNR
ELANE --> WT1
ELN --> PDGFRB
FGF2 --> PLIN2
FGF7 --> HES1
FN1 --> COL1A2
GDF15 --> SFRP4
HES1 --> TGFB2
HLA-DRA --> CCL21
ITGA3 --> FGF7
ITGAX --> ITGB1
ITGB1 --> ITGA3
ITGB6 --> FN1
JCHAIN --> MUC5AC
LAMP3 --> CDKN2A
LYZ --> ITGB6
MRC1 --> COL4A3
MUC5AC --> PPARG
PDGFRA --> ITGAX
PDGFRB --> LYZ
PLIN2 --> CCL18
PPARG --> JCHAIN
SFRP4 --> CEACAM6
TGFB2 --> PDGFRA
TGFB3 --> CSPG4
TREM2 --> ELANE
WT1 --> TREM2
