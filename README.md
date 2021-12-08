# Meylan et al. 2021 

This repository host the code used in our paper titled "Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer." currently in revision.

## Abstract 

Tertiary Lymphoid Structures and B cells positively impact cancer clinical outcome and response to immunotherapy. Using spatial transcriptomics and immunochemistry, we show that B cell maturation and plasma cell formation take place in TLS. In TLS+ tumors, IgG and IgA producing PCs disseminate into the tumor beds along fibroblastic tracks. B cell repertoire analysis revealed clonal expansion, diversification and selection in TLS and the presence of fully mature clonotypes at distance. We observed tumor cell-bound antibodies and demonstrated that TLS+ tumors exhibited not only high numbers of IgG-producing PCs but also high numbers of IgG-stained and apoptotic malignant cells, which suggests anti-tumor effector activity of these antibodies. Finally, therapeutic responses and progression-free survival correlated with IgG stained tumor cells in patients treated with immune check-point inhibitors. Altogether, these data demonstrate intra-tumoral generation of B cell immunity and antibody production that impact responses to immunotherapy.

![graphical abstract](https://user-images.githubusercontent.com/33417707/145216143-e9246525-c351-42ff-afd9-859e205836c0.png)

### Pre-processing of the Visium data 

Script detailing the pre-processing of Visium data are identified with the tag "00" and take as input the space ranger output (10X genomics).

### R data analysis  

R scripts used to generate the figures 1 and 2 are identified with their corresponding tag "01" ans "02".

### MiXCR clonotype analysis 

Shell scripts used to perform the MiXCR analysis were run on both bulk or spatial data and take as argument: (1) the input directory (2) the outpour directory (3) the mixcr executable (4) file names. 
