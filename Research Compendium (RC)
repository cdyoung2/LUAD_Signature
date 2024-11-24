#   Young et al Research Compendium
##  Corey Young, Courtney D. Dill, Eric B. Dammer
### January 17, 2023

____________________________________________________________

OVERVIEW

This package is a "Research Compendium" of analyses performed to compare and validate published, known, lung adenocarcinoma (LUAD) transcriptomic signatures versus the 8 gene survival predictor defined in this study by Young C, et al (2023).
Contents of this package for reproducing analysis which compares our LUAD survival gene signature described in our publication as described below, to known LUAD survival gene signatures.

The analysis uses TCGA STAR aligner software output gene-level unstranded read counts and unstranded FPKM data downloaded from the GDC portal (https://gdc.cancer.gov/) on 01-15-2023.
Clinical, exposure, and biospecimen-specific metadata for these cases were also downloaded from GDC.

More information about the TCGA mRNA expression pipeline summarizing mRNA gene-level counts, FPKM and other normalized abundance data can be found at https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm

____________________________________________________________

PACKAGE CONTENTS (OVERVIEW)

1) Curated Signature Gene Sets (gene symbol lists)
2) Survival Time Traits for TCGA LUAD patients
3) R code to reproduce the findings of our study, using GDC portal-downloaded TCGA LUAD gene expression quantitation, and items (1) and (2). 

____________________________________________________________

INSTRUCTIONS

To reproduce the analysis outlined in detail below, first extract the contents of this package on your Windows, Linux, or Mac. Then, follow steps I through III.

I. Download the GDC Data Transfer Tool for your platform (Mac, Windows, or Ubuntu/linux) from https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
   and run the gdc-client executable (January 2023 available version is v1.6.1) on the current manifest for TCGA-LUAD case samples (which can be downloaded from the above GDC portal link for filtered repository data).
   Then, store the gdc-client executable in the below subfolder (see step II.B.) of this research compendium.

II. DOWNLOAD DATA FOR LUAD PATIENTS IN THE TCGA COHORT BY EITHER:

   A. (Opting to download data for all available cohort patients at the time you are reproducing our analysis):
      The GDC portal provides public access to the above mentioned data, and the filters used to obtain a manifest of the data for the 557 TCGA cohort LUAD case samples (as of 01-15-2023) are applied in the following link, which can be pasted (as a single line) into a web browser:

https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22set_id%3ADMsMt4UB4aZtrit7NrNW%22%5D%7D%2C%22op%22%3A%22IN%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22STAR%20-%20Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%5D%7D

   - OR -

   B. (Choosing to use the provided copy of the GDC portal manifest as downloaded on September 19, 2024 to reproduce our LUAD survival gene signature comparative analysis):
      gdc_manifest.2024-09-19.txt is in the following subfolder of this research compendium:

      /2023wrapup/STAR Counts - Gene Expression Quant/

      Data in the above GDC manifest file is for 557 case samples, including 480 biological replicate (1 per patient case) primary LUAD tumor samples.
      557 subfolders will be created by the running the following command in the above subfolder with the provided manifest file.

   C. After either A OR B, run the following command from a console (Unix/Mac) or a command prompt (Windows), to download all data (specified in the manifest file from step A or B). The step B manifest is used below.

      gdc-client download -m gdc_manifest.2023-01-15.txt

III. RUN THE ANALYSIS

   A. Initiate an R console or RStudio session, and open in a text or code editor either of the following parallel pipeline .R scripts provided in this package:

      /2023wrapup/ROC+SignatureCorr-usingFPKM/FPKM_4signature_comparison_byCorr&ROC&Clustering-EBD.R
      /2023wrapup/ROC+SignatureCorr-usingCOUNTS/mRNA.Counts_4signature_comparison_byCorr&ROC&Clustering-EBD.R

   B. Edit the variable "rootdir" at the beginning of the script to specify the absolute path of this research compendium's extracted /2023wrapup/ folder.

   C. Run the .R script. Outputs you generate will write to /2023wrapup/. Compare these to the preprocessed pipeline output files in either of these subfolders:

      /2023wrapup/ROC+SignatureCorr-usingFPKM/
      /2023wrapup/ROC+SignatureCorr-usingCOUNTS/

____________________________________________________________

DETAILED DESCRIPTION OF PACKAGE CONTENTS

1) Survival Signature Gene Sets (symbol lists) from 3 previous studies were curated and information about how this was done is in the subfolder:
/2023wrapup/Curated Survival Signature Gene Sets/

The 3 studies for published gene signatures based on mRNA abundance that we curated to compare to ours follow:
	a) Shedden et al, Nat Med (2009), available online at https://www.nature.com/articles/nm.1790
	   459 genes covering 545 affymetrix array probes were already curated and available for download from http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6.html
	b) Soltis et al, Cell Rep Med (2022), available online at https://www.sciencedirect.com/science/article/pii/S2666379122003780
	   155 RNA-based metastasis-free survival (MFS) associated gene list [from Supplemental data Table 2C, available via the above link]
	c) Song et al, Sci Rep (2022), available online at https://www.nature.com/articles/s41598-022-14323-6.pdf
	   The positive and negative components of the signature in mRNA are taken from the hazard ratios of 11 inflammatory response genes making up the (IRG) signature shown in (Figure 2, panel D) of the above publication.


2) Survival time traits (OS, DSS, DFI, and PFI), and their corresponding binary values regarding censoring (vs. a corresponding relevant event occurring) at the time for all TCGA LUAD case samples were taken from the updated traits for >11,000 TCGA case samples as published in Liu, et al, Cell (2018).
   Liu, et al, Cell (2018) is available online at https://www.sciencedirect.com/science/article/pii/S0092867418302290

   Adapting and quoting from this reference's Notes on the previously published Supplemental Data: 
	   (i) Overall Survival (OS) - time in days to last_contact_days_to or death_days_to, whichever is larger.
	      Binary classification POSITIVE (value of 1) for death from any cause, NEGATIVE (value of 0, censored) for alive.
	   (ii) Disease Specific Survival (DSS) - time in days, last_contact_days_to or death_days_to, whichever is larger.
	      Binary classifier POSITIVE (value of 1) for patient whose vital_status was Dead and tumor_status was WITH TUMOR. If a patient died from the disease shown in field of cause_of_death, the status of DSS would be 1 for the patient.
	      Binary classifier NEGATIVE (value of 0, censored) for patient whose vital_status was Alive or whose vital_status was Dead and tumor_status was TUMOR FREE.
	      "This is not a 100% accurate definition but is the best we could do with this dataset. Technically a patient could be with tumor but died of a car accident and therefore incorrectly considered as an event."
	   (iii) Disease Free Interval (DFI) - time in days until new_tumor_event_dx_days_to (for POSITIVE events), or for censored cases, either last_contact_days_to or death_days_to, whichever is applicable.
	      Disease free status was defined by: first, treatment_outcome_first_course is "Complete Remission/Response"; if the tumor type doesn't have "treatment_outcome_first_course" then disease-free was defined by the value "R0" in the field of "residual_tumor";
	      otherwise, disease-free was defined by the value "negative" in the field of "margin_status". If the tumor type did not have any of these fields, then its DFI was NA.
	      Binary classifier POSITIVE (value of 1) for patient having new tumor event whether it is a local recurrence, distant metastasis, new primary tumor of the cancer, including cases with a new tumor event whose type is N/A.
	      Binary classifier NEGATIVE (value of 0, censored) for patients with disease free status during their followup period. 
	   (iv) Progression Free Interval (PFI) - time in days for a POSITIVE event regarding progression of disease in patient, or if negative (censored), the time recorded as either last_contact_days_to or death_days_to, whichever is applicable.
	      Progression is defined as patient having new tumor event whether it was a progression of disease, local recurrence, distant metastasis, new primary tumors (all sites), or died with the cancer without new tumor event, including cases with a new tumor event whose type is N/A.
	      For censored cases, recorded PFI time is either last_contact_days_to or death_days_to, whichever is applicable.
	      Binary classifier POSITIVE (value of 1) for the patient when trait metadata had a value for either new_tumor_event_dx_days_to or death_days_to, whichever is applicable.
	      Binary classifier NEGATIVE (value of 0, censored) for patients without progression during their followup period.

   NOTE: None of the 557 TCGA LUAD case samples were noted as having "Redaction" status.

   The above metadata for all TCGA LUAD cases as published by Liu, et al (2018), is in the following file of this package:
   /2023wrapup/Liu_et_al_2018_Endpoints_LUADonly.tsv


3) Following download of TCGA LUAD gene expression quantification data from the GDC portal as described above, the downloaded files can be loaded as input to R for further processing and visualization using  R code provided.
   The code is in the following .R files with additional detailed comments and instructions in this package:

   /2023wrapup/ROC+SignatureCorr-usingFPKM/FPKM_4signature_comparison_byCorr&ROC&Clustering-EBD.R
   /2023wrapup/ROC+SignatureCorr-usingCOUNTS/mRNA.Counts_4signature_comparison_byCorr&ROC&Clustering-EBD.R

   The R code is replicated to perform the same operations and visualization steps in each of the 2 above script files.
   The only difference between the 2 is whether [RNA-Seq unstranded mapped short read counts per gene] or [RNA-Seq FPKM gene-level normalized abundance] is used for calculation of survival gene siguatures.

____________________________________________________________

ANALYSIS PIPELINE STEP-BY-STEP DESCRIPTION

   Steps a through i, coded in the R script language, generate corresponding output files prefixed by 3x#, where x is a letter below representing the step in the pipeline, and # specifies sub-step corresponding to a small roman numeral.

	a) Import GDC portal-downloaded gene-level abundance (either counts or FPKM) into an R session.

	b) Assemble a corresponding matrix of sample-specific traits or metadata

	c) Visualize distributions of the full abundance data at different thresholds of filtering to remove rows with noise (low) abundance values, and compare to a prior gene subset with less noise.

	d) Calculate a prognostic ratio for each LUAD primary tumor sample biological replicate (N=480) for each of the 3 survival gene signatures above, and a fourth one defined in the current study using 8 gene product mRNA abundances, using ROC AUC optimization methodology.
	   This method was pioneered for obtaining a breast cancer survival biomarker in a previous study published in Dill CD, et al, iScience (2021), available at https://www.sciencedirect.com/science/article/pii/S2589004221004193
	   The full iterative optimization performed on TCGA LUAD data for this study identifying these 8 genes is available in a separate package.*

	e) Perform ROC curve analysis for each of the 4 signatures in the 480 primary tumor case sample biological replicates, determining prediction area under the curve (AUC) at 5 selected time points.
   	   Calculate related statistics for each signatures' ability to predict the below described 4 survival time and status traits:
	   (i) Overall Survival (OS)
	   (ii) Disease Specific Survival (DSS)
	   (iii) Disease Free Interval (DFI)
	   (iv) Progression Free Interval (PFI)

	   Prediction of status at a time point reflects known (or censored) outcome at that specific followup time after the tumor biopsy from which the quantification of mRNA occurred.
	   Every page of graphical output has each ROC curve for the signature's prediction of outcome at the following times:
	      12 months (green),
	      18 months (magenta),
	      3 years (red),
	      5 years (purple),
	      and 10 years (blue).

	   Output ROC curves (as a .PDF file) and a table (.CSV file) of all ROC statistics for each of the 4 above survival-related outcomes.

	f) Calculate gene signatures as a sum of abundance and signed log2(sum of abundance). In contrast to the prognostic ratio calculation, which weights each gene equally in the calculation, these calculations weight the genes by their relative abundance.
	   These alternate calculations of LUAD survival signatures are output in .CSV file tables in step h1 below, with the filename beginning "3h1."

	   NOTE:
	   This weighted calculation is similar to the gene expression signatures obtained for the signature gene lists entered into the USSC Xena Browser at https://xenabrowser.net/heatmap/ 
	   Briefly, on the Xena browser, select or type the following:
	 	Study: "GDC TCGA Lung Adenocarcinoma (LUAD)"
		First Variable: Genomic; Gene Expression; paste the signature gene list formatted as, e.g.:  "=ATP6V0E1+SVBP+HSDL1+UBTD1-GNPNAT1-XRCC2-TFAP2A-PPP1R13L"  (the 8 gene signature, from this study.)
		Second Variable: Genomic; Gene Expression; paste the signature gene list formatted as, e.g.:  "=ADM+NMI+PSEN1+PVR+SERPINE1+SHK1-GPC3-IL7R-NMUR1-PTPRE-SEMA4D"  (Song, et al, Sci Rep, 2022 signature as an example.)
		...additional genomic or phenotypic variables.

		Drag the box for the variable to sort by, to swap positions with the variable in position/box B. Find the formatted string for each of the 3 curated signatures in the above folder (as described in additional contents (1) above).

	g) Generate a set 3 of 8-page .PDF files, each containing all 4 signatures, calculated in 3 ways: (i) unweighted prognostic ratio, (ii) abundance sum, and (iii) signed log2(abundance sum). 
	   Each .PDF has 4 pages of centered (e.g. median-subtracted) and 4 pages of uncentered signature heatmaps, and each page includes all signature data for comparison, sorting the case samples (rows) of data for one of the 4 signatures (columns), from low to high.

	h) Compile (i) all traits (case sample metadata) used for the active analysis as a .CSV, including columns for the equal weight and summed abundance survival gene signatures.
	   Compile (ii) the GDC portal mRNA abundance data for all TCGA LUAD-related case samples as a .CSV file. No filtering of noise/missing values is applied.

	i) Calculate pairwise correlations of the survival gene signatures using both (i) Pearson (rho value) and (ii) bicor methods, and output each to a .CSV file.

____________________________________________________________
