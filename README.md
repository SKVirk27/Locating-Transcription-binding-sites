Transcription Binding Sites Identification
Description: This project converts count matrices into Position-Specific Scoring Matrices (PSSM) to locate transcription binding sites in a fasta file and identify the top 30 best hits.

Key Features:

Converts count matrices into PSSM.
Searches for binding sites in a given genome sequence.
Outputs include countmatrix.txt, scores.csv, and top30records.csv.
Required Libraries: os, pandas, numpy, pyfastx, Bio.motifs

Execution:

Ensure the required files (argR-counts-matrix.txt, E_coli_K12_MG1655.400_50.bz2, main.py) are in the same folder.
Open the project in PyCharm and set the working folder as Assignment4.
Run main.py to execute the script and identify the transcription binding sites.
