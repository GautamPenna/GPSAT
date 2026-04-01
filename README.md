# GPSAT - Genome/Protein Sequence Analysis Tool

This respository contains code files and a step-by-step documentation on how to use VCSAT for Viral Genome Sequence Analysis. Genome downloads from NCBI can be downloaded in a variety of formats, but to parse through the data and determine a *consensus sequence* is difficult. This tool guides its users through a series of steps, starting from downloading the coding regions of requested proteins on NCBI to graphical analysis when compared against the consensus or reference sequences. 

> *Consensus Sequence: each nucleotide position is determined by the most common nucleotide across all data points*

### ✍️ Author: Gautam Penna, [Ke Lab](https://www.zunlongke-lab.org/), University of Texas at Austin

**PURPOSE**
When present with the challenge of analyzing a viral genome, it becomes increasingly complicated when one realizes that no two-genomes are exactly the same. Mutations between every viral particle, even of the same genotype, subgroup or serotype, are inevitable. It is been agreed upon that the difference between genotypes is a genomic difference of 8% or greater, with subgroup variations between 4-8%. However, even within subgroup, nucleotide and even protein differences arise. The need for a consensus sequence is and the associated variability at each position is crucial for understanding the most mutable and constant regions of a viral genome and protein. Understanding this functionality allows for the creation of antibody and future drug targets against the virus. This application takes its users through a step-by-step methodology of determining the consensus sequence and variability of the protein at each position. We hope that you will be able to utilize this in any means possible to help your project.

**Methodology Overview**:
1.	Downloading sequences from NCBI
2.	Determining Protein Types and Selecting Key Protein Types
3.	Determining Longest Protein for every Accession Number (ID)
4.	Creating Consensus Determination Files and MAFFT Alignment
5.	Variability Calculation and Graphical Representation
6.	ORF Match File Creation and ORF Analysis
7.	Ascension Number Deletion/Genome Modifier
8.  Protein Translations and Votes File
9.  Color Coded Graph
10. Consensus Sequence Printer

## General Outline of Sequence Determination

The following steps determine the variability among the consensus sequence so it can be later used in AlphaFold and other 3d protein modeling softwares to determine structural variability of a virus. The first couple steps in this process involve downloading the required genomes from NCBI and determining the proteins required for analysis. The longest protein of each ID number (Accession Number) is used for further analysis since using incomplete DNA sequences will lead to multiple ‘-’ during the MAFFT alignment process. After creating the files required by the MAFFT software, variability between the determined consensus sequence and each individual is calculated. 

From here, it is important to analyze the ORF (open reading frame) of the start and stop codons. Finding ones that are in the same reading frame that also encompass the largest protein is crucial for further analysis. After determining the Accession Numbers that do not align with the ORF determined, deletion of these ID’s from the dataset and a run-through of steps 6-13 is required for accurate depiction of the protein.

**File Organization**

To view the code used for each step, reference the */code* directory in this repository. To view example files of each step, reference the */test_files* directory in the respository. 
The future_plans file present in the repository represents what still needs to be done in this project.

## Downloading Code for Usage

In order to use the code present /code file, download the respository to your personal computer. After doing so, open your terminal window and navigate to directory you want the code to be placed in. Once reached, type the following command: 

```

git clone https://github.com/GautamPenna/VCSAT

```

This should download the repistory to the directory you are in. Once downloaded, type *"`cd VCSAT`"* to go inside the repository. After doing so, type the following to start the command:

```py

python VCSAT.py

```

When updates are indicated, enter the VCSAT directory and type:

```

git pull

```

This insures that the latest updates to the code are downloaded to your computer without overlap.

Before running this code, make sure your computer has the numpy, pandas, matplotlib and seaborne python libraries installed. If not, follow the commands below to do so:

```

pip install numpy
pip install pandas
pip install matplotlib
pip install seaborne

```

To understand what each library is used for, please refer to the documentation below:
(1) [Numpy](https://numpy.org/doc/)
(2) [Pandas](https://pandas.pydata.org/docs/)
(3) [MatPlotLib](https://matplotlib.org/stable/index.html)
(4) [SeaBorne](https://seaborn.pydata.org/)

## (1) Downloading NCBI Viral Protein Genomes

This step guides the user through choosing the coding regions of the virus of interest. If already done, please move on to step 2.

Head over to the [NCBI Virus Database](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), and type the name of the virus of interest in the search bar.

After loading the page for the main databse, select the subset of datapoints based on the selection bar on the left. After doing so, select download on the top left-hand corner and follow the steps below: *"`Download -> Coding Regions -> All Records -> Default Naming Convention`"*. The output of this download is going to be a *"`.fasta`"* File. Open the file in *"`.txt`"* format and save that as a new file. For ease, name this *"`all_coding_regions.txt`"*. This is going to be the starting file for the rest of the functions.

## (2) Determining Protein Types Present and Selecting Key Protein Types

Run *"`python VCSAT.py`"* on your terminal window. When prompted to, type **01**. When prompted to input the specific file, type the entire pathway for the *"`all_coding_regions.txt`"* file. For example, if the file is present on your desktop **AND** if using a MacBook, the pathway would be *"`/Users/(username)/Desktop/all_coding_regions.txt`"*

The output of this command is going to be a list of all the various naming conventions for protein types that are present in the FASTA download. Copy and paste this list in a different document and pick out the protein types that you wish to include for further analysis. A file named ‘All_protein_entry.txt’ will be made. **DO NOT DELETE** this file if you wish to proceed with the code. 

For example, some the naming conventions present for the large surface protein for Hepatitis-B are *large S protein, truncated large surface protein, Large S protein, preS1/preS2/S protein, large S antigen, L-HBsAg*. Though many more are present, these are a few used. Make sure to select all possible naming conventions despite spelling errors, since this is what is used throughout the program. Make sure to list the selected protein types to the side before continuing the program.

Once the list is finalized, select **02** on the code. When prompted to, list out all the protein-type names to include. Make sure to paste them in the same format as they are present in the list. If spaces and other miscellaneous characters are present in the list, make sure to type it out in the exact same manner. Once finished, type done. The output file will be focused on just the protein types you chose.

An example of such a file can be found in *"`/test_files/(1)all_coding_regions/selectedpros.txt`"*.

Refer `/code/functions.py` and the `determining_protein_types` and `key_protein_lengths` function for technical code.

## (3) Determining Longest Protein Sequences

Picking the longest protein for each Accession ID is crucial since that is most representative of the entire protein length. For example, when dealing with the Large Surface Antigen of the Hep-B Virus (length = 1200 bp), entries close to this length will be chosen versus those with a larger difference. 

Run VCSAT.py code on your terminal window. When prompted to, pick option **03** and go along with the code. A file for each genotype will be created in your desired folder.

An example of such a file can be found in *"`/test_files/(2)pre_deletion_files/A_genotype_max_lengths.txt`"*.

## (4) Creating Files for Consensus Sequence Determination

After determining the longest genomic sequence of the protein to use for MAFFT alignment, it is time to determine the consensus sequence. We will be using the [MAFFT Alignment Website](https://mafft.cbrc.jp/alignment/server/index.html) for this feature. Before we do so, formatting the files in the required manner is necessary. Run *"`python VCSAT.py`"* on your terminal window. When prompted to, type **04**. Go along with the code, the output being files for the MAFFT program. The consensus determination files will be created in your desired directory. Make sure to run this code one time for each genotype. Open the link above and follow the instructions, using the file created as the input file. 

The bottom left pictures represent what should be selected with circles around the important portions in the menu options.

<img width="700" alt="image" src="https://github.com/user-attachments/assets/d0aaad7e-54d7-4e28-b6d8-db0d38c484d3" />

<img width="700" alt="image" src="https://github.com/user-attachments/assets/48c020d0-57ac-422a-be9f-1f45c362ec4a" />

Select Export Alignment in FASTA → Convert FASTA to .txt file. To be easy to access, it is recommended to name this in the following convention: `(genotype)_consensus_file.txt`. 

An example of such a file can be found in *"`/test_files/(2)pre_deletion_files/A_Determination.txt and ACon1.txt`"*.

## (5) Votes Calculation and Visualization

In order to determine the variability of each position, we will use the post-MAFFT file in order to keep the lengths of all genomes the same. Run *"`python VCSAT.py`"* on your terminal window. When prompted to, type **06**. This will create a new file with percent variability and the different numbers of A, T, G and C for each position. In order to visualize the variability, run VCSAT.py code on your terminal window again. When prompted to, type **07**. This should create a file that represents the variability of the consensus sequence at each position of the nucleotide.

An example of such a file can be found in *`/test_files/(2)pre_deletion_files/Avotes.txt`*.

An example of such a graphical representation can be found in *`/test_files/(4)graphs_pictures/A_pre_mod_graph.png`*.

## (6) ORF Analysis

After determining the variability of the entire consensus genome, it is important to find an ORF (Open Reading Frame) that includes both a start and stop codon. This truly represents a full protein sequence since it implies that the genetic code can be transcribed into a functional protein. In order to determine the most likely ORF, run *"`python VCSAT.py`"* on your terminal window. When prompted to, type **08**. Create the file in the desired directory and run VCSAT.py again. This time, pick option **(09)** and go along with the code. This allows you to analyze various options available for ORFs. When selecting likely candidate ORFs, make sure to account of the length of the desired protein and overall variability of the positions. More variability implies having a smaller sample size to perform analysis. 

An example of such a file can be found in *`/test_files/(2)pre_deletion_files/AORF.txt`*.

## (7) Post-ORF Data Cleaning

If at this stage of the code, there are likely Accession numbers that need to be deleted as well as genomes that could be modified due to the variability of the stop codon. To delete Accession numbers, run *"`python VCSAT.py`"* on your terminal window. When prompted to, type **05**. and go along with the code. The input for this code would be the max_lengths file that is the result at the end of Step 3. After deleting the said Accession numbers, run through steps 4-6. After getting the output of this step, run run *"`python VCSAT.py`"* on your terminal window. When prompted to, type **10**. and go along with the code. Edit the sequences are required and repeat Step 5.

An example of such a files can be found in *`/test_files/(3)pre_deletion_files`*.

## (8) Protein Translations and Protein Variability

After determining the full genome sequence of the desired proteins, it is important to translate the genome into a full amino acid sequence to determine the primary structure of the protein. Run *"`python VCSAT.py`"* on your terminal window. When prompted to, type **11**. This will prompt you to create to an alingment file that contains the protein translation of the DNA sequence. After doing so, type **12** when prompted to. This will allow you to create a votes file with variablities of the protein sequence. This is important since variabilities in the genome are not always present in the protein, espically when the variability is present in the 3rd nucleotide of a codon. 

## (9) Color-Coded Gradient Variability Graph

For a non-graphical visualization of the variability present in the genome and protein sequences, type **14** when prompted to do so. This will generate a gradient graph that displays the variability, with a dark-shade of red representing more variability.
 
