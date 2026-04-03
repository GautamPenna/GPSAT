# GPSAT - Genome/Protein Sequence Analysis Tool

GPSAT is a sequence analysis tool that biomedical reseachers can use to visualize positional variability and align sequences locally using MAFFT (Multiple Alignment using Fast Fourier Transform). This tool allows users to directly fetch sequences from UniProt and GenBank, load custom sequences, generate a consensus, identify conserved versus nonconserved regions and create beautiful figures that show positional variability.

> *Consensus Sequence: each nucleotide position is determined by the most common nucleotide across all data points*

### ✍️ Author: Gautam Penna, [Ke Lab](https://www.zunlongke-lab.org/), University of Texas at Austin

**PURPOSE**
When present with the challenge of analyzing genomes or proteomes, it becomes increasingly complicated when one realizes that no two sequences are exactly the same. Mutations between and within the same genotype, subgroup or serotype, are inevitable. The need for a consensus sequence is evident and the associated variability at each position is crucial for understanding the most mutable and constant regions of sequences. Understanding this functionality allows for the creation of antibody and future drug targets. This application takes its users through a step-by-step methodology of determining the consensus sequence and variability of the protein at each position. 

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
(5) [LogoMaker]

**Features**
Here is what you can do with GPSAT
(1) 
