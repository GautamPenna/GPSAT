# GPSAT GUI — Genomic & Proteomic Sequence Analysis Toolkit

**A desktop application for fetching, aligning, visualizing, and comparing viral and genomic sequences.**

Built by Gautam Penna - KeLab.

---

## Rationale

Studying viral and genomic sequences usually means using multiple separate programs — and most of them require the command line. That works for software-savvy researchers, but it shuts out a large part of the biology community.

GPSAT GUI exists to fix that. It brings the most common sequence analysis tasks — searching databases, aligning sequences, and visualizing results — into one simple desktop app. No coding required.

It is also designed with **classrooms in mind**. Students can export their results as a `.gpsat` file and hand them in. Instructors can load all submissions at once and compare them side by side. This makes GPSAT useful as both a research tool and a teaching companion.

---

## Features

- **Fetch sequences** from NCBI GenBank or UniProt by search term or accession ID
- **Load local FASTA files** to build a sequence pool
- **Align sequences** using MAFFT (multiple sequence alignment)
- **Compute consensus sequences** with per-position variability metrics
- **Visualize alignments** in a color-coded interactive grid
- **Compare student sequences** against a master reference
- **Generate sequence logos** using logomaker
- **Export results** as PNG images or `.gpsat` JSON files

---

## Requirements

### System
- **Python 3.7 or higher**
- **MAFFT** — external alignment tool, must be installed separately and available on your system PATH

### Python Packages
Listed in `requirements.txt`:
```
biopython>=1.80
requests>=2.28
logomaker>=0.8
matplotlib>=3.6
Pillow>=9.0
pandas>=1.5
```

---

## Installation

**Step 1 — Clone the repository**
```bash
git clone <repository-url>
cd GPSAT_GUI
```

**Step 2 — Create a virtual environment (recommended)**
```bash
python -m venv venv

# macOS / Linux:
source venv/bin/activate

# Windows:
venv\Scripts\activate
```

**Step 3 — Install Python dependencies**
```bash
pip install -r requirements.txt
```

**Step 4 — Install MAFFT**

MAFFT must be installed on your system and accessible from your terminal.

- **macOS** (with Homebrew):
  ```bash
  brew install mafft
  ```
- **Linux** (Debian/Ubuntu):
  ```bash
  sudo apt install mafft
  ```
- **Windows**: Download the installer from https://mafft.cbrc.jp/alignment/software/ and add MAFFT to your system PATH.

To verify MAFFT is installed correctly, run:
```bash
mafft --version
```

---

## Running the Application

```bash
python main.py
```

The GUI will launch with two tabs at the top.

---

## Usage

### Tab 1 — Prepare Data

Use this tab to build and align a sequence pool.

1. **Fetch sequences** using the GenBank or UniProt search panel, or load a local FASTA file.
2. **Filter** the pool by protein type if needed.
3. Click **Run Pipeline** to align sequences with MAFFT and compute a consensus.
4. The output is a `.gpsat` file — a portable snapshot of the consensus and its variability data.

### Tab 2 — Alignment Viewer

This tab has three sub-views, selectable at the top:

**Alignment**
- Load an aligned FASTA file to view the full alignment.
- Each column is color-coded by conservation level (orange = variable, white = conserved).
- Optionally load a UniProt reference sequence to include alongside the alignment.
- Export the view as a PNG or as a `.gpsat` consensus file.

**Compare**
- Load a master `.gpsat` file (e.g., the instructor's reference).
- Load one or more student `.gpsat` files.
- The tool aligns all sequences together and highlights insertions and positional differences.
- Useful for grading lab assignments or reviewing sequence variation across samples.

**Custom**
- Manually build a custom sequence pool from FASTA or database sources.
- Export it as a new `.gpsat` file for later use.

---

## File Formats

### `.gpsat` files

GPSAT uses a custom `.gpsat` format for storing consensus sequence data. These are plain JSON files.

```json
{
  "type": "gpsat_consensus",
  "version": 1,
  "name": "sequence_name",
  "consensus": "ACDEFG...",
  "length": 123,
  "positions": [
    {
      "position": 0,
      "consensus": "A",
      "variability": 0.12,
      "gap_count": 1,
      "total": 10,
      "counts": {"A": 8, "C": 2}
    }
  ]
}
```

- **`consensus`** — the most common residue at each position
- **`variability`** — a score from 0 (fully conserved) to 1 (maximally variable)
- **`counts`** — raw residue counts at each position

`.gpsat` files are the primary exchange format between the Prepare Data and Alignment Viewer tabs.

**Note:** Legacy `.vgat` files from earlier versions of the toolkit are still supported and will load correctly.

---

## Built With

| Tool | Purpose |
|------|---------|
| [Tkinter](https://docs.python.org/3/library/tkinter.html) | Desktop GUI framework |
| [Biopython](https://biopython.org/) | GenBank queries and sequence parsing |
| [MAFFT](https://mafft.cbrc.jp/) | Multiple sequence alignment |
| [logomaker](https://logomaker.readthedocs.io/) | Sequence logo generation |
| [matplotlib](https://matplotlib.org/) | Plotting backend |
| [Pillow](https://python-pillow.org/) | Image export |
| [UniProt REST API](https://www.uniprot.org/help/api) | Reference sequence lookup |
| [NCBI Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25499/) | GenBank sequence retrieval |

---

## Project Structure

```
GPSAT_GUI/
├── main.py            # Application entry point; initializes the GUI window and tabs
├── functions2o.py     # All business logic, UI components, and data processing
└── requirements.txt   # Python package dependencies
```

---

## Notes

- An active internet connection is required for GenBank and UniProt lookups.
- MAFFT must be installed and on your PATH for alignment features to work. The app will show an error if it cannot find MAFFT.
- Font size in the alignment viewer can be adjusted using the slider in the toolbar.
