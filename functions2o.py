import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import os
import re
import json
import subprocess
import tempfile
import threading
import io

import pandas as pd
import requests

try:
    from Bio import Entrez, SeqIO
    _BIO_AVAILABLE = True
except ImportError:
    _BIO_AVAILABLE = False

try:
    from PIL import Image as _PILImage, ImageDraw as _PILDraw, ImageFont as _PILFont
    _PIL_AVAILABLE = True
except ImportError:
    _PIL_AVAILABLE = False

_LOGO_IMPORT_ERROR = None
try:
    import logomaker
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    _LOGO_AVAILABLE = True
except Exception as _e:
    _LOGO_AVAILABLE = False
    _LOGO_IMPORT_ERROR = str(_e)

# ─── SHARED THEME COLOURS (UT Austin) ────────────────────────────────────────
BG       = '#FFFFFF'
PANEL    = '#FDF6EE'
ENTRY_BG = '#FFF8F2'
FG       = '#1A1A1A'
DIM_FG   = '#8A6040'
ACCENT   = '#BF5700'
SEL_FG   = '#FFFFFF'
BORDER   = '#E0C9B0'
INS_BG   = '#D0E8FF'
INS_FG   = '#1A3A6A'
UNI_BG   = '#E8F4E8'   # light green — UniProt reference row
UNI_FG   = '#1A4A1A'   # dark green text in UniProt row

_GAP_EXCLUSION_THRESHOLD = 0.30
_UNI_SENTINEL            = '__GPSAT_UNIPROT__'  # internal MAFFT name for UniProt row

# ─── FONT SCALING ─────────────────────────────────────────────────────────────
_font_size = 10
_scalable_widgets = []


def _fs(base):
    return max(6, int(base * _font_size / 10))


def _register_widget_font(widget, family, base_size, *mods):
    _scalable_widgets.append((widget, family, base_size, mods))


def apply_font_scale():
    for widget, family, base, mods in _scalable_widgets:
        try:
            widget.config(font=(family, _fs(base)) + mods)
        except Exception:
            pass


def _hex_to_rgb(hex_color: str) -> tuple:
    h = hex_color.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


# ─── CORE PARSING ────────────────────────────────────────────────────────────

def _parse_header(header_line):
    header = header_line.lstrip('>').strip()
    if re.match(r'^(sp|tr|gb|ref|emb|dbj|pir|prf|pdb)\|', header):
        parts  = header.split('|', 2)
        acc_id = parts[1] if len(parts) > 1 else header.split()[0]
        rest   = parts[2].strip() if len(parts) > 2 else ''
        recname = re.search(r'RecName:\s*Full=([^;]+)', rest)
        if recname:
            protein_name = recname.group(1).strip()
        else:
            protein_name = rest.split()[0] if rest else acc_id
        os_match = re.search(r'OS=([^=]+?)(?:\s+OX=|\s*$)', rest)
        organism = os_match.group(1).strip() if os_match else ''
        return acc_id, protein_name, organism

    space_idx = header.find(' ')
    if space_idx == -1:
        return header, header, ''
    acc_id = header[:space_idx]
    rest   = header[space_idx + 1:].strip()
    if rest.startswith('|'):
        rest = rest[1:].strip()
    org_match = re.search(r'\[([^\[\]]+)\]\s*$', rest)
    if org_match:
        organism     = org_match.group(1).strip()
        protein_name = rest[:org_match.start()].strip().rstrip(',').strip()
    else:
        organism     = ''
        protein_name = rest
    return acc_id, protein_name, organism


def parse_fasta_entries(input_path=None, fasta_text=None):
    """Parse from a file path or raw FASTA text string.
    Returns (entries_list, sorted_protein_types_list)."""
    entries       = []
    protein_types = set()
    current_header = None
    seq_parts      = []

    if fasta_text is not None:
        lines = fasta_text.splitlines()
    else:
        with open(input_path, 'r') as fh:
            lines = fh.readlines()

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith('>'):
            if current_header is not None:
                _commit_entry(current_header, seq_parts, entries, protein_types)
            current_header = stripped
            seq_parts      = []
        else:
            seq_parts.append(stripped)
    if current_header is not None and seq_parts:
        _commit_entry(current_header, seq_parts, entries, protein_types)
    return entries, sorted(protein_types)


def _commit_entry(header, seq_parts, entries, protein_types):
    sequence     = ''.join(seq_parts).upper()
    acc_id, protein_name, organism = _parse_header(header)
    entries.append({'id': acc_id, 'protein_name': protein_name,
                    'sequence': sequence, 'organism': organism,
                    'raw_header': header})
    protein_types.add(protein_name)


# ─── GENBANK / UNIPROT FETCHING ───────────────────────────────────────────────

_ENTREZ_EMAIL = 'gpsat@kelab.local'  # required by NCBI; update as needed


def _fetch_genbank(organism, gene, seq_type, max_results=200):
    """Fetch sequences from NCBI GenBank.
    seq_type: 'protein' or 'nucleotide'
    Returns list of dicts with keys: id, protein_name, sequence, organism, raw_header
    Raises RuntimeError on failure."""
    if not _BIO_AVAILABLE:
        raise RuntimeError(
            'Biopython is not installed.\n'
            'Install it with:  pip install biopython')
    Entrez.email = _ENTREZ_EMAIL
    db = 'protein' if seq_type == 'protein' else 'nuccore'
    query = f'{organism}[Organism] AND {gene}[Gene Name]'
    try:
        handle = Entrez.esearch(db=db, term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise RuntimeError(f'NCBI search failed:\n{e}')
    id_list = record.get('IdList', [])
    if not id_list:
        return []
    try:
        handle = Entrez.efetch(db=db, id=','.join(id_list),
                               rettype='fasta', retmode='text')
        fasta_text = handle.read()
        handle.close()
    except Exception as e:
        raise RuntimeError(f'NCBI fetch failed:\n{e}')
    entries, _ = parse_fasta_entries(fasta_text=fasta_text)
    return entries


def _fetch_uniprot_search(organism, gene, seq_type, max_results=200):
    """Fetch sequences from UniProt REST API.
    seq_type: 'protein' or 'nucleotide' (UniProt only has proteins)
    Returns list of dicts with keys: id, protein_name, sequence, organism, raw_header"""
    query_parts = []
    if organism:
        # Quote multi-word values so the API treats them as a phrase
        query_parts.append(f'organism_name:"{organism}"')
    if gene:
        query_parts.append(f'gene:"{gene}"')
    query = ' AND '.join(query_parts) if query_parts else '*'
    try:
        # Use params= so requests handles URL encoding correctly
        resp = requests.get(
            'https://rest.uniprot.org/uniprotkb/search',
            params={'query': query, 'format': 'fasta', 'size': max_results},
            headers={'User-Agent': 'GPSAT/1.0 (research tool)'},
            timeout=30)
        resp.raise_for_status()
    except Exception as e:
        raise RuntimeError(f'UniProt search failed:\n{e}')
    fasta_text = resp.text
    if not fasta_text.strip():
        return []
    if not fasta_text.strip().startswith('>'):
        raise RuntimeError(
            f'UniProt returned unexpected content (check organism/gene spelling):\n'
            f'{fasta_text[:300]}')
    entries, _ = parse_fasta_entries(fasta_text=fasta_text)
    return entries


def _fetch_uniprot_accession(query):
    """Fetch a single sequence from UniProt.
    Accepts: accession ID (P01871), entry name (IGHM_HUMAN), or gene name (IGHM).
    Returns (name, sequence_string) or raises RuntimeError."""
    query = query.strip()

    # Step 1: try direct accession / entry-name lookup
    url = f'https://rest.uniprot.org/uniprotkb/{query}.fasta'
    try:
        resp = requests.get(url, timeout=30,
                            headers={'User-Agent': 'GPSAT/1.0 (research tool)'})
        if resp.status_code == 200 and resp.text.strip().startswith('>'):
            entries, _ = parse_fasta_entries(fasta_text=resp.text.strip())
            if entries:
                e = entries[0]
                return e.get('protein_name') or query, e['sequence']
        if resp.status_code == 404:
            pass  # fall through to search
        elif resp.status_code != 400:
            # Unexpected HTTP error
            resp.raise_for_status()
    except RuntimeError:
        raise
    except Exception as e:
        raise RuntimeError(
            f'Network error contacting UniProt:\n{e}\n\n'
            f'Check your internet connection.')

    # Step 2: direct lookup failed (gene name / free text) — try search
    for field in (f'gene:"{query}"', f'protein_name:"{query}"'):
        try:
            r2 = requests.get(
                'https://rest.uniprot.org/uniprotkb/search',
                params={'query': field, 'format': 'fasta', 'size': 1},
                timeout=30)
            if r2.status_code == 200 and r2.text.strip().startswith('>'):
                entries, _ = parse_fasta_entries(fasta_text=r2.text.strip())
                if entries:
                    e = entries[0]
                    return e.get('protein_name') or query, e['sequence']
        except Exception:
            pass

    raise RuntimeError(
        f'Could not find "{query}" on UniProt.\n\n'
        f'Try entering:\n'
        f'  • An accession ID   e.g. P01871\n'
        f'  • An entry name     e.g. IGHM_HUMAN\n'
        f'  • A gene symbol     e.g. IGHM')


# ─── VOTES / CONSENSUS ───────────────────────────────────────────────────────

_NUCLEOTIDES = list('ATGCU')
_AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')


def _detect_seq_type(genomes):
    nucl = set('ATGCUN-')
    total = 0; nucl_count = 0
    for seq in genomes[:20]:
        sample = seq[:300].replace('-', '').replace('N', '').replace('U', '')
        total      += len(sample)
        nucl_count += sum(1 for c in sample if c in nucl)
    if total == 0:
        return 'nucleotide'
    return 'nucleotide' if (nucl_count / total) > 0.85 else 'protein'


def _votes_in_memory(aligned_seqs):
    """Compute consensus, variability, and per-position counts from an in-memory
    list of (name, aligned_seq) tuples.  Returns (consensus_str, variability_list,
    counts_data_list)."""
    genomes = [seq for _, seq in aligned_seqs]
    if not genomes:
        return '', [], []
    length      = min(len(s) for s in genomes)
    seq_type    = _detect_seq_type(genomes)
    tracked     = _NUCLEOTIDES if seq_type == 'nucleotide' else _AMINO_ACIDS
    tracked_set = set(tracked)
    all_counts  = []
    for i in range(length):
        counts = {c: 0 for c in tracked}
        counts['Gap'] = 0; counts['Other'] = 0
        for seq in genomes:
            ch = seq[i] if i < len(seq) else None
            if ch in ('-', 'N', 'X', None, ''):
                counts['Gap'] += 1
            elif ch in tracked_set:
                counts[ch] += 1
            else:
                counts['Other'] += 1
        all_counts.append(counts)
    raw_consensus = []
    for counts in all_counts:
        candidates = {k: v for k, v in counts.items()
                      if k not in ('Gap', 'Other') and v > 0}
        raw_consensus.append(max(candidates, key=candidates.get) if candidates else '-')
    final_consensus = []
    for i in range(length):
        con    = raw_consensus[i]
        counts = all_counts[i]
        if con == '-':
            valid   = {k: v for k, v in counts.items()
                       if k not in ('Gap', 'Other') and v > 0}
            flanked = (0 < i < length - 1 and
                       raw_consensus[i - 1] != '-' and raw_consensus[i + 1] != '-')
            if valid and flanked:
                con = max(valid, key=valid.get)
        final_consensus.append(con)
    consensus_str = ''.join(final_consensus)
    variability   = []
    counts_data   = []
    for i in range(length):
        counts    = all_counts[i]
        total     = sum(counts.values())
        gap_count = counts.get('Gap', 0)
        con       = final_consensus[i]
        var       = 0.0
        if total > 0:
            eff = (total - gap_count
                   if gap_count / total > _GAP_EXCLUSION_THRESHOLD else total)
            if eff > 0:
                var = (eff - counts.get(con, 0)) / eff * 100
        variability.append(var)
        counts_data.append({'position':    i,
                             'consensus':   con,
                             'variability': var,
                             'gap_count':   gap_count,
                             'total':       total,
                             'counts':      {k: v for k, v in counts.items()
                                             if k not in ('Gap', 'Other') and v > 0}})
    return consensus_str, variability, counts_data


def _remove_gap_consensus_cols(aligned_seqs, consensus, variability, counts_data):
    """Drop columns where consensus == '-'."""
    keep = [i for i, c in enumerate(consensus) if c != '-']
    if not keep:
        return [], '', [], []
    new_seqs = [(name, ''.join(seq[i] if i < len(seq) else '-' for i in keep))
                for name, seq in aligned_seqs]
    new_consensus    = ''.join(consensus[i]    for i in keep)
    new_variability  = [variability[i]         for i in keep]
    new_counts_data  = []
    for new_i, old_i in enumerate(keep):
        d = dict(counts_data[old_i])
        d['position'] = new_i
        new_counts_data.append(d)
    return new_seqs, new_consensus, new_variability, new_counts_data


def _build_master_numbering(master_aligned_seq):
    labels      = []
    master_pos  = 0
    ins_count   = 0
    for ch in master_aligned_seq:
        if ch != '-':
            master_pos += 1
            ins_count   = 0
            labels.append(str(master_pos))
        else:
            ins_count += 1
            labels.append(f'{master_pos}+{ins_count}')
    return labels


def _counts_data_to_df(counts_data):
    if not counts_data:
        return None
    rows = []
    for d in counts_data:
        row = {'Position':           d['position'],
               'Consensus':          d['consensus'],
               'Percent_Variability': d['variability'],
               'Gap_count':          d['gap_count'],
               'Total':              d['total']}
        for letter, count in d.get('counts', {}).items():
            row[f'{letter}_count'] = count
        rows.append(row)
    return pd.DataFrame(rows)


def _most_representative_sequence(seqs, consensus):
    """Return (name, similarity_fraction) for the sequence most similar to the
    consensus at non-gap consensus positions."""
    best_name  = None
    best_score = -1.0
    for name, seq in seqs:
        matches  = sum(1 for s, c in zip(seq, consensus)
                       if c not in ('-', 'N', 'X') and s == c)
        non_gap  = sum(1 for c in consensus if c not in ('-', 'N', 'X'))
        score    = matches / non_gap if non_gap > 0 else 0.0
        if score > best_score:
            best_score = score
            best_name  = name
    return best_name, best_score


def _cross_variability(aligned_seqs, length):
    variability = []; consensus = []
    for i in range(length):
        residues = [s[i] for _, s in aligned_seqs
                    if i < len(s) and s[i] not in ('-', 'N', 'X')]
        total  = len(aligned_seqs)
        n_gap  = total - len(residues)
        if not residues:
            variability.append(100.0); consensus.append('-'); continue
        effective = len(residues) if n_gap / total > _GAP_EXCLUSION_THRESHOLD else total
        con_char  = max(set(residues), key=residues.count)
        con_count = residues.count(con_char)
        var = (effective - con_count) / effective * 100.0 if effective > 0 else 0.0
        variability.append(var); consensus.append(con_char)
    return variability, ''.join(consensus)


# ─── Amino-acid property groups for conservation colouring ───────────────────

# Each frozenset is a group of amino acids that are biochemically similar.
# A column is "similar" if ALL non-gap residues belong to at least one group.
_PROP_GROUPS = (
    frozenset('AVILMFWP'),   # hydrophobic / nonpolar
    frozenset('FWYH'),        # aromatic (incl. His)
    frozenset('KRH'),         # basic / positive
    frozenset('DE'),          # acidic / negative
    frozenset('STNQCY'),     # polar uncharged
    frozenset('ST'),          # hydroxyl
    frozenset('NQ'),          # amide
    frozenset('AVILM'),       # aliphatic hydrophobic
    frozenset('ILV'),         # branched-chain aliphatic
)


def _classify_columns(all_rows, length):
    """Return a list of 'identical' | 'similar' | 'variable' for each column.

    'identical' — every sequence has the exact same residue (no gaps allowed).
    'similar'   — all non-gap residues belong to the same biochemical group.
    'variable'  — everything else.
    """
    classes = []
    total_count = len(all_rows)
    for col in range(length):
        aa = []
        has_gap  = False
        gap_count = 0
        for _, seq in all_rows:
            if col < len(seq):
                c = seq[col].upper()
                if c in ('-', '.'):
                    has_gap = True
                    gap_count += 1
                elif c not in ('X', 'N', 'B', 'Z'):
                    aa.append(c)
        # Majority gaps → treat as variable (no blue box)
        if gap_count > total_count / 2:
            classes.append('variable')
            continue
        aa_set = frozenset(aa)
        if not has_gap and len(aa_set) == 1:
            classes.append('identical')
        elif aa_set and any(aa_set <= g for g in _PROP_GROUPS):
            classes.append('similar')
        else:
            classes.append('variable')
    return classes


# ─── Compare file loading ────────────────────────────────────────────────────

def _load_cmp_seq(path):
    """Return the sequence string from a .gpsat/.vgat or FASTA file."""
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.gpsat', '.vgat'):
        with open(path) as fh:
            data = json.load(fh)
        if data.get('type') not in ('gpsat_consensus', 'vgat_consensus'):
            raise ValueError('Not a valid .gpsat consensus file.')
        return data['consensus']
    else:
        entries, _ = parse_fasta_entries(input_path=path)
        if not entries:
            raise ValueError('No sequences found in file.')
        return entries[0]['sequence']


# ─── MAFFT ───────────────────────────────────────────────────────────────────

def _run_mafft(sequences):
    fasta_lines = ''.join(f'>{name}\n{seq}\n' for name, seq in sequences)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta',
                                     delete=False) as tmp:
        tmp.write(fasta_lines)
        tmp_path = tmp.name
    try:
        result = subprocess.run(
            ['mafft', '--auto', '--quiet', tmp_path],
            capture_output=True, text=True, timeout=300)
    except FileNotFoundError:
        raise RuntimeError(
            'MAFFT executable not found.\n\n'
            'Install MAFFT and make sure it is on your PATH:\n'
            '  • macOS:  brew install mafft\n'
            '  • Linux:  sudo apt install mafft\n'
            '  • Windows: https://mafft.cbrc.jp/alignment/software/')
    except subprocess.TimeoutExpired:
        raise RuntimeError('MAFFT timed out after 300 s.')
    finally:
        os.unlink(tmp_path)
    if result.returncode != 0:
        raise RuntimeError(f'MAFFT returned an error:\n{result.stderr[:400]}')
    aligned = []; cur_name = None; cur_parts = []
    for line in result.stdout.splitlines():
        if line.startswith('>'):
            if cur_name is not None:
                aligned.append((cur_name, ''.join(cur_parts).upper()))
            cur_name = line[1:].strip(); cur_parts = []
        else:
            cur_parts.append(line.strip())
    if cur_name is not None:
        aligned.append((cur_name, ''.join(cur_parts).upper()))
    return aligned


# ─── PIL font loader ──────────────────────────────────────────────────────────

_MONO_FONT_PATHS = [
    '/System/Library/Fonts/SFNSMono.ttf',
    '/System/Library/Fonts/Supplemental/Courier New.ttf',
    '/System/Library/Fonts/Supplemental/Andale Mono.ttf',
    'C:/Windows/Fonts/cour.ttf',
    '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf',
]


def _load_mono_font(size):
    if not _PIL_AVAILABLE:
        return None
    for path in _MONO_FONT_PATHS:
        try:
            return _PILFont.truetype(path, size)
        except Exception:
            continue
    try:
        return _PILFont.load_default(size=size)
    except TypeError:
        return _PILFont.load_default()


# ─── TAB HELPERS ──────────────────────────────────────────────────────────────

def _styled_listbox(parent, **kw):
    lb = tk.Listbox(parent, bg=ENTRY_BG, fg=FG,
                    selectbackground=ACCENT, selectforeground=SEL_FG,
                    activestyle='none', borderwidth=0, highlightthickness=1,
                    highlightcolor=BORDER, highlightbackground=BORDER,
                    font=('Helvetica', _fs(10)), **kw)
    _register_widget_font(lb, 'Helvetica', 10)
    return lb


def _styled_text(parent, **kw):
    w = tk.Text(parent, bg=ENTRY_BG, fg=FG, insertbackground=FG,
                selectbackground=ACCENT, selectforeground=SEL_FG,
                borderwidth=0, highlightthickness=1,
                highlightcolor=BORDER, highlightbackground=BORDER,
                font=('Courier', _fs(10)), **kw)
    _register_widget_font(w, 'Courier', 10)
    return w


def _section_label(parent, text):
    return ttk.Label(parent, text=text, style='Section.TLabel')


# ─── GPSAT FILE HELPERS ───────────────────────────────────────────────────────

def _write_gpsat(path, name, consensus, variability, counts_data):
    """Save a gpsat_consensus JSON file."""
    positions = []
    for i, d in enumerate(counts_data):
        positions.append({
            'position':    i,
            'consensus':   d.get('consensus', consensus[i] if i < len(consensus) else '-'),
            'variability': d.get('variability', variability[i] if i < len(variability) else 100.0),
            'gap_count':   d.get('gap_count', 0),
            'total':       d.get('total', 0),
            'counts':      d.get('counts', {})
        })
    payload = {'type': 'gpsat_consensus', 'version': 1,
               'name': name, 'consensus': consensus,
               'length': len(consensus), 'positions': positions}
    with open(path, 'w') as fh:
        json.dump(payload, fh, indent=2)


def _write_gpsat_single_seq(path, name, seq_name, sequence):
    """Save a single (non-consensus) sequence as a gpsat_consensus file.
    Variability is set to 0 (fully conserved) since it's one sequence."""
    length = len(sequence.replace('-', ''))
    raw = sequence.replace('-', '')
    positions = [{'position': i, 'consensus': c,
                  'variability': 0.0, 'gap_count': 0, 'total': 1,
                  'counts': {c: 1}}
                 for i, c in enumerate(raw)]
    payload = {'type': 'gpsat_consensus', 'version': 1,
               'name': seq_name or name,
               'consensus': raw,
               'length': length, 'positions': positions}
    with open(path, 'w') as fh:
        json.dump(payload, fh, indent=2)


# ════════════════════════════════════════════════════════════════════════════════
#  TAB 1: PREPARE DATA
# ════════════════════════════════════════════════════════════════════════════════

class ProteinFilterTab(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, style='Panel.TFrame')
        self._pool = []          # list of entry dicts in the sequence pool
        self._fetch_results = [] # list of entry dicts from last fetch
        self._build()

    # ── Build UI ─────────────────────────────────────────────────────────────

    def _build(self):
        # ── Scrollable container so content stays accessible at any font size ──
        _scroll_host = tk.Frame(self, bg=PANEL)
        _scroll_host.pack(fill='both', expand=True)
        _vsb = ttk.Scrollbar(_scroll_host, orient='vertical')
        _vsb.pack(side='right', fill='y')
        _cv = tk.Canvas(_scroll_host, bg=PANEL, highlightthickness=0,
                        yscrollcommand=_vsb.set)
        _cv.pack(side='left', fill='both', expand=True)
        _vsb.config(command=_cv.yview)
        outer = ttk.Frame(_cv, style='Panel.TFrame')
        _win = _cv.create_window((14, 10), window=outer, anchor='nw')

        def _on_inner_resize(event):
            _cv.configure(scrollregion=_cv.bbox('all'))
        outer.bind('<Configure>', _on_inner_resize)

        def _on_canvas_resize(event):
            # Keep outer frame as wide as the canvas minus left+right margins
            _cv.itemconfig(_win, width=max(1, event.width - 28))
        _cv.bind('<Configure>', _on_canvas_resize)

        def _tab1_wheel(event):
            if event.num == 4:
                units = -1
            elif event.num == 5:
                units = 1
            elif event.delta:
                units = -1 if event.delta > 0 else 1
            else:
                return
            _cv.yview_scroll(units, 'units')
        for _w in (_cv, outer):
            _w.bind('<MouseWheel>', _tab1_wheel)
            _w.bind('<Button-4>',   _tab1_wheel)
            _w.bind('<Button-5>',   _tab1_wheel)

        # all content below is packed into `outer` (no padding needed — handled above)

        # ── FETCH SECTION ────────────────────────────────────────────────────
        _section_label(outer, 'Fetch Sequences from Database:').pack(
            anchor='w', pady=(0, 4))

        fetch_grid = ttk.Frame(outer, style='Panel.TFrame')
        fetch_grid.pack(fill='x', pady=(0, 4))

        ttk.Label(fetch_grid, text='Organism:', style='Panel.TLabel').grid(
            row=0, column=0, sticky='w', padx=(0, 4))
        self._org_var = tk.StringVar()
        ttk.Entry(fetch_grid, textvariable=self._org_var,
                  width=22).grid(row=0, column=1, sticky='ew', padx=(0, 12))

        ttk.Label(fetch_grid, text='Gene / Protein:', style='Panel.TLabel').grid(
            row=0, column=2, sticky='w', padx=(0, 4))
        self._gene_var = tk.StringVar()
        ttk.Entry(fetch_grid, textvariable=self._gene_var,
                  width=22).grid(row=0, column=3, sticky='ew')
        fetch_grid.columnconfigure(1, weight=1)
        fetch_grid.columnconfigure(3, weight=1)

        type_row = ttk.Frame(outer, style='Panel.TFrame')
        type_row.pack(fill='x', pady=(0, 6))
        ttk.Label(type_row, text='Type:', style='Panel.TLabel').pack(
            side='left', padx=(0, 8))
        self._seq_type_var = tk.StringVar(value='protein')
        ttk.Radiobutton(type_row, text='Protein',
                        variable=self._seq_type_var,
                        value='protein').pack(side='left', padx=(0, 12))
        ttk.Radiobutton(type_row, text='Nucleotide',
                        variable=self._seq_type_var,
                        value='nucleotide').pack(side='left')

        fetch_btn_row = ttk.Frame(outer, style='Panel.TFrame')
        fetch_btn_row.pack(fill='x', pady=(0, 4))
        ttk.Button(fetch_btn_row, text='Fetch GenBank',
                   command=self._fetch_genbank,
                   style='Accent.TButton').pack(side='left', padx=(0, 8))
        ttk.Button(fetch_btn_row, text='Fetch UniProt',
                   command=self._fetch_uniprot).pack(side='left', padx=(0, 12))
        self._fetch_status = ttk.Label(fetch_btn_row, text='', style='Dim.TLabel')
        self._fetch_status.pack(side='left')

        # Fetch results list
        ttk.Label(outer, text='Fetch results — check entries to add to pool:',
                  style='Dim.TLabel').pack(anchor='w', pady=(0, 2))
        res_frame = ttk.Frame(outer, style='Panel.TFrame')
        res_frame.pack(fill='x', pady=(0, 4))
        self._result_lb = _styled_listbox(res_frame,
                                          selectmode='multiple', height=5)
        res_sb = ttk.Scrollbar(res_frame, orient='vertical',
                               command=self._result_lb.yview)
        self._result_lb.config(yscrollcommand=res_sb.set)
        self._result_lb.pack(side='left', fill='both', expand=True)
        res_sb.pack(side='right', fill='y')

        res_btn = ttk.Frame(outer, style='Panel.TFrame')
        res_btn.pack(fill='x', pady=(0, 4))
        ttk.Button(res_btn, text='Select All Results',
                   command=self._select_all_results).pack(side='left', padx=(0, 6))
        ttk.Button(res_btn, text='Add Selected to Pool',
                   command=self._add_results_to_pool).pack(side='left')

        ttk.Separator(outer, orient='horizontal').pack(fill='x', pady=(4, 6))

        # ── LOCAL FILE SECTION ───────────────────────────────────────────────
        _section_label(outer, 'OR Load Local FASTA File:').pack(
            anchor='w', pady=(0, 4))
        file_row = ttk.Frame(outer, style='Panel.TFrame')
        file_row.pack(fill='x', pady=(0, 4))
        self._input_path = tk.StringVar()
        ttk.Entry(file_row, textvariable=self._input_path,
                  state='readonly').pack(side='left', fill='x', expand=True)
        ttk.Button(file_row, text='Browse',
                   command=self._browse_input).pack(side='left', padx=(6, 0))
        ttk.Button(outer, text='Load FASTA into Pool',
                   command=self._load_fasta,
                   style='Accent.TButton').pack(pady=(4, 4))

        ttk.Separator(outer, orient='horizontal').pack(fill='x', pady=(4, 6))

        # ── SEQUENCE POOL ─────────────────────────────────────────────────────
        pool_hdr = ttk.Frame(outer, style='Panel.TFrame')
        pool_hdr.pack(fill='x', pady=(0, 2))
        _section_label(pool_hdr, 'Sequence Pool:').pack(side='left')
        self._pool_count_lbl = ttk.Label(pool_hdr, text='0 sequences',
                                          style='Dim.TLabel')
        self._pool_count_lbl.pack(side='left', padx=(10, 0))
        ttk.Button(pool_hdr, text='Clear Pool',
                   command=self._clear_pool).pack(side='right')

        pool_frame = ttk.Frame(outer, style='Panel.TFrame')
        pool_frame.pack(fill='x', pady=(0, 6))
        self._pool_lb = _styled_listbox(pool_frame, height=4)
        pool_sb = ttk.Scrollbar(pool_frame, orient='vertical',
                                command=self._pool_lb.yview)
        self._pool_lb.config(yscrollcommand=pool_sb.set)
        self._pool_lb.pack(side='left', fill='both', expand=True)
        pool_sb.pack(side='right', fill='y')

        ttk.Separator(outer, orient='horizontal').pack(fill='x', pady=(4, 6))

        # ── PROTEIN TYPE FILTER ───────────────────────────────────────────────
        _section_label(outer, 'Filter by Protein Type:').pack(
            anchor='w', pady=(0, 2))
        ttk.Label(outer,
                  text='Hold Ctrl / Cmd to select multiple types. '
                       'Leave empty to use all.',
                  style='Dim.TLabel').pack(anchor='w', pady=(0, 4))

        lb_container = ttk.Frame(outer, style='Panel.TFrame')
        lb_container.pack(fill='x', pady=(0, 4))
        self.protein_lb = _styled_listbox(lb_container,
                                          selectmode='multiple', height=5)
        type_sb = ttk.Scrollbar(lb_container, orient='vertical',
                                command=self.protein_lb.yview)
        self.protein_lb.config(yscrollcommand=type_sb.set, state='disabled')
        self.protein_lb.pack(side='left', fill='both', expand=True)
        type_sb.pack(side='right', fill='y')

        type_btn_row = ttk.Frame(outer, style='Panel.TFrame')
        type_btn_row.pack(fill='x', pady=(0, 4))
        ttk.Button(type_btn_row, text='Select All',
                   command=self._select_all_types).pack(side='left', padx=(0, 6))
        ttk.Button(type_btn_row, text='Deselect All',
                   command=self._deselect_all_types).pack(side='left')
        self.count_label = ttk.Label(type_btn_row, text='', style='Dim.TLabel')
        self.count_label.pack(side='right')
        self.protein_lb.bind('<<ListboxSelect>>', self._update_type_count)

        ttk.Button(outer, text='Refresh Protein Types from Pool',
                   command=self._refresh_protein_types).pack(pady=(0, 6))

        ttk.Separator(outer, orient='horizontal').pack(fill='x', pady=(4, 6))

        # ── PIPELINE ──────────────────────────────────────────────────────────
        _section_label(outer,
                        'Run Pipeline  →  aligned.fasta + .gpsat').pack(
            anchor='w', pady=(0, 2))
        ttk.Label(outer,
                  text='Filters selected protein types → MAFFT alignment '
                       '(saves aligned.fasta) → consensus in memory → saves .gpsat',
                  style='Dim.TLabel', wraplength=580).pack(
            anchor='w', pady=(0, 6))

        dir_row = ttk.Frame(outer, style='Panel.TFrame')
        dir_row.pack(fill='x', pady=(0, 4))
        _section_label(dir_row, 'Output Directory:').pack(
            side='left', padx=(0, 6))
        self._pipe_dir = tk.StringVar()
        ttk.Entry(dir_row, textvariable=self._pipe_dir,
                  state='readonly').pack(side='left', fill='x', expand=True)
        ttk.Button(dir_row, text='Browse',
                   command=self._browse_pipeline_dir).pack(
            side='left', padx=(6, 0))

        pipe_btn_row = ttk.Frame(outer, style='Panel.TFrame')
        pipe_btn_row.pack(fill='x', pady=(6, 4))
        self._pipe_btn = ttk.Button(pipe_btn_row,
                                     text='Run Pipeline → .gpsat',
                                     command=self._run_pipeline,
                                     style='Accent.TButton')
        self._pipe_btn.pack(side='left')
        self._pipe_status = ttk.Label(pipe_btn_row, text='Ready',
                                       style='Dim.TLabel')
        self._pipe_status.pack(side='left', padx=(12, 0))

        log_frame = ttk.Frame(outer, style='Panel.TFrame')
        log_frame.pack(fill='both', expand=True, pady=(0, 4))
        self._pipe_log = _styled_text(log_frame, height=6, state='disabled',
                                       wrap='word')
        log_sb = ttk.Scrollbar(log_frame, orient='vertical',
                                command=self._pipe_log.yview)
        self._pipe_log.config(yscrollcommand=log_sb.set)
        self._pipe_log.pack(side='left', fill='both', expand=True)
        log_sb.pack(side='right', fill='y')

    # ── Fetch GenBank ─────────────────────────────────────────────────────────

    def _fetch_genbank(self):
        org  = self._org_var.get().strip()
        gene = self._gene_var.get().strip()
        if not org and not gene:
            messagebox.showerror('Error',
                                 'Enter an organism name and/or gene name.')
            return
        self._fetch_status.config(text='Fetching from GenBank…')
        self._result_lb.delete(0, tk.END)
        self._fetch_results.clear()
        seq_type = self._seq_type_var.get()

        def worker():
            try:
                entries = _fetch_genbank(org, gene, seq_type)
                self.after(0, lambda e=entries: self._fetch_done(e, 'GenBank'))
            except Exception as exc:
                self.after(0, lambda m=str(exc): self._fetch_error(m))

        threading.Thread(target=worker, daemon=True).start()

    def _fetch_uniprot(self):
        org  = self._org_var.get().strip()
        gene = self._gene_var.get().strip()
        if not org and not gene:
            messagebox.showerror('Error',
                                 'Enter an organism name and/or gene name.')
            return
        self._fetch_status.config(text='Fetching from UniProt…')
        self._result_lb.delete(0, tk.END)
        self._fetch_results.clear()
        seq_type = self._seq_type_var.get()

        def worker():
            try:
                entries = _fetch_uniprot_search(org, gene, seq_type)
                self.after(0, lambda e=entries: self._fetch_done(e, 'UniProt'))
            except Exception as exc:
                self.after(0, lambda m=str(exc): self._fetch_error(m))

        threading.Thread(target=worker, daemon=True).start()

    def _fetch_done(self, entries, source):
        self._fetch_results = entries
        self._result_lb.delete(0, tk.END)
        for e in entries:
            org  = e.get('organism', '')
            name = e.get('protein_name', e['id'])
            label = (f"{e['id']}  |  {name[:40]}  |  {org[:30]}  "
                     f"  ({len(e['sequence'])} aa/nt)")
            self._result_lb.insert(tk.END, label)
        self._fetch_status.config(
            text=f'{len(entries)} result(s) from {source}')

    def _fetch_error(self, msg):
        self._fetch_status.config(text='Error — see dialog')
        messagebox.showerror('Fetch Error', msg)

    def _select_all_results(self):
        self._result_lb.select_set(0, tk.END)

    def _add_results_to_pool(self):
        idxs = self._result_lb.curselection()
        if not idxs:
            messagebox.showinfo('Nothing selected',
                                'Select rows in the results list first.')
            return
        added = 0
        existing_ids = {e['id'] for e in self._pool}
        for i in idxs:
            if i < len(self._fetch_results):
                entry = self._fetch_results[i]
                if entry['id'] not in existing_ids:
                    self._pool.append(entry)
                    existing_ids.add(entry['id'])
                    added += 1
        self._refresh_pool_ui()
        self._refresh_protein_types()
        messagebox.showinfo('Added', f'{added} sequence(s) added to pool.')

    # ── Local file ───────────────────────────────────────────────────────────

    def _browse_input(self):
        path = filedialog.askopenfilename(
            title='Select FASTA / Genome File',
            filetypes=[('Sequence Files', '*.txt *.fasta *.fa *.faa *.fna'),
                       ('All Files', '*.*')])
        if path:
            self._input_path.set(path)

    def _load_fasta(self):
        path = self._input_path.get()
        if not path or not os.path.isfile(path):
            messagebox.showerror('Error', 'Please select a valid input file first.')
            return
        try:
            entries, _ = parse_fasta_entries(input_path=path)
        except Exception as e:
            messagebox.showerror('Parse Error', str(e)); return
        if not entries:
            messagebox.showinfo('No Sequences',
                                'No sequence entries found in the file.')
            return
        existing_ids = {e['id'] for e in self._pool}
        added = 0
        for entry in entries:
            if entry['id'] not in existing_ids:
                self._pool.append(entry)
                existing_ids.add(entry['id'])
                added += 1
        self._refresh_pool_ui()
        self._refresh_protein_types()
        messagebox.showinfo('Loaded',
                            f'{added} new sequences added to pool.\n'
                            f'Pool total: {len(self._pool)}')

    # ── Pool management ───────────────────────────────────────────────────────

    def _refresh_pool_ui(self):
        self._pool_lb.delete(0, tk.END)
        for e in self._pool:
            org  = e.get('organism', '')
            name = e.get('protein_name', e['id'])
            self._pool_lb.insert(
                tk.END,
                f"{e['id']}  |  {name[:35]}  |  {org[:25]}  "
                f"({len(e['sequence'])} nt/aa)")
        self._pool_count_lbl.config(
            text=f'{len(self._pool)} sequence(s)')

    def _clear_pool(self):
        if not self._pool:
            return
        if messagebox.askyesno('Clear Pool', 'Remove all sequences from the pool?'):
            self._pool.clear()
            self._refresh_pool_ui()
            self.protein_lb.config(state='normal')
            self.protein_lb.delete(0, tk.END)
            self.protein_lb.config(state='disabled')
            self.count_label.config(text='')

    # ── Protein types ─────────────────────────────────────────────────────────

    def _refresh_protein_types(self):
        types = sorted({e['protein_name'] for e in self._pool})
        self.protein_lb.config(state='normal')
        self.protein_lb.delete(0, tk.END)
        for t in types:
            self.protein_lb.insert(tk.END, t)
        self.count_label.config(
            text=f'0 / {len(types)} selected')

    def _select_all_types(self):
        self.protein_lb.select_set(0, tk.END); self._update_type_count()

    def _deselect_all_types(self):
        self.protein_lb.selection_clear(0, tk.END); self._update_type_count()

    def _update_type_count(self, _event=None):
        n = len(self.protein_lb.curselection())
        total = self.protein_lb.size()
        self.count_label.config(text=f'{n} / {total} selected')

    # ── Pipeline ─────────────────────────────────────────────────────────────

    def _browse_pipeline_dir(self):
        path = filedialog.askdirectory(title='Select Output Directory')
        if path: self._pipe_dir.set(path)

    def _pipeline_log(self, msg, clear=False):
        def _update():
            self._pipe_log.config(state='normal')
            if clear:
                self._pipe_log.delete(1.0, tk.END)
            self._pipe_log.insert(tk.END, msg + '\n')
            self._pipe_log.see(tk.END)
            self._pipe_log.config(state='disabled')
        self.after(0, _update)

    def _run_pipeline(self):
        if not self._pool:
            messagebox.showerror('No Data', 'Add sequences to the pool first.')
            return
        out_dir = self._pipe_dir.get().strip()
        if not out_dir:
            messagebox.showerror('No Directory', 'Choose an output directory first.')
            return

        # Determine which sequences to use
        selected_idx = self.protein_lb.curselection()
        if selected_idx:
            selected_types = {self.protein_lb.get(i) for i in selected_idx}
            filtered = [e for e in self._pool
                        if e['protein_name'] in selected_types]
        else:
            filtered = list(self._pool)

        if not filtered:
            messagebox.showerror('Empty Filter',
                                 'No sequences match the selected protein types.')
            return

        os.makedirs(out_dir, exist_ok=True)
        # Derive stem from gene name or fallback
        gene = self._gene_var.get().strip()
        raw_stem = gene if gene else 'output'
        stem = re.sub(r'[^\w\-]', '_', raw_stem).strip('_') or 'output'

        self._pipe_btn.config(state='disabled', text='Running…')
        self._pipe_status.config(text='Starting…')
        self._pipeline_log('', clear=True)
        self._pipeline_log(f'Output directory : {out_dir}')
        self._pipeline_log(f'File stem        : {stem}')
        self._pipeline_log(f'Sequences        : {len(filtered)}')
        self._pipeline_log('─' * 50)

        def worker():
            try:
                result = self._pipeline_worker(out_dir, filtered, stem)
                self.after(0, lambda r=result: self._pipeline_done(r))
            except Exception as exc:
                self.after(0, lambda m=str(exc): self._pipeline_error(m))

        threading.Thread(target=worker, daemon=True).start()

    def _pipeline_worker(self, out_dir, filtered, stem):
        def log(msg):    self._pipeline_log(msg)
        def status(msg): self.after(0, lambda m=msg: self._pipe_status.config(text=m))

        # Step 1: Build unaligned FASTA in memory
        status('Step 1/3 — preparing sequences…')
        seqs = [(f"{e['id']}|{e.get('protein_name','')}|{e.get('organism','')}",
                 e['sequence']) for e in filtered]
        log(f'[1/3] Prepared {len(seqs)} sequences')

        # Step 2: MAFFT alignment
        status('Step 2/3 — running MAFFT (may take a moment)…')
        log(f'[2/3] Running MAFFT on {len(seqs)} sequences…')
        aligned = _run_mafft(seqs)
        aligned_path = os.path.join(out_dir, f'{stem}_aligned.fasta')
        with open(aligned_path, 'w') as fh:
            for name, seq in aligned:
                fh.write(f'>{name}\n{seq}\n')
        log(f'[2/3] Aligned FASTA    → {os.path.basename(aligned_path)}'
            f'  ({len(aligned[0][1]) if aligned else 0} positions)')

        # Step 3: Votes in memory → .gpsat
        status('Step 3/3 — computing consensus & saving .gpsat…')
        consensus, variability, counts_data = _votes_in_memory(aligned)
        aligned2, consensus2, variability2, counts_data2 = _remove_gap_consensus_cols(
            aligned, consensus, variability, counts_data)

        gpsat_path = os.path.join(out_dir, f'{stem}_consensus.gpsat')
        _write_gpsat(gpsat_path, stem, consensus2, variability2, counts_data2)
        log(f'[3/3] Consensus .gpsat  → {os.path.basename(gpsat_path)}')

        return {'consensus': consensus2,
                'n_seqs': len(aligned),
                'aligned_path': aligned_path,
                'gpsat_path': gpsat_path,
                'stem': stem,
                'out_dir': out_dir}

    def _pipeline_done(self, result):
        self._pipe_btn.config(state='normal', text='Run Pipeline → .gpsat')
        self._pipe_status.config(text='Done ✓')
        consensus = result['consensus']
        preview   = (consensus[:80] + '…') if len(consensus) > 80 else consensus
        self._pipeline_log('─' * 50)
        self._pipeline_log(f'Consensus length : {len(consensus)} positions')
        self._pipeline_log(f'Consensus        : {preview}')
        self._pipeline_log('')
        self._pipeline_log(f'Files saved to:  {result["out_dir"]}')
        self._pipeline_log(f'  Aligned FASTA : {os.path.basename(result["aligned_path"])}')
        self._pipeline_log(f'  Consensus .gpsat : {os.path.basename(result["gpsat_path"])}')
        self._pipeline_log('')
        self._pipeline_log('→ Open Alignment Viewer and load the aligned FASTA to explore.')

    def _pipeline_error(self, msg):
        self._pipe_btn.config(state='normal', text='Run Pipeline → .gpsat')
        self._pipe_status.config(text='Error')
        self._pipeline_log('─' * 50)
        self._pipeline_log(f'ERROR: {msg}')


# ════════════════════════════════════════════════════════════════════════════════
#  TAB 2: ALIGNMENT VIEWER
# ════════════════════════════════════════════════════════════════════════════════

class AlignmentViewerTab(ttk.Frame):

    _B_CW      = 14
    _B_CH      = 16
    _B_LW      = 162
    _B_WRAP_CW = 20
    _B_WRAP_CH = 24
    _B_CMP_CW  = 14
    _B_CMP_CH  = 16
    _B_RULER_H = 20

    @property
    def CW(self):      return max(8,  int(self._B_CW      * _font_size / 10))
    @property
    def CH(self):      return max(10, int(self._B_CH      * _font_size / 10))
    @property
    def LW(self):      return max(80, int(self._B_LW      * _font_size / 10))
    @property
    def WRAP_CW(self): return max(12, int(self._B_WRAP_CW * _font_size / 10))
    @property
    def WRAP_CH(self): return max(14, int(self._B_WRAP_CH * _font_size / 10))
    @property
    def CMP_CW(self):  return max(8,  int(self._B_CMP_CW  * _font_size / 10))
    @property
    def CMP_CH(self):  return max(10, int(self._B_CMP_CH  * _font_size / 10))
    @property
    def RULER_H(self): return max(12, int(self._B_RULER_H * _font_size / 10))

    def __init__(self, master):
        super().__init__(master, style='Panel.TFrame')

        # ── Alignment / Wrapped data ──────────────────────────────────────────
        self._seqs_all             = []
        self._pre_uniprot_seqs_all = []   # stable snapshot before any UniProt run
        self._seqs_raw_cov    = []
        self._seqs            = []
        self._consensus       = ''
        self._variability     = []
        self._counts_data     = []
        self._counts_df       = None
        self._length          = 0
        self._orig_seqs_all   = []
        self._orig_consensus  = ''
        self._orig_variability = []
        self._orig_counts_data = []
        self._orig_counts_df  = None
        self._orig_length     = 0

        # UniProt reference
        self._uniprot_seq          = None    # raw (ungapped) sequence string
        self._uniprot_name         = None
        self._uniprot_aligned      = None    # aligned sequence (after MAFFT with pool)
        self._trim_start           = None    # inclusive column index (start trim)
        self._trim_end             = None    # exclusive column index (end trim)
        self._trim_to_uniprot_var  = tk.BooleanVar(value=True)
        self._match_uniprot_var    = tk.BooleanVar(value=False)

        # LogoMaker toggle state for alignment view
        self._logo_visible_align  = False
        self._logo_canvas_align   = None
        self._logo_fig_align      = None
        self._logo_start_var      = tk.StringVar(value='')
        self._logo_end_var        = tk.StringVar(value='')

        # ── Compare view data ─────────────────────────────────────────────────
        self._cmp_seqs         = []
        self._cmp_variability  = []
        self._cmp_consensus    = ''
        self._cmp_length       = 0
        self._cmp_master_name  = None
        self._cmp_master_nums  = []
        self._cmp_student_paths = []    # list of (display_name, path, user_name)
        self._logo_visible_cmp  = False
        self._logo_canvas_cmp   = None
        self._logo_fig_cmp      = None
        self._logo_cmp_start_var = tk.StringVar(value='')
        self._logo_cmp_end_var   = tk.StringVar(value='')

        # ── Custom / .gpsat factory data ───────────────────────────────────────
        self._custom_entries     = []   # (name, raw_seq)

        # ── Thread state ──────────────────────────────────────────────────────
        self._cov_job         = None
        self._mafft_busy      = False
        self._uniprot_pending = False   # retry incorporation after current MAFFT

        # Layout dicts
        self._wrap_layout = {}
        self._cmp_layout  = {}

        # View mode: 'alignment' or 'wrapped' within the Alignment sub-tab
        self._align_view_mode = 'alignment'

        self._build()

    # ════════════════════════════════════════════════════════════════════════
    #  BUILD
    # ════════════════════════════════════════════════════════════════════════

    def _build(self):
        # ── Sub-tab selector ─────────────────────────────────────────────────
        sel = ttk.Frame(self, style='Panel.TFrame')
        sel.pack(fill='x', padx=14, pady=(10, 2))

        lbl = tk.Label(sel, text='VIEW:', bg=PANEL, fg=DIM_FG,
                       font=('Helvetica', _fs(8), 'bold'))
        lbl.pack(side='left', padx=(0, 10))
        _register_widget_font(lbl, 'Helvetica', 8, 'bold')

        self._view_mode = tk.StringVar(value='alignment')
        self._view_btns = {}
        for value, label in [('alignment', 'Alignment'),
                              ('compare',   'Compare'),
                              ('custom',    'Custom Sequences')]:
            b = ttk.Button(sel, text=label,
                           command=lambda v=value: self._switch_view(v))
            b.pack(side='left', padx=(0, 4))
            self._view_btns[value] = b

        ttk.Separator(self, orient='horizontal').pack(fill='x', padx=14, pady=(4, 2))

        self._status = ttk.Label(
            self, text='Select a view above, then load your data.',
            style='Dim.TLabel')
        self._status.pack(anchor='w', padx=14, pady=(0, 2))

        # ── Per-view control panels ───────────────────────────────────────────
        self._ctrl_host = tk.Frame(self, bg=PANEL)
        self._ctrl_host.pack(fill='x', padx=14, pady=(2, 0))

        self._ctrl_panels = {
            'alignment': self._build_alignment_ctrl(self._ctrl_host),
            'compare':   self._build_compare_ctrl(self._ctrl_host),
            'custom':    self._build_custom_ctrl(self._ctrl_host),
        }

        # ── Canvas container ──────────────────────────────────────────────────
        self._canvas_host = tk.Frame(self, bg=BG)
        self._canvas_host.pack(fill='both', expand=True,
                               padx=14, pady=(4, 14))

        self._canvas_frames = {
            'alignment': self._build_alignment_canvas(self._canvas_host),
            'compare':   self._build_compare_canvas(self._canvas_host),
            'custom':    self._build_custom_canvas(self._canvas_host),
        }

        self._switch_view('alignment')

    # ── Alignment control panel ───────────────────────────────────────────────

    def _build_alignment_ctrl(self, parent):
        f = tk.Frame(parent, bg=PANEL)

        # Row 1: file loading
        row1 = ttk.Frame(f, style='Panel.TFrame')
        row1.pack(fill='x', pady=(6, 0))
        _section_label(row1, 'Aligned FASTA:').pack(side='left', padx=(0, 4))
        self.align_path = tk.StringVar()
        ttk.Entry(row1, textvariable=self.align_path,
                  state='readonly', width=30).pack(
            side='left', fill='x', expand=True, padx=(0, 4))
        ttk.Button(row1, text='Browse',
                   command=self._browse_align).pack(side='left', padx=(0, 8))
        ttk.Button(row1, text='Import .gpsat',
                   command=self._import_gpsat).pack(side='left')
        ttk.Label(f, text='Note: submit .gpsat files only',
                  style='Dim.TLabel').pack(anchor='w', padx=2, pady=(1, 0))

        # Row 2: UniProt reference
        row2 = ttk.Frame(f, style='Panel.TFrame')
        row2.pack(fill='x', pady=(6, 0))
        _section_label(row2, 'UniProt ref (accession):').pack(
            side='left', padx=(0, 4))
        self._uniprot_acc_var = tk.StringVar()
        _uni_entry = ttk.Entry(row2, textvariable=self._uniprot_acc_var, width=18)
        _uni_entry.pack(side='left', padx=(0, 4))
        _uni_entry.bind('<Return>', lambda e: self._load_uniprot_ref())
        ttk.Button(row2, text='Load UniProt',
                   command=self._load_uniprot_ref).pack(side='left', padx=(0, 8))
        self._uniprot_status = ttk.Label(row2, text='', style='Dim.TLabel')
        self._uniprot_status.pack(side='left')
        ttk.Button(row2, text='Remove UniProt',
                   command=self._remove_uniprot_ref).pack(side='left', padx=(8, 0))
        ttk.Checkbutton(
            row2, text='Trim to UniProt length',
            variable=self._trim_to_uniprot_var,
            style='Panel.TCheckbutton'
        ).pack(side='left', padx=(12, 0))
        ttk.Checkbutton(
            row2, text='No gaps in UniProt',
            variable=self._match_uniprot_var,
            command=self._on_match_uniprot_toggle,
            style='Panel.TCheckbutton'
        ).pack(side='left', padx=(8, 0))

        # Row 3: view toggle + sliders
        row3 = ttk.Frame(f, style='Panel.TFrame')
        row3.pack(fill='x', pady=(6, 0))
        _section_label(row3, 'View:').pack(side='left', padx=(0, 4))
        self._align_view_btn = ttk.Button(
            row3, text='⇄ Switch to Wrapped',
            command=self._toggle_align_view)
        self._align_view_btn.pack(side='left', padx=(0, 16))

        _section_label(row3, 'Coverage filter:').pack(side='left', padx=(0, 4))
        self._cov_var = tk.DoubleVar(value=0.0)
        self._cov_lbl = ttk.Label(row3, text='Off', style='Dim.TLabel', width=16)
        self._cov_lbl.pack(side='left', padx=(0, 4))
        ttk.Scale(row3, from_=0, to=100, orient='horizontal',
                  variable=self._cov_var, length=130,
                  command=self._on_cov_change).pack(side='left', padx=(0, 16))

        # (variability-N slider removed — use coverage filter instead)

        # Row 4: action buttons
        row4 = ttk.Frame(f, style='Panel.TFrame')
        row4.pack(pady=(8, 6))
        _section_label(row4, 'Start #:').pack(side='left', padx=(0, 4))
        self._start_aa_var = tk.StringVar(value='1')
        ttk.Entry(row4, textvariable=self._start_aa_var,
                  width=6).pack(side='left', padx=(0, 10))
        ttk.Button(row4, text='Render',
                   command=self._render,
                   style='Accent.TButton').pack(side='left', padx=(0, 6))
        ttk.Button(row4, text='Most Representative',
                   command=self._find_most_representative).pack(
            side='left', padx=(0, 6))
        self._logo_btn_align = ttk.Button(
            row4, text='Show Sequence Logo ▼',
            command=self._toggle_logo_align)
        self._logo_btn_align.pack(side='left', padx=(0, 6))
        ttk.Button(row4, text='Export SVG',
                   command=self._export_png_align).pack(side='left', padx=(0, 6))
        ttk.Button(row4, text='Export .gpsat ▼',
                   command=self._export_gpsat_dialog).pack(
            side='left', padx=(0, 6))

        # Logo control row — shown only when logo is visible
        self._logo_ctrl_frm_align = ttk.Frame(f, style='Panel.TFrame')
        # not packed initially; shown in _toggle_logo_align
        _section_label(self._logo_ctrl_frm_align, 'Logo positions:').pack(
            side='left', padx=(0, 4))
        ttk.Entry(self._logo_ctrl_frm_align, textvariable=self._logo_start_var,
                  width=6).pack(side='left', padx=(0, 2))
        ttk.Label(self._logo_ctrl_frm_align, text='to',
                  style='Panel.TLabel').pack(side='left', padx=(0, 2))
        ttk.Entry(self._logo_ctrl_frm_align, textvariable=self._logo_end_var,
                  width=6).pack(side='left', padx=(0, 6))
        ttk.Button(self._logo_ctrl_frm_align, text='Update',
                   command=self._refresh_logo_align).pack(side='left', padx=(0, 8))
        ttk.Button(self._logo_ctrl_frm_align, text='Export Logo PNG',
                   command=self._export_logo_align).pack(side='left')

        return f

    # ── Compare control panel ─────────────────────────────────────────────────

    def _build_compare_ctrl(self, parent):
        f = tk.Frame(parent, bg=PANEL)

        # Master file
        mrow = ttk.Frame(f, style='Panel.TFrame')
        mrow.pack(fill='x', pady=(6, 0))
        _section_label(mrow, 'Master Reference (.gpsat / .fasta):').pack(
            side='left', padx=(0, 6))
        self._cmp_master_display = ttk.Label(mrow, text='(none)',
                                              style='Dim.TLabel')
        self._cmp_master_display.pack(side='left', padx=(0, 8))
        ttk.Button(mrow, text='Browse Master',
                   command=self._browse_cmp_master).pack(side='left')

        # Student files listbox
        _section_label(f, 'Student Sequences (.gpsat / .fasta):').pack(
            anchor='w', pady=(6, 2))
        ttk.Label(f, text='Note: submit .gpsat files only',
                  style='Dim.TLabel').pack(anchor='w', padx=2, pady=(0, 2))
        sb_frame = ttk.Frame(f, style='Panel.TFrame')
        sb_frame.pack(fill='x')
        self._cmp_student_lb = _styled_listbox(sb_frame,
                                               selectmode='single', height=4)
        sb_sb = ttk.Scrollbar(sb_frame, orient='vertical',
                              command=self._cmp_student_lb.yview)
        self._cmp_student_lb.config(yscrollcommand=sb_sb.set)
        self._cmp_student_lb.pack(side='left', fill='both', expand=True)
        sb_sb.pack(side='right', fill='y')

        sbtns = ttk.Frame(f, style='Panel.TFrame')
        sbtns.pack(fill='x', pady=(4, 0))
        ttk.Button(sbtns, text='Add Student File',
                   command=self._add_cmp_student).pack(side='left', padx=(0, 6))
        ttk.Button(sbtns, text='Remove Selected',
                   command=self._remove_cmp_student).pack(side='left')

        # Action buttons
        btn = ttk.Frame(f, style='Panel.TFrame')
        btn.pack(pady=(8, 6))
        _section_label(btn, 'Start #:').pack(side='left', padx=(0, 4))
        self._cmp_start_aa_var = tk.StringVar(value='1')
        ttk.Entry(btn, textvariable=self._cmp_start_aa_var,
                  width=6).pack(side='left', padx=(0, 10))
        ttk.Button(btn, text='Compare & Render',
                   command=self._run_compare,
                   style='Accent.TButton').pack(side='left', padx=(0, 6))
        self._logo_btn_cmp = ttk.Button(
            btn, text='Show Sequence Logo ▼',
            command=self._toggle_logo_cmp)
        self._logo_btn_cmp.pack(side='left', padx=(0, 6))
        ttk.Button(btn, text='Export SVG',
                   command=self._export_png_compare).pack(
            side='left', padx=(0, 6))

        # Logo control row — shown only when logo is visible
        self._logo_ctrl_frm_cmp = ttk.Frame(f, style='Panel.TFrame')
        # not packed initially; shown in _toggle_logo_cmp
        _section_label(self._logo_ctrl_frm_cmp, 'Logo positions:').pack(
            side='left', padx=(0, 4))
        ttk.Entry(self._logo_ctrl_frm_cmp, textvariable=self._logo_cmp_start_var,
                  width=6).pack(side='left', padx=(0, 2))
        ttk.Label(self._logo_ctrl_frm_cmp, text='to',
                  style='Panel.TLabel').pack(side='left', padx=(0, 2))
        ttk.Entry(self._logo_ctrl_frm_cmp, textvariable=self._logo_cmp_end_var,
                  width=6).pack(side='left', padx=(0, 6))
        ttk.Button(self._logo_ctrl_frm_cmp, text='Update',
                   command=self._refresh_logo_cmp).pack(side='left', padx=(0, 8))
        ttk.Button(self._logo_ctrl_frm_cmp, text='Export Logo PNG',
                   command=self._export_logo_cmp).pack(side='left')

        return f

    # ── Custom control panel ──────────────────────────────────────────────────

    def _build_custom_ctrl(self, parent):
        f = tk.Frame(parent, bg=PANEL)

        ttk.Label(f,
                  text='Add sequences to the pool, align via MAFFT, then export '
                       'as a .gpsat file for use in the Compare view.',
                  style='Dim.TLabel', wraplength=600).pack(
            anchor='w', pady=(4, 6))

        # Name + paste
        add_row = ttk.Frame(f, style='Panel.TFrame')
        add_row.pack(fill='x', pady=(0, 2))
        _section_label(add_row, 'Name:').pack(side='left', padx=(0, 4))
        self._cust_name_var = tk.StringVar()
        ttk.Entry(add_row, textvariable=self._cust_name_var,
                  width=18).pack(side='left', padx=(0, 10))
        ttk.Button(add_row, text='Add Sequence',
                   command=self._custom_add_from_entry).pack(side='left')

        ttk.Label(f, text='Paste sequence below then click Add Sequence:',
                  style='Dim.TLabel').pack(anchor='w', pady=(2, 2))
        seq_frame = ttk.Frame(f, style='Panel.TFrame')
        seq_frame.pack(fill='x', pady=(0, 6))
        self._cust_paste_text = _styled_text(seq_frame, height=3, wrap='char')
        paste_sb = ttk.Scrollbar(seq_frame, orient='vertical',
                                  command=self._cust_paste_text.yview)
        self._cust_paste_text.config(yscrollcommand=paste_sb.set)
        self._cust_paste_text.pack(side='left', fill='both', expand=True)
        paste_sb.pack(side='right', fill='y')

        # Upload buttons
        upload_row = ttk.Frame(f, style='Panel.TFrame')
        upload_row.pack(fill='x', pady=(0, 6))
        ttk.Button(upload_row, text='Upload FASTA',
                   command=self._custom_upload_fasta).pack(
            side='left', padx=(0, 6))
        ttk.Button(upload_row, text='Upload .gpsat',
                   command=self._custom_upload_gpsat).pack(
            side='left', padx=(0, 6))
        ttk.Button(upload_row, text='Upload other file',
                   command=self._custom_upload_other).pack(side='left')

        # Pool list
        _section_label(f, 'Sequence pool:').pack(anchor='w', pady=(0, 2))
        lb_f = ttk.Frame(f, style='Panel.TFrame')
        lb_f.pack(fill='x', pady=(0, 4))
        self._cust_lb = _styled_listbox(lb_f, height=5)
        cust_sb = ttk.Scrollbar(lb_f, orient='vertical',
                                command=self._cust_lb.yview)
        self._cust_lb.config(yscrollcommand=cust_sb.set)
        self._cust_lb.pack(side='left', fill='both', expand=True)
        cust_sb.pack(side='right', fill='y')

        btn = ttk.Frame(f, style='Panel.TFrame')
        btn.pack(pady=(4, 6))
        ttk.Button(btn, text='Align via MAFFT & Export .gpsat',
                   command=self._custom_align_export,
                   style='Accent.TButton').pack(side='left', padx=(0, 8))
        ttk.Button(btn, text='Remove Selected',
                   command=self._custom_remove).pack(side='left', padx=(0, 6))
        ttk.Button(btn, text='Clear All',
                   command=self._custom_clear).pack(side='left')

        # Status
        self._cust_status = ttk.Label(f, text='', style='Dim.TLabel')
        self._cust_status.pack(anchor='w', pady=(0, 4))

        return f

    # ── Per-view canvas frames ────────────────────────────────────────────────

    def _build_alignment_canvas(self, parent):
        outer = tk.Frame(parent, bg=BG)

        # The grid area (Alignment or Wrapped depending on toggle)
        grid = tk.Frame(outer, bg=BG)
        grid.pack(fill='both', expand=True)

        # ── Alignment grid sub-frame ─────────────────────────────────────────
        self._align_frm = tk.Frame(grid, bg=BG)
        self._align_frm.pack(fill='both', expand=True)

        self._lbl_cv = tk.Canvas(self._align_frm, width=self.LW, bg=BG,
                                  highlightthickness=0)
        self._lbl_cv.pack(side='left', fill='y')
        right = tk.Frame(self._align_frm, bg=BG)
        right.pack(side='left', fill='both', expand=True)
        xsb = ttk.Scrollbar(right, orient='horizontal')
        self._ysb = ttk.Scrollbar(right, orient='vertical',
                                   command=self._sync_yview)
        self._seq_cv = tk.Canvas(right, bg=BG, highlightthickness=0,
                                  xscrollcommand=xsb.set,
                                  yscrollcommand=self._on_yscroll_update)
        xsb.config(command=self._seq_cv.xview)
        xsb.pack(side='bottom', fill='x')
        self._ysb.pack(side='right', fill='y')
        self._seq_cv.pack(side='left', fill='both', expand=True)
        self._tip = tk.Label(self._seq_cv, bg='#FDF0E4', fg=FG,
                              font=('Helvetica', _fs(9)), justify='left',
                              relief='solid', borderwidth=1, padx=6, pady=3)
        _register_widget_font(self._tip, 'Helvetica', 9)
        self._seq_cv.bind('<Motion>', self._on_hover)
        self._seq_cv.bind('<Leave>',  lambda _e: self._tip.place_forget())
        for cv in (self._seq_cv, self._lbl_cv):
            cv.bind('<MouseWheel>', self._on_wheel)
            cv.bind('<Button-4>',   self._on_wheel)
            cv.bind('<Button-5>',   self._on_wheel)

        # ── Wrapped sub-frame ────────────────────────────────────────────────
        self._wrap_frm = tk.Frame(grid, bg=BG)
        # not packed initially; shown by toggle

        wrap_ysb = ttk.Scrollbar(self._wrap_frm, orient='vertical')
        self._wrap_cv = tk.Canvas(self._wrap_frm, bg=BG, highlightthickness=0,
                                   yscrollcommand=wrap_ysb.set)
        wrap_ysb.config(command=self._wrap_cv.yview)
        wrap_ysb.pack(side='right', fill='y')
        self._wrap_cv.pack(side='left', fill='both', expand=True)
        self._wrap_tip = tk.Label(self._wrap_cv, bg='#FDF0E4', fg=FG,
                                   font=('Helvetica', _fs(9)), justify='left',
                                   relief='solid', borderwidth=1, padx=6, pady=3)
        _register_widget_font(self._wrap_tip, 'Helvetica', 9)
        self._wrap_cv.bind('<Configure>', self._on_wrap_resize)
        self._wrap_cv.bind('<Motion>',    self._on_wrap_hover)
        self._wrap_cv.bind('<Leave>',     lambda _e: self._wrap_tip.place_forget())
        self._wrap_cv.bind('<MouseWheel>', self._on_wrap_wheel)
        self._wrap_cv.bind('<Button-4>',   self._on_wrap_wheel)
        self._wrap_cv.bind('<Button-5>',   self._on_wrap_wheel)

        # ── LogoMaker panel (hidden initially) ───────────────────────────────
        self._logo_frm_align = tk.Frame(outer, bg=BG,
                                         relief='sunken', borderwidth=1)
        # not packed initially

        return outer

    def _build_compare_canvas(self, parent):
        outer = tk.Frame(parent, bg=BG)

        cmp_frm = tk.Frame(outer, bg=BG)
        cmp_frm.pack(fill='both', expand=True)

        cmp_ysb = ttk.Scrollbar(cmp_frm, orient='vertical')
        self._cmp_cv = tk.Canvas(cmp_frm, bg=BG, highlightthickness=0,
                                  yscrollcommand=cmp_ysb.set)
        cmp_ysb.config(command=self._cmp_cv.yview)
        cmp_ysb.pack(side='right', fill='y')
        self._cmp_cv.pack(side='left', fill='both', expand=True)
        self._cmp_tip = tk.Label(self._cmp_cv, bg='#FDF0E4', fg=FG,
                                  font=('Helvetica', _fs(9)), justify='left',
                                  relief='solid', borderwidth=1, padx=6, pady=3)
        _register_widget_font(self._cmp_tip, 'Helvetica', 9)
        self._cmp_cv.bind('<Configure>', self._on_cmp_resize)
        self._cmp_cv.bind('<Motion>',    self._on_cmp_hover)
        self._cmp_cv.bind('<Leave>',     lambda _e: self._cmp_tip.place_forget())
        self._cmp_cv.bind('<MouseWheel>', self._on_cmp_wheel)
        self._cmp_cv.bind('<Button-4>',   self._on_cmp_wheel)
        self._cmp_cv.bind('<Button-5>',   self._on_cmp_wheel)

        # LogoMaker panel for compare
        self._logo_frm_cmp = tk.Frame(outer, bg=BG,
                                       relief='sunken', borderwidth=1)

        return outer

    def _build_custom_canvas(self, parent):
        # Custom tab has no canvas — it's pure controls.
        # Return a placeholder frame.
        frm = tk.Frame(parent, bg=BG)
        ttk.Label(frm,
                  text='Use the controls above to build your sequence pool '
                       'and export a .gpsat file.',
                  style='Dim.TLabel').pack(padx=20, pady=20)
        return frm

    # ════════════════════════════════════════════════════════════════════════
    #  VIEW SWITCHING
    # ════════════════════════════════════════════════════════════════════════

    def _switch_view(self, mode=None):
        if mode is None:
            mode = self._view_mode.get()
        self._view_mode.set(mode)

        for v, btn in self._view_btns.items():
            btn.config(style='Accent.TButton' if v == mode else 'TButton')

        for v, panel in self._ctrl_panels.items():
            if v == mode:
                panel.pack(fill='x')
            else:
                panel.pack_forget()

        for v, frm in self._canvas_frames.items():
            if v == mode:
                frm.pack(fill='both', expand=True)
            else:
                frm.pack_forget()

        if mode == 'alignment' and self._length > 0:
            self._redraw_align()
        elif mode == 'compare' and self._cmp_length > 0:
            self._draw_compare()

    def _redraw_current(self):
        mode = self._view_mode.get()
        if mode == 'alignment':
            self._redraw_align()
        elif mode == 'compare':
            self._draw_compare()

    def _redraw_align(self):
        if self._align_view_mode == 'alignment':
            self._draw()
        else:
            self._draw_wrapped()

    def _toggle_align_view(self):
        if self._align_view_mode == 'alignment':
            self._align_view_mode = 'wrapped'
            self._align_view_btn.config(text='⇄ Switch to Alignment')
            self._align_frm.pack_forget()
            self._wrap_frm.pack(fill='both', expand=True)
            if self._length > 0:
                self._draw_wrapped()
        else:
            self._align_view_mode = 'alignment'
            self._align_view_btn.config(text='⇄ Switch to Wrapped')
            self._wrap_frm.pack_forget()
            self._align_frm.pack(fill='both', expand=True)
            if self._length > 0:
                self._draw()

    # ════════════════════════════════════════════════════════════════════════
    #  DATA LOADING — ALIGNMENT
    # ════════════════════════════════════════════════════════════════════════

    def _browse_align(self):
        path = filedialog.askopenfilename(
            filetypes=[('FASTA Files', '*.fa *.fasta *.txt'),
                       ('All Files', '*.*')])
        if path: self.align_path.set(path)

    def _render(self):
        align = self.align_path.get()
        if not align or not os.path.isfile(align):
            messagebox.showerror('Error', 'Select a valid alignment FASTA.')
            return
        try:
            self._load_fasta_data(align)
            self._cov_var.set(0.0)
            self._cov_lbl.config(text='Off')
            self._seqs = list(self._seqs_all)
            # If UniProt ref is loaded, re-incorporate it
            if self._uniprot_seq:
                self._incorporate_uniprot()
            else:
                self._redraw_align()
        except Exception as e:
            messagebox.showerror('Render Error', str(e))

    def _load_fasta_data(self, align_path):
        entries, _ = parse_fasta_entries(input_path=align_path)
        self._seqs_all = [(e['id'], e['sequence']) for e in entries]
        self._length   = (min(len(s) for _, s in self._seqs_all)
                          if self._seqs_all else 0)

        self._seqs_raw_cov = []
        for name, aligned_seq in self._seqs_all:
            n_total   = len(aligned_seq)
            n_residue = n_total - aligned_seq.count('-')
            cov       = n_residue / n_total if n_total > 0 else 0.0
            self._seqs_raw_cov.append((name, aligned_seq.replace('-', ''), cov))

        consensus, variability, counts_data = _votes_in_memory(self._seqs_all)
        aligned2, c2, v2, cd2 = _remove_gap_consensus_cols(
            self._seqs_all, consensus, variability, counts_data)

        self._seqs_all             = aligned2
        self._pre_uniprot_seqs_all = list(aligned2)   # stable base for UniProt re-runs
        self._seqs        = list(aligned2)
        self._consensus   = c2
        self._variability = v2
        self._counts_data = cd2
        self._counts_df   = _counts_data_to_df(cd2)
        self._length      = len(c2)

        self._orig_seqs_all    = list(self._seqs_all)
        self._orig_consensus   = self._consensus
        self._orig_variability = list(self._variability)
        self._orig_counts_data = list(self._counts_data)
        self._orig_counts_df   = self._counts_df
        self._orig_length      = self._length

        # Update variability slider max
        n = len(self._seqs_all)
        self._status.config(
            text=f'Loaded {n} sequence{"s" if n != 1 else ""} · '
                 f'{self._length} positions')

    def _import_gpsat(self):
        path = filedialog.askopenfilename(
            title='Import consensus file',
            filetypes=[('GPSAT consensus file', '*.gpsat *.vgat'),
                       ('All Files', '*.*')])
        if not path: return
        try:
            with open(path, 'r') as fh:
                data = json.load(fh)
        except Exception as e:
            messagebox.showerror('Import Error', f'Could not read file:\n{e}')
            return
        if data.get('type') not in ('gpsat_consensus', 'vgat_consensus'):
            messagebox.showerror('Import Error',
                                 'File is not a valid GPSAT consensus file.')
            return
        self._consensus   = data['consensus']
        self._length      = data['length']
        positions         = data.get('positions', [])
        self._variability = [p.get('variability', 100.0) for p in positions]
        self._counts_data = positions
        rows = []
        for p in positions:
            row = {'Position':            p['position'],
                   'Consensus':           p.get('consensus', '-'),
                   'Percent_Variability': p.get('variability', 100.0),
                   'Gap_count':           p.get('gap_count', 0),
                   'Total':               p.get('total', 0)}
            for letter, count in p.get('counts', {}).items():
                row[f'{letter}_count'] = count
            rows.append(row)
        self._counts_df       = pd.DataFrame(rows) if rows else None
        self._seqs_all        = []
        self._seqs            = []
        self._uniprot_aligned = None   # no sequences to align UniProt against
        self._trim_start      = None
        self._trim_end        = None
        self._redraw_align()
        uni_note = '  · UniProt loaded (load FASTA to display)' if self._uniprot_seq else ''
        self._status.config(
            text=f'Imported: {os.path.basename(path)}  ·  {self._length} positions{uni_note}')

    # ── UniProt reference ─────────────────────────────────────────────────────

    def _load_uniprot_ref(self):
        acc = self._uniprot_acc_var.get().strip()
        if not acc:
            messagebox.showerror('Error', 'Enter a UniProt accession.')
            return
        self._uniprot_status.config(text=f'Looking up "{acc}"…')

        def worker():
            try:
                name, seq = _fetch_uniprot_accession(acc)
                self.after(0, lambda n=name, s=seq: self._uniprot_ref_done(n, s))
            except Exception as exc:
                self.after(0, lambda m=str(exc): self._uniprot_ref_error(m))

        threading.Thread(target=worker, daemon=True).start()

    def _uniprot_ref_done(self, name, seq):
        self._uniprot_seq  = seq
        self._uniprot_name = name
        self._uniprot_status.config(
            text=f'Loaded: {name[:40]}  ({len(seq)} aa)')
        if self._seqs_all:
            self._incorporate_uniprot()

    def _uniprot_ref_error(self, msg):
        self._uniprot_status.config(text='Error — see dialog')
        messagebox.showerror('UniProt Error', msg)

    def _incorporate_uniprot(self):
        """Re-run MAFFT with UniProt sequence added, then hard-trim to UniProt length."""
        if not self._uniprot_seq or not self._seqs_all:
            return
        if self._mafft_busy:
            self._uniprot_pending = True   # retry once current MAFFT finishes
            return
        self._uniprot_pending = False
        self._mafft_busy = True
        self._status.config(text='Re-aligning with UniProt reference…')
        seqs = [(n, s.replace('-', '')) for n, s in self._pre_uniprot_seqs_all]
        # Use sentinel to avoid collisions with data sequence names
        seqs_with_uni = [(_UNI_SENTINEL, self._uniprot_seq)] + seqs

        def worker():
            try:
                aligned = _run_mafft(seqs_with_uni)
                self.after(0, lambda a=aligned: self._uniprot_mafft_done(a))
            except Exception as exc:
                self.after(0, lambda m=str(exc): self._uniprot_mafft_error(m))

        threading.Thread(target=worker, daemon=True).start()

    def _uniprot_mafft_done(self, aligned):
        self._mafft_busy = False
        # Find aligned UniProt row via sentinel (immune to name collisions)
        uni_aligned = next(
            (s for n, s in aligned if n == _UNI_SENTINEL), None)
        if uni_aligned is None:
            messagebox.showerror('Error', 'Could not find UniProt row in alignment.')
            return

        # Update sequences (all rows except the UniProt sentinel)
        other = [(n, s) for n, s in aligned if n != _UNI_SENTINEL]
        consensus, variability, counts_data = _votes_in_memory(other)

        # Remove columns where ALL data sequences AND UniProt have gaps.
        # If UniProt has a residue at a position, keep it even if consensus is '-'.
        keep = [i for i, c in enumerate(consensus)
                if c != '-' or (i < len(uni_aligned) and uni_aligned[i] != '-')]
        if keep:
            other       = [(n, ''.join(s[i] if i < len(s) else '-' for i in keep))
                           for n, s in other]
            uni_aligned = ''.join(uni_aligned[i] if i < len(uni_aligned) else '-'
                                  for i in keep)
            consensus   = ''.join(consensus[i]   for i in keep)
            variability = [variability[i]         for i in keep]
            counts_data = [{**counts_data[old_i], 'position': new_i}
                           for new_i, old_i in enumerate(keep)]
        if self._trim_to_uniprot_var.get():
            # Trim to the first AND last non-gap in UniProt aligned row
            nongap_cols = [i for i, c in enumerate(uni_aligned) if c != '-']
            if nongap_cols:
                first_nongap = nongap_cols[0]
                last_nongap  = nongap_cols[-1]
            else:
                first_nongap, last_nongap = 0, len(uni_aligned) - 1
            self._trim_start = first_nongap
            self._trim_end   = last_nongap + 1   # exclusive
        else:
            self._trim_start = None
            self._trim_end   = None

        ts = self._trim_start if self._trim_start is not None else 0
        te = self._trim_end   if self._trim_end   is not None else len(consensus)

        # Apply start + end trim to all data
        other       = [(n, s[ts:te]) for n, s in other]
        uni_aligned = uni_aligned[ts:te]
        consensus   = consensus[ts:te]
        variability = variability[ts:te]
        counts_data = [{**d, 'position': d['position'] - ts}
                       for d in counts_data if ts <= d['position'] < te]

        # Remove all columns where UniProt still has a gap (complete-match mode)
        if self._match_uniprot_var.get() and '-' in uni_aligned:
            keep = [i for i, c in enumerate(uni_aligned) if c != '-']
            if keep:
                col_remap   = {old: new for new, old in enumerate(keep)}
                other       = [(n, ''.join(s[i] if i < len(s) else '-' for i in keep))
                               for n, s in other]
                uni_aligned = ''.join(uni_aligned[i] for i in keep)
                consensus   = ''.join(consensus[i] if i < len(consensus) else '-'
                                      for i in keep)
                variability = [variability[i] for i in keep if i < len(variability)]
                counts_data = [{**d, 'position': col_remap[d['position']]}
                               for d in counts_data if d['position'] in col_remap]

        self._uniprot_aligned = uni_aligned
        self._seqs_all    = list(other)
        self._seqs        = list(other)
        self._consensus   = consensus
        self._variability = variability
        self._counts_data = counts_data
        self._counts_df   = _counts_data_to_df(self._counts_data)
        self._length      = len(consensus)

        trim_note = (f'trimmed to {len(consensus)} positions'
                     if (self._trim_start is not None or self._trim_end is not None)
                     else f'full alignment ({len(consensus)} positions)')
        self._status.config(
            text=f'UniProt ref loaded · {trim_note} '
                 f'· {len(self._seqs)} sequences')
        self._redraw_align()

    def _uniprot_mafft_error(self, msg):
        self._mafft_busy      = False
        self._uniprot_pending = False
        self._uniprot_status.config(text='Error')
        messagebox.showerror('MAFFT Error', msg)

    def _on_match_uniprot_toggle(self):
        """Re-run MAFFT+UniProt pipeline when the 'No gaps in UniProt' toggle changes."""
        if self._uniprot_seq and self._seqs_all:
            self._incorporate_uniprot()

    def _remove_uniprot_ref(self):
        """Clear the loaded UniProt reference and restore the alignment without it."""
        if not self._uniprot_seq:
            return
        self._uniprot_seq     = None
        self._uniprot_name    = None
        self._uniprot_aligned = None
        self._uniprot_acc_var.set('')
        self._uniprot_status.config(text='')
        self._trim_start = None
        self._trim_end   = None
        pct = self._cov_var.get()
        if pct >= 1.0 and self._seqs_raw_cov:
            # Coverage filter active — re-fire it; _incorporate_uniprot won't run
            # because self._uniprot_seq is now None
            self._cov_fire()
        else:
            # No coverage filter — restore original pre-UniProt alignment
            self._seqs_all    = list(self._orig_seqs_all)
            self._seqs        = list(self._orig_seqs_all)
            self._consensus   = self._orig_consensus
            self._variability = list(self._orig_variability)
            self._counts_data = list(self._orig_counts_data)
            self._counts_df   = self._orig_counts_df
            self._length      = self._orig_length
            self._redraw_align()
            if self._seqs:
                self._status.config(
                    text=f'{len(self._seqs)} sequences · {self._length} positions')

    # ── Coverage filter ───────────────────────────────────────────────────────

    def _on_cov_change(self, value=None):
        pct = self._cov_var.get()
        if pct < 1.0:
            self._cov_lbl.config(text='Off')
        else:
            n_kept = sum(1 for _, _, cov in self._seqs_raw_cov
                         if cov >= pct / 100.0)
            self._cov_lbl.config(text=f'≥{pct:.0f}%  ({n_kept})')
        if self._cov_job is not None:
            self.after_cancel(self._cov_job)
        self._cov_job = self.after(700, self._cov_fire)

    def _cov_fire(self):
        self._cov_job = None
        pct = self._cov_var.get()

        if pct < 1.0:
            self._seqs_all             = list(self._orig_seqs_all)
            self._pre_uniprot_seqs_all = list(self._orig_seqs_all)
            self._seqs        = list(self._orig_seqs_all)
            self._consensus   = self._orig_consensus
            self._variability = list(self._orig_variability)
            self._counts_data = list(self._orig_counts_data)
            self._counts_df   = self._orig_counts_df
            self._length      = self._orig_length
            self._trim_start  = None
            self._trim_end    = None
            self._uniprot_aligned = None
            if self._uniprot_seq:
                self._incorporate_uniprot()
            else:
                self._redraw_align()
                self._status.config(
                    text=f'{len(self._seqs)} sequences · {self._length} positions  '
                         f'(coverage filter off)')
            return

        if not self._seqs_raw_cov:
            return
        if self._mafft_busy:
            return

        min_cov = pct / 100.0
        filtered = [(name, seq)
                    for name, seq, cov in self._seqs_raw_cov
                    if cov >= min_cov]

        if len(filtered) < 2:
            self._status.config(
                text=f'Coverage ≥{pct:.0f}%: only {len(filtered)} sequence(s) — '
                     f'threshold too high.')
            return

        self._mafft_busy = True
        self._status.config(
            text=f'Realigning {len(filtered)} sequences via MAFFT …')
        threading.Thread(target=self._cov_mafft_worker,
                         args=(filtered,), daemon=True).start()

    def _cov_mafft_worker(self, seqs):
        try:
            aligned = _run_mafft(seqs)
            consensus, variability, counts_data = _votes_in_memory(aligned)
            a2, c2, v2, cd2 = _remove_gap_consensus_cols(
                aligned, consensus, variability, counts_data)
            self.after(0, self._cov_mafft_done, a2, c2, v2, cd2)
        except Exception as exc:
            self.after(0, self._cov_mafft_error, str(exc))

    def _cov_mafft_done(self, aligned, consensus, variability, counts_data):
        self._mafft_busy           = False
        self._seqs_all             = aligned
        self._pre_uniprot_seqs_all = list(aligned)   # stable base for UniProt re-runs
        self._seqs        = list(aligned)
        self._consensus   = consensus
        self._variability = variability
        self._counts_data = counts_data
        self._counts_df   = _counts_data_to_df(counts_data)
        self._length      = len(consensus)
        self._trim_start  = None
        self._trim_end    = None
        pct = self._cov_var.get()
        self._status.config(
            text=f'{len(aligned)} sequences · {len(consensus)} positions  '
                 f'(coverage ≥ {pct:.0f}%)')
        if self._uniprot_seq:
            self._uniprot_aligned = None
            self._incorporate_uniprot()
        else:
            self._redraw_align()

    # ── Variability-N slider ──────────────────────────────────────────────────

    def _on_var_n_change(self, value=None):
        n = self._var_n_var.get()
        total = len(self._seqs_all)
        if n == 0 or n >= total:
            self._var_n_lbl.config(text='All seqs')
        else:
            self._var_n_lbl.config(text=f'First {n}')
        if hasattr(self, '_var_n_job'):
            self.after_cancel(self._var_n_job)
        self._var_n_job = self.after(700, self._var_n_fire)

    def _var_n_fire(self):
        if self._length == 0:
            return
        n = self._var_n_var.get()
        total = len(self._seqs_all)
        use_seqs = self._seqs_all[:n] if (0 < n < total) else self._seqs_all
        consensus, variability, counts_data = _votes_in_memory(use_seqs)
        aligned2, c2, v2, cd2 = _remove_gap_consensus_cols(
            use_seqs, consensus, variability, counts_data)
        self._consensus   = c2
        self._variability = v2
        self._counts_data = cd2
        self._counts_df   = _counts_data_to_df(cd2)
        self._length      = len(c2)
        self._redraw_align()

    # ════════════════════════════════════════════════════════════════════════
    #  COLOUR HELPERS
    # ════════════════════════════════════════════════════════════════════════

    @staticmethod
    def _var_to_color(variability):
        c = 1.0 - min(max(variability, 0.0), 100.0) / 100.0
        r = int(255 - 64  * c)
        g = int(255 - 168 * c)
        b = int(255 - 255 * c)
        return f'#{r:02x}{g:02x}{b:02x}'

    # ════════════════════════════════════════════════════════════════════════
    #  DRAWING — ALIGNMENT VIEW
    # ════════════════════════════════════════════════════════════════════════

    def _draw(self):
        if self._length == 0:
            return
        # Build rows: CONSENSUS first, then UNIPROT if loaded, then sequences
        rows = self._build_display_rows()
        self._draw_alignment_grid(
            self._seq_cv, self._lbl_cv, rows,
            self._consensus, self._variability, self._length)
        n = len(self._seqs)
        self._status.config(
            text=f'{n} sequence{"s" if n != 1 else ""} · {self._length} positions  '
                 '| orange=conserved  white=variable  peach=mismatch')

    def _build_display_rows(self):
        """Build the list of (label, seq, is_special) rows for the grid.
        is_special: 'consensus' | 'uniprot' | None"""
        rows = [('CONSENSUS', self._consensus, 'consensus')]
        if self._uniprot_aligned is not None:
            rows.append((self._uniprot_name or 'UNIPROT', self._uniprot_aligned, 'uniprot'))
        for name, seq in self._seqs:
            rows.append((name, seq, None))
        return rows

    def _draw_alignment_grid(self, seq_cv, lbl_cv, rows,
                              consensus, variability, length):
        CW, CH, LW, RH = self.CW, self.CH, self.LW, self.RULER_H
        n_pos  = length
        n_rows = len(rows)

        lbl_cv.config(width=LW)
        seq_cv.delete('all')
        lbl_cv.delete('all')

        total_w = n_pos * CW
        total_h = RH + n_rows * CH
        seq_cv.config(scrollregion=(0, 0, total_w, total_h), yscrollincrement=CH)
        lbl_cv.config(scrollregion=(0, 0, LW,      total_h), yscrollincrement=CH)

        # Ruler
        lbl_cv.create_rectangle(0, 0, LW, RH, fill=PANEL, outline='')
        lbl_cv.create_text(6, RH // 2, text='Position', anchor='w',
                            fill=DIM_FG, font=('Helvetica', _fs(8), 'italic'))
        seq_cv.create_rectangle(0, 0, total_w, RH, fill=PANEL, outline='')
        for col in range(n_pos):
            pos1 = col + 1
            if pos1 == 1 or pos1 % 10 == 0:
                xc = col * CW + CW // 2
                seq_cv.create_text(xc, RH // 2, text=str(pos1),
                                    anchor='center', fill=DIM_FG,
                                    font=('Helvetica', _fs(7)))
                seq_cv.create_line(xc, RH - 3, xc, RH,
                                    fill=BORDER, width=1)
        seq_cv.create_line(0, RH, total_w, RH, fill=BORDER, width=1)
        lbl_cv.create_line(0, RH, LW,     RH, fill=BORDER, width=1)

        col_colors = [self._var_to_color(variability[c])
                      for c in range(n_pos)]

        for row_i, (name, seq, special) in enumerate(rows):
            y0   = RH + row_i * CH
            yc   = y0 + CH // 2

            if special == 'consensus':
                lbl_bg, lbl_fg = PANEL, ACCENT
                font_style = 'bold'
            elif special == 'uniprot':
                lbl_bg, lbl_fg = UNI_BG, UNI_FG
                font_style = 'bold'
            else:
                lbl_bg, lbl_fg = BG, DIM_FG
                font_style = 'normal'

            lbl_cv.create_rectangle(0, y0, LW, y0 + CH, fill=lbl_bg, outline='')
            disp = (name[:20] + '…') if len(name) > 21 else name
            lbl_cv.create_text(6, yc, text=disp, anchor='w', fill=lbl_fg,
                                font=('Courier', _fs(8), font_style))

            for col in range(min(n_pos, len(seq))):
                x0  = col * CW
                ch  = seq[col]
                var = variability[col] if col < len(variability) else 100.0
                bg  = col_colors[col]
                con = consensus[col] if col < len(consensus) else ''
                conservation = 1.0 - var / 100.0
                txt_fg = '#ffffff' if conservation > 0.60 else '#1A1A1A'

                if special == 'consensus':
                    cell_bg, cell_fg = bg, txt_fg
                elif special == 'uniprot':
                    if ch in ('-', 'N', 'X'):
                        cell_bg, cell_fg = UNI_BG, BORDER
                    else:
                        cell_bg, cell_fg = UNI_BG, UNI_FG
                elif ch in ('-', 'N', 'X'):
                    cell_bg, cell_fg = BG, BORDER
                elif ch == con:
                    cell_bg, cell_fg = bg, txt_fg
                else:
                    cell_bg, cell_fg = '#FFE8D8', '#A33000'

                seq_cv.create_rectangle(x0, y0, x0 + CW, y0 + CH,
                                         fill=cell_bg, outline='')
                seq_cv.create_text(x0 + CW // 2, yc, text=ch,
                                    anchor='center', fill=cell_fg,
                                    font=('Courier', _fs(7)))

        # Divider after consensus row
        seq_cv.create_line(0, RH + CH, total_w, RH + CH,
                            fill=ACCENT, width=1)
        lbl_cv.create_line(0, RH + CH, LW,      RH + CH,
                            fill=ACCENT, width=1)

    # ── Alignment hover ───────────────────────────────────────────────────────

    def _on_hover(self, event):
        self._grid_hover(event, self._seq_cv, self._tip,
                         self._consensus, self._variability,
                         self._length, self._counts_df)

    def _grid_hover(self, event, cv, tip,
                    consensus, variability, length, counts_df):
        CW, RH = self.CW, self.RULER_H
        cx = cv.canvasx(event.x)
        cy = cv.canvasy(event.y)
        if cy < RH:
            tip.place_forget(); return
        col = int(cx / CW)
        if col < 0 or col >= length:
            tip.place_forget(); return
        var = variability[col] if col < len(variability) else 0.0
        con = consensus[col]   if col < len(consensus)   else '?'
        lines = [f'Position {col + 1}',
                 f'Consensus : {con}',
                 f'Conservation : {100 - var:.1f}%',
                 f'Variability  : {var:.1f}%']
        if counts_df is not None and col < len(counts_df):
            r = counts_df.iloc[col]
            cnt_cols = [c for c in counts_df.columns
                        if c.endswith('_count')
                        and c not in ('Gap_count', 'Other_count')]
            top = sorted([(c.replace('_count', ''), int(r[c]))
                           for c in cnt_cols if r[c] > 0],
                          key=lambda x: -x[1])[:6]
            if top:
                lines.append('Counts : ' +
                              '  '.join(f'{c}:{n}' for c, n in top))
            lines.append(f'Gaps : {int(r.get("Gap_count", 0))} / '
                         f'{int(r.get("Total", 0))}')
        tip.config(text='\n'.join(lines))
        tip.place(x=event.x + 14, y=event.y + 14)

    # ── Scroll helpers ────────────────────────────────────────────────────────

    def _sync_yview(self, *args):
        self._seq_cv.yview(*args)
        self._lbl_cv.yview(*args)

    def _on_yscroll_update(self, first, last):
        self._ysb.set(first, last)
        self._lbl_cv.yview_moveto(first)

    def _on_wheel(self, event):
        if event.num == 4:
            units = -1
        elif event.num == 5:
            units = 1
        elif event.delta:
            # macOS: delta is small (1–20 px); Windows: multiples of 120
            units = -1 if event.delta > 0 else 1
        else:
            return
        self._seq_cv.yview_scroll(units, 'units')
        self._lbl_cv.yview_scroll(units, 'units')

    # ════════════════════════════════════════════════════════════════════════
    #  DRAWING — WRAPPED VIEW
    # ════════════════════════════════════════════════════════════════════════

    def _draw_wrapped(self):
        if self._length == 0:
            return
        CW, CH, LW, RH = self.WRAP_CW, self.WRAP_CH, self.LW, self.RULER_H
        BLOCK_GAP = 10
        canvas_w = self._wrap_cv.winfo_width()
        if canvas_w <= 1: canvas_w = 700
        cols_per_block = max(1, (canvas_w - LW) // CW)
        n_pos  = self._length
        all_rows = self._build_display_rows()
        n_rpb  = len(all_rows)
        block_h = n_rpb * CH
        n_blocks = (n_pos + cols_per_block - 1) // cols_per_block
        total_h  = n_blocks * (RH + block_h + BLOCK_GAP)
        self._wrap_cv.delete('all')
        self._wrap_cv.config(scrollregion=(0, 0, canvas_w, total_h), yscrollincrement=CH)
        self._wrap_layout = {
            'cols_per_block': cols_per_block, 'n_all_rows': n_rpb,
            'block_h': block_h, 'block_gap': BLOCK_GAP,
            'ruler_h': RH, 'CW': CW, 'CH': CH, 'LW': LW, 'n_pos': n_pos}
        try:
            start_aa = max(1, int(self._start_aa_var.get()))
        except (ValueError, AttributeError):
            start_aa = 1
        col_colors = [self._var_to_color(v) for v in self._variability[:n_pos]]
        for block_i in range(n_blocks):
            start_col = block_i * cols_per_block
            end_col   = min(start_col + cols_per_block, n_pos)
            block_y   = block_i * (RH + block_h + BLOCK_GAP)
            self._wrap_cv.create_rectangle(0, block_y, canvas_w,
                                            block_y + RH, fill=PANEL, outline='')
            self._wrap_cv.create_text(
                6, block_y + RH // 2, text='Pos',
                anchor='w', fill=DIM_FG,
                font=('Helvetica', _fs(8), 'italic'))
            for col in range(start_col, end_col):
                pos1 = col + start_aa
                if col == start_col or pos1 % 10 == 0:
                    xc = LW + (col - start_col) * CW + CW // 2
                    self._wrap_cv.create_text(
                        xc, block_y + RH // 2, text=str(pos1),
                        anchor='center', fill=DIM_FG,
                        font=('Helvetica', _fs(7)))
                    self._wrap_cv.create_line(
                        xc, block_y + RH - 3, xc, block_y + RH,
                        fill=BORDER, width=1)
            self._wrap_cv.create_line(0, block_y + RH, canvas_w,
                                       block_y + RH, fill=BORDER, width=1)
            for row_i, (name, seq, special) in enumerate(all_rows):
                y0   = block_y + RH + row_i * CH
                yc   = y0 + CH // 2
                if special == 'consensus':
                    lbl_bg, lbl_fg, font_style = PANEL, ACCENT, 'bold'
                elif special == 'uniprot':
                    lbl_bg, lbl_fg, font_style = UNI_BG, UNI_FG, 'bold'
                else:
                    lbl_bg, lbl_fg, font_style = BG, DIM_FG, 'normal'
                self._wrap_cv.create_rectangle(0, y0, LW, y0 + CH,
                                                fill=lbl_bg, outline='')
                disp = (name[:20] + '…') if len(name) > 21 else name
                self._wrap_cv.create_text(
                    6, yc, text=disp, anchor='w', fill=lbl_fg,
                    font=('Courier', _fs(9), font_style))
                for col in range(start_col, end_col):
                    if col >= len(seq): continue
                    x0  = LW + (col - start_col) * CW
                    ch  = seq[col]
                    var = self._variability[col] if col < len(self._variability) else 100.0
                    bg  = col_colors[col]
                    con = self._consensus[col] if col < len(self._consensus) else ''
                    conservation = 1.0 - var / 100.0
                    txt_fg = '#ffffff' if conservation > 0.60 else '#1A1A1A'
                    if special == 'consensus':
                        cell_bg, cell_fg = bg, txt_fg
                    elif special == 'uniprot':
                        cell_bg, cell_fg = (UNI_BG, BORDER) if ch in ('-', 'N', 'X') else (UNI_BG, UNI_FG)
                    elif ch in ('-', 'N', 'X'):
                        cell_bg, cell_fg = BG, BORDER
                    elif ch == con:
                        cell_bg, cell_fg = bg, txt_fg
                    else:
                        cell_bg, cell_fg = '#FFE8D8', '#A33000'
                    self._wrap_cv.create_rectangle(x0, y0, x0 + CW, y0 + CH,
                                                    fill=cell_bg, outline='')
                    self._wrap_cv.create_text(x0 + CW // 2, yc, text=ch,
                                               anchor='center', fill=cell_fg,
                                               font=('Courier', _fs(10)))
            self._wrap_cv.create_line(
                0, block_y + RH + CH, canvas_w, block_y + RH + CH,
                fill=ACCENT, width=1)
        n = len(self._seqs)
        cov = self._cov_var.get()
        cov_str = f'  ·  cov ≥ {cov:.0f}%' if cov >= 1 else ''
        self._status.config(
            text=f'{n} seq · {n_pos} pos · {n_blocks} block(s){cov_str}')

    def _on_wrap_resize(self, event):
        if hasattr(self, '_wrap_resize_job'):
            self.after_cancel(self._wrap_resize_job)
        self._wrap_resize_job = self.after(
            80, lambda: self._draw_wrapped()
            if self._length > 0 and self._align_view_mode == 'wrapped'
            else None)

    def _on_wrap_hover(self, event):
        layout = self._wrap_layout
        if not layout: return
        cpb   = layout['cols_per_block']
        bh    = layout['block_h']; bg_gap = layout['block_gap']
        RH    = layout['ruler_h']
        CW    = layout['CW'];  LW = layout['LW']
        n_pos = layout['n_pos']
        cx = self._wrap_cv.canvasx(event.x)
        cy = self._wrap_cv.canvasy(event.y)
        bt  = RH + bh + bg_gap
        bi  = int(cy / bt)
        yib = cy - bi * bt
        if yib < RH or yib >= RH + bh or cx < LW:
            self._wrap_tip.place_forget(); return
        cib = int((cx - LW) / CW)
        col = bi * cpb + cib
        if col >= n_pos or cib >= cpb:
            self._wrap_tip.place_forget(); return
        var = self._variability[col] if col < len(self._variability) else 0.0
        con = self._consensus[col]   if col < len(self._consensus)   else '?'
        lines = [f'Position {col + 1}', f'Consensus : {con}',
                 f'Conservation : {100 - var:.1f}%',
                 f'Variability  : {var:.1f}%']
        if self._counts_df is not None and col < len(self._counts_df):
            r = self._counts_df.iloc[col]
            cnt_cols = [c for c in self._counts_df.columns
                        if c.endswith('_count')
                        and c not in ('Gap_count', 'Other_count')]
            top = sorted([(c.replace('_count', ''), int(r[c]))
                           for c in cnt_cols if r[c] > 0],
                          key=lambda x: -x[1])[:6]
            if top:
                lines.append('Counts : ' +
                              '  '.join(f'{c}:{n}' for c, n in top))
        self._wrap_tip.config(text='\n'.join(lines))
        self._wrap_tip.place(x=event.x + 14, y=event.y + 14)

    def _on_wrap_wheel(self, event):
        if event.num == 4:
            units = -1
        elif event.num == 5:
            units = 1
        elif event.delta:
            units = -1 if event.delta > 0 else 1
        else:
            return
        self._wrap_cv.yview_scroll(units, 'units')

    # ════════════════════════════════════════════════════════════════════════
    #  MOST REPRESENTATIVE SEQUENCE — panel dialog
    # ════════════════════════════════════════════════════════════════════════

    def _find_most_representative(self):
        if not self._seqs or not self._consensus:
            messagebox.showerror('No data', 'Render an alignment first.')
            return
        name, score = _most_representative_sequence(self._seqs, self._consensus)
        if name is None:
            messagebox.showinfo('Result', 'No sequences to compare.')
            return

        # Get the raw (ungapped) sequence for export
        raw_seq = next((s.replace('-', '') for n, s in self._seqs if n == name), '')
        pct     = score * 100
        non_gap = sum(1 for c in self._consensus if c not in ('-', 'N', 'X'))
        matches = round(score * non_gap)

        self._status.config(
            text=f'Most Representative: {name[:50]}  ({pct:.1f}% match)')

        # Open dialog window
        dlg = tk.Toplevel(self)
        dlg.title('Most Representative Sequence')
        dlg.geometry('620x420')
        dlg.configure(bg=PANEL)
        dlg.resizable(True, True)

        ttk.Label(dlg,
                  text=f'Most representative sequence: {name}',
                  style='Section.TLabel').pack(anchor='w', padx=14, pady=(12, 2))
        ttk.Label(dlg,
                  text=f'Similarity to consensus: {pct:.1f}%  '
                       f'({matches} / {non_gap} positions)',
                  style='Dim.TLabel').pack(anchor='w', padx=14, pady=(0, 8))

        seq_frame = ttk.Frame(dlg, style='Panel.TFrame')
        seq_frame.pack(fill='both', expand=True, padx=14, pady=(0, 8))
        seq_text = _styled_text(seq_frame, wrap='char', state='normal')
        seq_sb   = ttk.Scrollbar(seq_frame, orient='vertical',
                                  command=seq_text.yview)
        seq_text.config(yscrollcommand=seq_sb.set)
        seq_text.pack(side='left', fill='both', expand=True)
        seq_sb.pack(side='right', fill='y')
        seq_text.insert(tk.END, raw_seq)
        seq_text.config(state='disabled')

        btn_row = ttk.Frame(dlg, style='Panel.TFrame')
        btn_row.pack(pady=(0, 12))

        def _copy():
            dlg.clipboard_clear()
            dlg.clipboard_append(raw_seq)
            dlg.update()

        def _export():
            path = filedialog.asksaveasfilename(
                title='Export sequence as .gpsat',
                defaultextension='.gpsat',
                filetypes=[('GPSAT file', '*.gpsat'), ('All Files', '*.*')])
            if not path: return
            seq_name = os.path.splitext(os.path.basename(path))[0]
            try:
                _write_gpsat_single_seq(path, seq_name, name, raw_seq)
                messagebox.showinfo('Saved', f'Sequence saved:\n{path}', parent=dlg)
            except Exception as e:
                messagebox.showerror('Error', str(e), parent=dlg)

        ttk.Button(btn_row, text='Copy Sequence',
                   command=_copy).pack(side='left', padx=(0, 8))
        ttk.Button(btn_row, text='Export .gpsat',
                   command=_export).pack(side='left', padx=(0, 8))
        ttk.Button(btn_row, text='Close',
                   command=dlg.destroy).pack(side='left')

    # ════════════════════════════════════════════════════════════════════════
    #  EXPORT GPSAT — dialog to choose which sequence
    # ════════════════════════════════════════════════════════════════════════

    def _export_gpsat_dialog(self):
        if not self._consensus and not self._seqs:
            messagebox.showerror('No data', 'Render an alignment first.')
            return

        dlg = tk.Toplevel(self)
        dlg.title('Export .gpsat')
        dlg.geometry('480x320')
        dlg.configure(bg=PANEL)
        dlg.resizable(True, True)

        ttk.Label(dlg, text='Select sequence to export:',
                  style='Section.TLabel').pack(anchor='w', padx=14, pady=(12, 4))

        lb_frame = ttk.Frame(dlg, style='Panel.TFrame')
        lb_frame.pack(fill='both', expand=True, padx=14, pady=(0, 8))
        lb = _styled_listbox(lb_frame, selectmode='single')
        lb_sb = ttk.Scrollbar(lb_frame, orient='vertical', command=lb.yview)
        lb.config(yscrollcommand=lb_sb.set)
        lb.pack(side='left', fill='both', expand=True)
        lb_sb.pack(side='right', fill='y')

        # Populate: consensus first, then individual sequences
        lb.insert(tk.END, 'CONSENSUS  (computed consensus sequence)')
        for name, seq in self._seqs:
            preview = seq.replace('-', '')[:20]
            lb.insert(tk.END, f'{name}  [{preview}…]')

        def _do_export():
            sel = lb.curselection()
            if not sel:
                messagebox.showinfo('Select', 'Choose a sequence first.',
                                    parent=dlg)
                return
            idx = sel[0]
            path = filedialog.asksaveasfilename(
                title='Save .gpsat file', defaultextension='.gpsat',
                filetypes=[('GPSAT file', '*.gpsat'), ('All Files', '*.*')],
                parent=dlg)
            if not path: return
            file_name = os.path.splitext(os.path.basename(path))[0]
            try:
                if idx == 0:
                    # Export consensus
                    _write_gpsat(path, file_name, self._consensus,
                                self._variability, self._counts_data)
                else:
                    # Export individual sequence.
                    # self._seqs always reflects the current alignment state —
                    # i.e. it is already trimmed/column-filtered if "Trim to UniProt
                    # length" or "No gaps in UniProt" are active.  Removing gaps
                    # gives the sequence exactly as used in the displayed alignment.
                    seq_name, raw_aligned = self._seqs[idx - 1]
                    raw = raw_aligned.replace('-', '')
                    _write_gpsat_single_seq(path, file_name, seq_name, raw)
                messagebox.showinfo('Saved', f'Saved:\n{path}', parent=dlg)
                dlg.destroy()
            except Exception as e:
                messagebox.showerror('Error', str(e), parent=dlg)

        btn_row = ttk.Frame(dlg, style='Panel.TFrame')
        btn_row.pack(pady=(0, 12))
        ttk.Button(btn_row, text='Export Selected',
                   command=_do_export,
                   style='Accent.TButton').pack(side='left', padx=(0, 8))
        ttk.Button(btn_row, text='Cancel',
                   command=dlg.destroy).pack(side='left')

    # ── Import .gpsat already defined above

    # ════════════════════════════════════════════════════════════════════════
    #  LOGOMAKER TOGGLE — ALIGNMENT
    # ════════════════════════════════════════════════════════════════════════

    def _toggle_logo_align(self):
        if not _LOGO_AVAILABLE:
            detail = f'\n\nError: {_LOGO_IMPORT_ERROR}' if _LOGO_IMPORT_ERROR else ''
            messagebox.showerror(
                'logomaker not available',
                'Run the app from your conda terminal:\n'
                '  conda activate <env> && python main.py\n\n'
                'Or install into the active Python:\n'
                '  pip install logomaker matplotlib' + detail)
            return
        if not self._seqs and not self._consensus:
            messagebox.showerror('No data', 'Render an alignment first.')
            return
        if self._logo_visible_align:
            self._hide_logo(self._logo_frm_align, '_logo_canvas_align',
                            '_logo_fig_align')
            self._logo_visible_align = False
            self._logo_btn_align.config(text='Show Sequence Logo ▼')
            self._logo_ctrl_frm_align.pack_forget()
        else:
            start, end = self._parse_logo_range()
            self._show_logo(self._logo_frm_align, '_logo_canvas_align',
                            '_logo_fig_align',
                            self._seqs, self._consensus,
                            start_pos=start, end_pos=end)
            self._logo_visible_align = True
            self._logo_btn_align.config(text='Hide Sequence Logo ▲')
            self._logo_ctrl_frm_align.pack(fill='x', pady=(4, 0))

    def _toggle_logo_cmp(self):
        if not _LOGO_AVAILABLE:
            detail = f'\n\nError: {_LOGO_IMPORT_ERROR}' if _LOGO_IMPORT_ERROR else ''
            messagebox.showerror(
                'logomaker not available',
                'Run the app from your conda terminal:\n'
                '  conda activate <env> && python main.py\n\n'
                'Or install into the active Python:\n'
                '  pip install logomaker matplotlib' + detail)
            return
        if not self._cmp_seqs:
            messagebox.showerror('No data', 'Run a comparison first.')
            return
        if self._logo_visible_cmp:
            self._hide_logo(self._logo_frm_cmp, '_logo_canvas_cmp',
                            '_logo_fig_cmp')
            self._logo_visible_cmp = False
            self._logo_btn_cmp.config(text='Show Sequence Logo ▼')
            self._logo_ctrl_frm_cmp.pack_forget()
        else:
            start, end = self._parse_logo_cmp_range()
            student_rows = [(n, s) for n, s in self._cmp_seqs
                            if n != self._cmp_master_name]
            all_seqs = ([(self._cmp_master_name or 'MASTER',
                          self._cmp_consensus)] + student_rows)
            self._show_logo(self._logo_frm_cmp, '_logo_canvas_cmp',
                            '_logo_fig_cmp',
                            all_seqs, self._cmp_consensus,
                            start_pos=start, end_pos=end)
            self._logo_visible_cmp = True
            self._logo_btn_cmp.config(text='Hide Sequence Logo ▲')
            self._logo_ctrl_frm_cmp.pack(fill='x', pady=(4, 0))

    def _show_logo(self, frame, canvas_attr, fig_attr, seqs, consensus,
                   start_pos=None, end_pos=None):
        """Render a sequence logo using logomaker and embed it in frame."""
        try:
            fig = _render_logo(seqs, consensus, start_pos=start_pos, end_pos=end_pos)
        except Exception as e:
            messagebox.showerror('Logo Error', str(e))
            return
        # Destroy previous canvas if any
        prev = getattr(self, canvas_attr, None)
        if prev is not None:
            try:
                prev.get_tk_widget().destroy()
            except Exception:
                pass
        prev_fig = getattr(self, fig_attr, None)
        if prev_fig is not None:
            try:
                plt.close(prev_fig)
            except Exception:
                pass
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)
        setattr(self, canvas_attr, canvas)
        setattr(self, fig_attr, fig)
        frame.pack(fill='x', padx=0, pady=(4, 0))

    def _hide_logo(self, frame, canvas_attr, fig_attr):
        prev = getattr(self, canvas_attr, None)
        if prev is not None:
            try:
                prev.get_tk_widget().destroy()
            except Exception:
                pass
        prev_fig = getattr(self, fig_attr, None)
        if prev_fig is not None:
            try:
                plt.close(prev_fig)
            except Exception:
                pass
        setattr(self, canvas_attr, None)
        setattr(self, fig_attr, None)
        frame.pack_forget()

    def _parse_logo_range(self):
        """Return (start_pos, end_pos) as 1-indexed ints from the entry fields, or (None, None)."""
        try:
            s = int(self._logo_start_var.get().strip()) if self._logo_start_var.get().strip() else None
        except ValueError:
            s = None
        try:
            e = int(self._logo_end_var.get().strip()) if self._logo_end_var.get().strip() else None
        except ValueError:
            e = None
        return s, e

    def _refresh_logo_align(self):
        """Re-render the alignment logo with the current position range."""
        if not self._logo_visible_align:
            return
        start, end = self._parse_logo_range()
        self._show_logo(self._logo_frm_align, '_logo_canvas_align',
                        '_logo_fig_align',
                        self._seqs, self._consensus,
                        start_pos=start, end_pos=end)

    def _export_logo_align(self):
        """Save the current alignment logo as a PNG file."""
        fig = getattr(self, '_logo_fig_align', None)
        if fig is None:
            messagebox.showerror('No logo', 'Show the sequence logo first.')
            return
        path = filedialog.asksaveasfilename(
            title='Export Logo PNG', defaultextension='.png',
            filetypes=[('PNG image', '*.png'), ('All Files', '*.*')])
        if not path:
            return
        try:
            fig.savefig(path, dpi=150, bbox_inches='tight')
            messagebox.showinfo('Saved', f'Logo saved:\n{path}')
        except Exception as e:
            messagebox.showerror('Export Error', str(e))

    # ── Compare logo helpers ──────────────────────────────────────────────────

    def _parse_logo_cmp_range(self):
        """Return (start_pos, end_pos) as 1-indexed ints from the compare logo
        entry fields, or (None, None)."""
        try:
            s = int(self._logo_cmp_start_var.get().strip()) if self._logo_cmp_start_var.get().strip() else None
        except ValueError:
            s = None
        try:
            e = int(self._logo_cmp_end_var.get().strip()) if self._logo_cmp_end_var.get().strip() else None
        except ValueError:
            e = None
        return s, e

    def _refresh_logo_cmp(self):
        """Re-render the compare logo with the current position range."""
        if not self._logo_visible_cmp:
            return
        start, end = self._parse_logo_cmp_range()
        student_rows = [(n, s) for n, s in self._cmp_seqs
                        if n != self._cmp_master_name]
        all_seqs = [(self._cmp_master_name or 'MASTER',
                     self._cmp_consensus)] + student_rows
        self._show_logo(self._logo_frm_cmp, '_logo_canvas_cmp', '_logo_fig_cmp',
                        all_seqs, self._cmp_consensus,
                        start_pos=start, end_pos=end)

    def _export_logo_cmp(self):
        """Save the current compare logo as a PNG file."""
        fig = getattr(self, '_logo_fig_cmp', None)
        if fig is None:
            messagebox.showerror('No logo', 'Show the sequence logo first.')
            return
        path = filedialog.asksaveasfilename(
            title='Export Logo PNG', defaultextension='.png',
            filetypes=[('PNG image', '*.png'), ('All Files', '*.*')])
        if not path:
            return
        try:
            fig.savefig(path, dpi=150, bbox_inches='tight')
            messagebox.showinfo('Saved', f'Logo saved:\n{path}')
        except Exception as e:
            messagebox.showerror('Export Error', str(e))

    # ════════════════════════════════════════════════════════════════════════
    #  EXPORT PNG
    # ════════════════════════════════════════════════════════════════════════

    def _ask_start_aa(self, default=1):
        """Modal dialog for start position number."""
        result = [None]
        dlg = tk.Toplevel(self)
        dlg.title('Start Position')
        dlg.transient(self)          # attach to parent window
        dlg.resizable(False, False)
        dlg.configure(bg=PANEL)

        ttk.Label(dlg, text='Starting residue number:',
                  style='Section.TLabel').pack(padx=20, pady=(14, 4))

        # No StringVar — read directly from widget to avoid update-lag bugs
        entry = ttk.Entry(dlg, width=12, justify='center')
        entry.insert(0, str(default))
        entry.pack(padx=20, pady=(0, 10))

        def _ok(event=None):
            try:
                v = int(entry.get().strip())
                if v >= 1:
                    result[0] = v
                    dlg.destroy()
            except ValueError:
                pass

        def _cancel(event=None):
            dlg.destroy()

        entry.bind('<Return>',   _ok)
        entry.bind('<KP_Enter>', _ok)     # numpad Enter
        entry.bind('<Escape>',   _cancel)

        btn_row = ttk.Frame(dlg, style='Panel.TFrame')
        btn_row.pack(pady=(0, 12))
        ttk.Button(btn_row, text='OK',     command=_ok,
                   style='Accent.TButton').pack(side='left', padx=(0, 6))
        ttk.Button(btn_row, text='Cancel', command=_cancel).pack(side='left')

        # Lay out, center, then grab + force focus so keystrokes land in entry
        dlg.update_idletasks()
        px = self.winfo_rootx() + (self.winfo_width()  - dlg.winfo_width())  // 2
        py = self.winfo_rooty() + (self.winfo_height() - dlg.winfo_height()) // 2
        dlg.geometry(f'+{px}+{py}')
        dlg.grab_set()               # grab AFTER geometry so window is visible
        dlg.lift()                   # bring above main window
        entry.focus_force()          # force keyboard focus into the entry
        entry.select_range(0, 'end') # pre-select existing text

        self.wait_window(dlg)
        return result[0]

    def _export_png_align(self):
        if not self._seqs and not self._consensus:
            messagebox.showerror('No data', 'Render an alignment first.')
            return
        try:
            start_aa = max(1, int(self._start_aa_var.get()))
        except (ValueError, AttributeError):
            start_aa = 1
        cov = int(self._cov_var.get())
        view = self._align_view_mode
        default_name = f'alignment_{view}_cov{cov}pct.svg'
        path = filedialog.asksaveasfilename(
            title='Export as SVG', defaultextension='.svg',
            initialfile=default_name,
            filetypes=[('SVG image', '*.svg'), ('All Files', '*.*')])
        if not path: return
        rows = self._build_display_rows()
        seqs_for_export = [(name, seq) for name, seq, _ in rows[1:]]
        try:
            svg = self._build_svg_wrapped(
                seqs_for_export, self._consensus, self._variability,
                pos_labels=None,
                uniprot_name=self._uniprot_name if self._uniprot_aligned else None,
                start_aa=start_aa)
            with open(path, 'w', encoding='utf-8') as fh:
                fh.write(svg)
            messagebox.showinfo('Saved', f'SVG saved:\n{path}')
        except Exception as e:
            messagebox.showerror('Export Error', str(e))

    def _export_png_compare(self):
        if not self._cmp_seqs:
            messagebox.showerror('No data', 'Run a comparison first.')
            return
        try:
            start_aa = max(1, int(self._cmp_start_aa_var.get()))
        except (ValueError, AttributeError):
            start_aa = 1
        path = filedialog.asksaveasfilename(
            title='Export comparison as SVG', defaultextension='.svg',
            filetypes=[('SVG image', '*.svg'), ('All Files', '*.*')])
        if not path: return
        master_seq = self._cmp_consensus
        student_rows = [(n, s) for n, s in self._cmp_seqs
                        if n != self._cmp_master_name]
        seqs = student_rows
        try:
            svg = self._build_svg_wrapped(
                seqs, master_seq, self._cmp_variability,
                pos_labels=self._cmp_master_nums or None,
                master_name=self._cmp_master_name,
                start_aa=start_aa,
                master_seq=master_seq)
            with open(path, 'w', encoding='utf-8') as fh:
                fh.write(svg)
            messagebox.showinfo('Saved', f'SVG saved:\n{path}')
        except Exception as e:
            messagebox.showerror('Export Error', str(e))

    def _build_png_wrapped(self, seqs, consensus, variability,
                            pos_labels=None, uniprot_name=None,
                            master_name=None, start_aa=1, master_seq=None):
        SCALE = 3
        CW = 14 * SCALE; CH = 16 * SCALE; LW = 200 * SCALE
        RH = 20 * SCALE; BLOCK_GAP = 8 * SCALE; TOP_PAD = 8 * SCALE
        COLS_PER_BLOCK = 75
        length = len(consensus)
        if pos_labels is None:
            pos_labels = [str(i + start_aa) for i in range(length)]

        top_label = master_name or 'CONSENSUS'
        all_rows   = [(top_label, consensus)] + list(seqs)
        n_rpb      = len(all_rows)
        block_h    = n_rpb * CH
        n_blocks   = (length + COLS_PER_BLOCK - 1) // COLS_PER_BLOCK
        img_w      = LW + COLS_PER_BLOCK * CW
        img_h      = TOP_PAD + n_blocks * (RH + block_h + BLOCK_GAP)

        img  = _PILImage.new('RGB', (img_w, img_h), _hex_to_rgb(BG))
        draw = _PILDraw.Draw(img)
        font_seq   = _load_mono_font(10 * SCALE)
        font_lbl   = _load_mono_font(9  * SCALE)
        font_ruler = _load_mono_font(8  * SCALE)

        col_colors = [_hex_to_rgb(self._var_to_color(variability[i]))
                      for i in range(length)]

        for block_i in range(n_blocks):
            start_col = block_i * COLS_PER_BLOCK
            end_col   = min(start_col + COLS_PER_BLOCK, length)
            block_y   = TOP_PAD + block_i * (RH + block_h + BLOCK_GAP)

            draw.rectangle([0, block_y, img_w, block_y + RH],
                            fill=_hex_to_rgb(PANEL))
            for col in range(start_col, end_col):
                label = pos_labels[col] if col < len(pos_labels) else str(col + 1)
                try:
                    pos_num = int(label)
                    show = (col == start_col or pos_num % 10 == 0)
                except ValueError:
                    show = False
                if show:
                    xc = LW + (col - start_col) * CW + CW // 2
                    draw.text((max(0, xc - 10 * SCALE), block_y + 4 * SCALE),
                               label, fill=_hex_to_rgb(DIM_FG), font=font_ruler)

            for row_i, (name, seq) in enumerate(all_rows):
                y0   = block_y + RH + row_i * CH
                is_c = (row_i == 0)
                lbl_bg = _hex_to_rgb(PANEL if is_c else BG)
                lbl_fg = _hex_to_rgb(ACCENT if is_c else DIM_FG)
                draw.rectangle([0, y0, LW, y0 + CH], fill=lbl_bg)
                disp = (name[:24] + '…') if len(name) > 25 else name
                draw.text((6 * SCALE, y0 + 3 * SCALE), disp,
                           fill=lbl_fg, font=font_lbl)

                for col in range(start_col, end_col):
                    if col >= len(seq): continue
                    x0  = LW + (col - start_col) * CW
                    ch  = seq[col]
                    var = variability[col] if col < len(variability) else 100.0
                    bg  = col_colors[col]
                    con = consensus[col] if col < len(consensus) else ''
                    conservation = 1.0 - var / 100.0
                    txt_fg = (255, 255, 255) if conservation > 0.60 else (26, 26, 26)
                    master_char = master_seq[col] if (master_seq and col < len(master_seq)) else None
                    is_insertion = (master_char == '-') if master_char is not None else False
                    if is_c:
                        if is_insertion and ch == '-':
                            cell_bg = _hex_to_rgb(INS_BG); cell_fg = _hex_to_rgb(INS_FG)
                        else:
                            cell_bg = bg; cell_fg = txt_fg
                    elif is_insertion:
                        if ch in ('-', 'N', 'X'):
                            cell_bg = _hex_to_rgb(BG); cell_fg = _hex_to_rgb(BORDER)
                        else:
                            cell_bg = _hex_to_rgb(INS_BG); cell_fg = _hex_to_rgb(INS_FG)
                    elif ch in ('-', 'N', 'X'):
                        cell_bg = _hex_to_rgb(BG)
                        cell_fg = _hex_to_rgb(BORDER)
                    elif ch == con:
                        cell_bg = bg; cell_fg = txt_fg
                    else:
                        cell_bg = _hex_to_rgb('#FFE8D8')
                        cell_fg = _hex_to_rgb('#A33000')
                    draw.rectangle([x0, y0, x0 + CW, y0 + CH], fill=cell_bg)
                    draw.text((x0 + 2 * SCALE, y0 + 3 * SCALE), ch,
                               fill=cell_fg, font=font_seq)

            div_y = block_y + RH + CH
            draw.line([(0, div_y), (img_w, div_y)],
                      fill=_hex_to_rgb(ACCENT), width=SCALE)

        return img

    def _build_svg_wrapped(self, seqs, consensus, variability,
                            pos_labels=None, uniprot_name=None,
                            master_name=None, start_aa=1, master_seq=None):
        """Build a scalable SVG alignment image styled like the reference."""
        FONT_SIZE    = 13
        CW           = 7.8        # Courier New char width ≈ 0.6 × font-size
        CH           = 18         # row height (line height)
        BASELINE     = 13         # text baseline from row top
        LW           = 220.0      # label column width (~28 chars)
        RULER_H      = 22         # ruler row height
        RULER_BASE   = 15         # ruler text baseline from ruler top
        BLOCK_GAP    = 42         # vertical gap between wrapped blocks
        TOP_PAD      = 14         # padding above first block
        COLS_PER_ROW = 75

        def esc(s):
            return (str(s)
                    .replace('&', '&amp;')
                    .replace('<', '&lt;')
                    .replace('>', '&gt;'))

        length = len(consensus)
        # Normalise pos_labels: generate or rebase so they start at start_aa
        if pos_labels is None:
            pos_labels = [str(i + start_aa) for i in range(length)]
        elif start_aa != 1:
            # Skip insertion labels like '0+1'; find the first plain integer
            _first = next(
                (int(p) for p in pos_labels if str(p).lstrip('-').isdigit()),
                None)
            if _first is not None:
                _shift = start_aa - _first
                pos_labels = [
                    str(int(p) + _shift) if str(p).lstrip('-').isdigit() else p
                    for p in pos_labels
                ]

        top_label  = master_name or 'CONSENSUS'
        all_rows   = [(top_label, consensus)] + list(seqs)
        n_rows     = len(all_rows)
        col_classes = _classify_columns(all_rows, length)

        n_blocks  = (length + COLS_PER_ROW - 1) // COLS_PER_ROW
        block_h   = n_rows * CH
        img_w     = LW + COLS_PER_ROW * CW
        img_h     = TOP_PAD + n_blocks * (RULER_H + block_h + BLOCK_GAP)

        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{img_w:.1f}" height="{img_h:.1f}">',
            f'<style>text{{font-family:"Courier New",Courier,monospace;'
            f'font-size:{FONT_SIZE}px}}</style>',
            f'<rect width="{img_w:.1f}" height="{img_h:.1f}" fill="{BG}"/>',
        ]

        for block_i in range(n_blocks):
            start_col = block_i * COLS_PER_ROW
            end_col   = min(start_col + COLS_PER_ROW, length)
            by        = TOP_PAD + block_i * (RULER_H + block_h + BLOCK_GAP)

            # Ruler background
            parts.append(f'<rect x="0" y="{by:.1f}" width="{img_w:.1f}" '
                         f'height="{RULER_H}" fill="{PANEL}"/>')

            # Ruler numbers — centered over their column
            # Pre-find how many cols until the first 10-multiple in this block
            _first_ten_offset = None
            for _c in range(start_col, end_col):
                _l = pos_labels[_c] if _c < len(pos_labels) else str(_c + 1)
                try:
                    if int(_l) % 10 == 0:
                        _first_ten_offset = _c - start_col
                        break
                except ValueError:
                    pass

            for col in range(start_col, end_col):
                lbl = pos_labels[col] if col < len(pos_labels) else str(col + 1)
                try:
                    n = int(lbl)
                    if col == start_col:
                        if n % 10 == 0:
                            # Start col is itself a 10-multiple — always show it
                            show = True
                        else:
                            # Only show if the nearest 10-multiple is far enough
                            # away to avoid overlap
                            show = (_first_ten_offset is None
                                    or _first_ten_offset >= 5)
                    else:
                        show = (n % 10 == 0)
                except ValueError:
                    show = False
                if show:
                    xc = LW + (col - start_col + 0.5) * CW
                    # Number label
                    parts.append(
                        f'<text x="{xc:.1f}" y="{by + RULER_BASE:.1f}" '
                        f'fill="{DIM_FG}" font-size="{FONT_SIZE - 2}px" '
                        f'text-anchor="middle">{esc(lbl)}</text>')
                    # Tick mark at the bottom of the ruler row
                    parts.append(
                        f'<line x1="{xc:.1f}" y1="{by + RULER_H - 4:.1f}" '
                        f'x2="{xc:.1f}" y2="{by + RULER_H:.1f}" '
                        f'stroke="{DIM_FG}" stroke-width="1"/>')

            # Sequence rows
            for row_i, (name, seq) in enumerate(all_rows):
                ry   = by + RULER_H + row_i * CH
                is_c = (row_i == 0)
                lbg  = PANEL if is_c else BG
                lfg  = ACCENT if is_c else DIM_FG

                # Label cell
                parts.append(f'<rect x="0" y="{ry:.1f}" width="{LW:.1f}" '
                             f'height="{CH}" fill="{lbg}"/>')
                disp_name = (name[:27] + '\u2026') if len(name) > 28 else name
                parts.append(f'<text x="5" y="{ry + BASELINE:.1f}" '
                             f'fill="{lfg}">{esc(disp_name)}</text>')

                # Characters
                for col in range(start_col, end_col):
                    if col >= len(seq):
                        continue
                    ch     = seq[col]
                    mc     = master_seq[col] if (master_seq and col < len(master_seq)) else None
                    is_ins = (mc == '-')     if mc is not None else False
                    cc     = col_classes[col]
                    disp   = '.' if ch == '-' else ch

                    x    = LW + (col - start_col) * CW
                    xmid = x + CW / 2

                    if ch == '-':
                        cell_bg = None; cell_fg = BORDER
                    elif cc == 'identical':
                        cell_bg = ACCENT; cell_fg = '#FFFFFF'
                    elif cc == 'similar':
                        cell_bg = None; cell_fg = ACCENT
                    else:
                        cell_bg = None; cell_fg = FG

                    if cell_bg:
                        parts.append(
                            f'<rect x="{x:.1f}" y="{ry:.1f}" '
                            f'width="{CW:.1f}" height="{CH}" fill="{cell_bg}"/>')
                    parts.append(
                        f'<text x="{xmid:.1f}" y="{ry + BASELINE:.1f}" '
                        f'text-anchor="middle" fill="{cell_fg}">{esc(disp)}</text>')

            # Blue box borders around runs of similar-group (non-identical) columns
            i = start_col
            while i < end_col:
                if col_classes[i] == 'similar':
                    j = i
                    while j < end_col and col_classes[j] == 'similar':
                        j += 1
                    bx = LW + (i - start_col) * CW
                    bw = min((j - i) * CW, img_w - bx)
                    parts.append(
                        f'<rect x="{bx:.1f}" y="{by + RULER_H:.1f}" '
                        f'width="{bw:.1f}" height="{block_h:.1f}" '
                        f'fill="none" stroke="#2255AA" stroke-width="1"/>')
                    i = j
                else:
                    i += 1

            # Thin accent divider after consensus/master row (drawn last = on top)
            div_y = by + RULER_H + CH
            parts.append(
                f'<line x1="0" y1="{div_y:.1f}" x2="{img_w:.1f}" y2="{div_y:.1f}" '
                f'stroke="{ACCENT}" stroke-width="0.5"/>')

        parts.append('</svg>')
        return '\n'.join(parts)

    # ════════════════════════════════════════════════════════════════════════
    #  COMPARE VIEW — loading
    # ════════════════════════════════════════════════════════════════════════

    def _browse_cmp_master(self):
        path = filedialog.askopenfilename(
            title='Select Master file',
            filetypes=[('Sequence files', '*.gpsat *.vgat *.fa *.fasta *.faa *.fna *.txt'),
                       ('All Files', '*.*')])
        if not path: return
        name = simpledialog.askstring(
            'Master Name',
            'Enter a display name for the master sequence:',
            initialvalue=os.path.splitext(os.path.basename(path))[0])
        if name is None: return
        self._cmp_master_path_str = path
        self._cmp_master_user_name = name.strip() or os.path.splitext(
            os.path.basename(path))[0]
        self._cmp_master_display.config(
            text=f'{self._cmp_master_user_name}  ({os.path.basename(path)})')

    def _add_cmp_student(self):
        path = filedialog.askopenfilename(
            title='Select Student file',
            filetypes=[('Sequence files', '*.gpsat *.vgat *.fa *.fasta *.faa *.fna *.txt'),
                       ('All Files', '*.*')])
        if not path: return
        name = simpledialog.askstring(
            'Student Name',
            'Enter a display name for this student sequence:',
            initialvalue=os.path.splitext(os.path.basename(path))[0])
        if name is None: return
        display_name = name.strip() or os.path.splitext(os.path.basename(path))[0]
        # Avoid duplicates
        if path not in [p for _, p, _ in self._cmp_student_paths]:
            self._cmp_student_paths.append((display_name, path, display_name))
            self._cmp_student_lb.insert(
                tk.END,
                f'{display_name}  ({os.path.basename(path)})')

    def _remove_cmp_student(self):
        sel = list(self._cmp_student_lb.curselection())
        for i in reversed(sel):
            self._cmp_student_lb.delete(i)
            del self._cmp_student_paths[i]

    def _run_compare(self):
        master_path = getattr(self, '_cmp_master_path_str', '').strip()
        if not master_path or not os.path.isfile(master_path):
            messagebox.showerror('Error', 'Select a master file first.')
            return
        if not self._cmp_student_paths:
            messagebox.showerror('Error', 'Add at least one student file.')
            return
        if self._mafft_busy:
            return

        master_name = getattr(self, '_cmp_master_user_name',
                              os.path.splitext(os.path.basename(master_path))[0])
        try:
            master_seq = _load_cmp_seq(master_path)
        except Exception as e:
            messagebox.showerror('Error', f'Could not read master file:\n{e}')
            return

        student_seqs = []
        for display_name, spath, user_name in self._cmp_student_paths:
            try:
                student_seqs.append((user_name, _load_cmp_seq(spath)))
            except Exception as e:
                messagebox.showerror('Error', f'Could not read {spath}:\n{e}')
                return

        all_seqs = [(master_name, master_seq)] + student_seqs
        self._mafft_busy = True
        self._status.config(
            text=f'Aligning {len(all_seqs)} consensus sequences via MAFFT …')

        threading.Thread(target=self._compare_mafft_worker,
                         args=(all_seqs, master_name), daemon=True).start()

    def _compare_mafft_worker(self, all_seqs, master_name):
        try:
            aligned = _run_mafft(all_seqs)
            length  = min(len(s) for _, s in aligned)
            variability, consensus = _cross_variability(aligned, length)
            master_aligned = next(
                (s for n, s in aligned if n == master_name), aligned[0][1])
            master_nums = _build_master_numbering(master_aligned)
            self.after(0, self._compare_mafft_done,
                       aligned, variability, consensus, length,
                       master_name, master_nums)
        except Exception as exc:
            self.after(0, self._cov_mafft_error, str(exc))

    def _compare_mafft_done(self, aligned, variability, consensus, length,
                             master_name, master_nums):
        self._mafft_busy      = False
        self._cmp_seqs        = aligned
        self._cmp_variability = variability
        self._cmp_length      = length
        self._cmp_master_name = master_name
        self._cmp_master_nums = master_nums
        master_aligned = next(
            (s for n, s in aligned if n == master_name), aligned[0][1])
        self._cmp_consensus = master_aligned
        self._switch_view('compare')
        self._status.config(
            text=f'Compare: master "{master_name}" vs. '
                 f'{len(aligned) - 1} student(s) · {length} columns')

    # ════════════════════════════════════════════════════════════════════════
    #  DRAWING — COMPARE VIEW
    # ════════════════════════════════════════════════════════════════════════

    def _draw_compare(self):
        if self._cmp_length == 0:
            return
        CW  = self.CMP_CW;  CH  = self.CMP_CH
        LW  = self.LW;       RH  = self.RULER_H
        BLOCK_GAP = 10

        canvas_w = self._cmp_cv.winfo_width()
        if canvas_w <= 1: canvas_w = 700
        cols_per_block = max(1, (canvas_w - LW) // CW)
        n_pos  = self._cmp_length
        master_seq = self._cmp_consensus

        student_rows = [(n, s) for n, s in self._cmp_seqs
                        if n != self._cmp_master_name]
        all_rows = [(self._cmp_master_name or 'MASTER', master_seq)] + student_rows
        n_rpb    = len(all_rows)
        block_h  = n_rpb * CH
        n_blocks = (n_pos + cols_per_block - 1) // cols_per_block
        total_h  = n_blocks * (RH + block_h + BLOCK_GAP)

        self._cmp_cv.delete('all')
        self._cmp_cv.config(scrollregion=(0, 0, canvas_w, total_h))
        self._cmp_layout = {
            'cols_per_block': cols_per_block, 'n_all_rows': n_rpb,
            'block_h': block_h, 'block_gap': BLOCK_GAP,
            'ruler_h': RH, 'CW': CW, 'CH': CH, 'LW': LW, 'n_pos': n_pos}

        try:
            cmp_start_aa = max(1, int(self._cmp_start_aa_var.get()))
        except (ValueError, AttributeError):
            cmp_start_aa = 1
        nums = list(self._cmp_master_nums) if self._cmp_master_nums else \
               [str(i + 1) for i in range(n_pos)]
        if cmp_start_aa != 1:
            # Skip insertion labels like '0+1'; find the first plain integer
            _first = next(
                (int(n) for n in nums if str(n).lstrip('-').isdigit()), None)
            if _first is not None:
                _shift = cmp_start_aa - _first
                nums = [str(int(n) + _shift) if str(n).lstrip('-').isdigit()
                        else n for n in nums]
        col_colors = [self._var_to_color(v)
                      for v in self._cmp_variability[:n_pos]]

        for block_i in range(n_blocks):
            start_col = block_i * cols_per_block
            end_col   = min(start_col + cols_per_block, n_pos)
            block_y   = block_i * (RH + block_h + BLOCK_GAP)

            self._cmp_cv.create_rectangle(0, block_y, canvas_w,
                                           block_y + RH, fill=PANEL, outline='')
            self._cmp_cv.create_text(
                6, block_y + RH // 2, text='Pos (master)',
                anchor='w', fill=DIM_FG,
                font=('Helvetica', _fs(7), 'italic'))
            for col in range(start_col, end_col):
                label = nums[col] if col < len(nums) else str(col + 1)
                try:
                    pos_num = int(label)
                    show = (col == start_col or pos_num % 10 == 0)
                except ValueError:
                    show = False
                if show:
                    xc = LW + (col - start_col) * CW + CW // 2
                    self._cmp_cv.create_text(
                        xc, block_y + RH // 2, text=label,
                        anchor='center', fill=DIM_FG,
                        font=('Helvetica', _fs(7)))
                    self._cmp_cv.create_line(
                        xc, block_y + RH - 3, xc, block_y + RH,
                        fill=BORDER, width=1)
            self._cmp_cv.create_line(0, block_y + RH, canvas_w,
                                      block_y + RH, fill=BORDER, width=1)

            for row_i, (name, seq) in enumerate(all_rows):
                y0   = block_y + RH + row_i * CH
                yc   = y0 + CH // 2
                is_m = (row_i == 0)

                lbl_bg = PANEL if is_m else BG
                lbl_fg = ACCENT if is_m else DIM_FG
                self._cmp_cv.create_rectangle(0, y0, LW, y0 + CH,
                                               fill=lbl_bg, outline='')
                disp = (name[:20] + '…') if len(name) > 21 else name
                self._cmp_cv.create_text(
                    6, yc, text=disp, anchor='w', fill=lbl_fg,
                    font=('Courier', _fs(8), 'bold' if is_m else 'normal'))

                for col in range(start_col, end_col):
                    if col >= len(seq): continue
                    x0  = LW + (col - start_col) * CW
                    ch  = seq[col]
                    var = self._cmp_variability[col] if col < len(self._cmp_variability) else 100.0
                    bg  = col_colors[col]
                    master_char = master_seq[col] if col < len(master_seq) else '-'
                    is_insertion = (master_char == '-')
                    conservation = 1.0 - var / 100.0
                    txt_fg = '#ffffff' if conservation > 0.60 else '#1A1A1A'

                    if is_m:
                        if ch == '-':
                            cell_bg, cell_fg = INS_BG, INS_FG
                        else:
                            cell_bg, cell_fg = bg, txt_fg
                    elif is_insertion:
                        if ch in ('-', 'N', 'X'):
                            cell_bg, cell_fg = BG, BORDER
                        else:
                            cell_bg, cell_fg = INS_BG, INS_FG
                    elif ch in ('-', 'N', 'X'):
                        cell_bg, cell_fg = BG, BORDER
                    elif ch == master_char:
                        cell_bg, cell_fg = bg, txt_fg
                    else:
                        cell_bg, cell_fg = '#FFE8D8', '#A33000'

                    self._cmp_cv.create_rectangle(
                        x0, y0, x0 + CW, y0 + CH, fill=cell_bg, outline='')
                    self._cmp_cv.create_text(
                        x0 + CW // 2, yc, text=ch, anchor='center',
                        fill=cell_fg, font=('Courier', _fs(7)))

            self._cmp_cv.create_line(
                0, block_y + RH + CH, canvas_w, block_y + RH + CH,
                fill=ACCENT, width=1)

    def _on_cmp_resize(self, event):
        if hasattr(self, '_cmp_resize_job'):
            self.after_cancel(self._cmp_resize_job)
        self._cmp_resize_job = self.after(
            80, lambda: self._draw_compare()
            if self._cmp_length > 0 and self._view_mode.get() == 'compare'
            else None)

    def _on_cmp_wheel(self, event):
        if event.num == 4:
            units = -1
        elif event.num == 5:
            units = 1
        elif event.delta:
            units = -1 if event.delta > 0 else 1
        else:
            return
        self._cmp_cv.yview_scroll(units, 'units')

    def _on_cmp_hover(self, event):
        layout = self._cmp_layout
        if not layout: return
        cpb   = layout['cols_per_block']
        bh    = layout['block_h']; bg_gap = layout['block_gap']
        RH    = layout['ruler_h']
        CW    = layout['CW']; LW = layout['LW']
        n_pos = layout['n_pos']
        cx = self._cmp_cv.canvasx(event.x)
        cy = self._cmp_cv.canvasy(event.y)
        bt  = RH + bh + bg_gap
        bi  = int(cy / bt)
        yib = cy - bi * bt
        if yib < RH or yib >= RH + bh or cx < LW:
            self._cmp_tip.place_forget(); return
        cib = int((cx - LW) / CW)
        col = bi * cpb + cib
        if col >= n_pos or cib >= cpb:
            self._cmp_tip.place_forget(); return
        var   = self._cmp_variability[col] if col < len(self._cmp_variability) else 0.0
        label = self._cmp_master_nums[col]  if col < len(self._cmp_master_nums) else str(col + 1)
        master_char = self._cmp_consensus[col] if col < len(self._cmp_consensus) else '-'
        lines = [f'Master position : {label}',
                 f'Master residue  : {master_char}',
                 f'Cross-variability : {var:.1f}%']
        if master_char == '-':
            lines.append('(Insertion — not in master)')
        student_rows = [(n, s) for n, s in self._cmp_seqs
                        if n != self._cmp_master_name]
        for sname, sseq in student_rows[:8]:
            ch = sseq[col] if col < len(sseq) else '-'
            lines.append(f'  {sname[:18]}: {ch}')
        if len(student_rows) > 8:
            lines.append(f'  … +{len(student_rows) - 8} more')
        self._cmp_tip.config(text='\n'.join(lines))
        self._cmp_tip.place(x=event.x + 14, y=event.y + 14)

    # ════════════════════════════════════════════════════════════════════════
    #  CUSTOM SEQUENCES — .gpsat FACTORY
    # ════════════════════════════════════════════════════════════════════════

    def _custom_add_from_entry(self):
        name = self._cust_name_var.get().strip()
        raw  = self._cust_paste_text.get('1.0', tk.END).strip()
        seq  = raw.upper().replace(' ', '').replace('\n', '').replace('\r', '')
        if not name:
            messagebox.showerror('Error', 'Enter a sequence name.'); return
        if not seq:
            messagebox.showerror('Error', 'Paste a sequence first.'); return
        # Strip FASTA header if pasted
        if seq.startswith('>'):
            lines = seq.split('\n')
            seq = ''.join(l for l in lines if not l.startswith('>')).strip()
        valid = set('ACDEFGHIKLMNPQRSTVWYATGCUN-')
        bad   = [c for c in seq if c not in valid]
        if bad:
            messagebox.showerror(
                'Invalid characters',
                f'Sequence contains unrecognised characters: '
                f'{", ".join(set(bad))}\n'
                'Use standard amino-acid or nucleotide letters.')
            return
        self._custom_entries.append((name, seq))
        preview = seq[:30] + ('…' if len(seq) > 30 else '')
        self._cust_lb.insert(tk.END, f'{name}  ({len(seq)} aa/nt)  {preview}')
        self._cust_name_var.set('')
        self._cust_paste_text.delete('1.0', tk.END)

    def _custom_upload_fasta(self):
        path = filedialog.askopenfilename(
            title='Upload FASTA',
            filetypes=[('FASTA Files', '*.fa *.fasta *.faa *.txt'),
                       ('All Files', '*.*')])
        if not path: return
        try:
            entries, _ = parse_fasta_entries(input_path=path)
        except Exception as e:
            messagebox.showerror('Error', str(e)); return
        added = 0
        for e in entries:
            name = e.get('protein_name') or e['id']
            seq  = e['sequence']
            self._custom_entries.append((name, seq))
            preview = seq[:30] + ('…' if len(seq) > 30 else '')
            self._cust_lb.insert(
                tk.END, f'{name}  ({len(seq)} aa/nt)  {preview}')
            added += 1
        self._cust_status.config(text=f'{added} sequence(s) added from FASTA.')

    def _custom_upload_gpsat(self):
        path = filedialog.askopenfilename(
            title='Upload .gpsat file',
            filetypes=[('GPSAT file', '*.gpsat *.vgat'), ('All Files', '*.*')])
        if not path: return
        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception as e:
            messagebox.showerror('Error', f'Could not read file:\n{e}'); return
        if data.get('type') not in ('gpsat_consensus', 'vgat_consensus'):
            messagebox.showerror('Error', 'Not a valid .gpsat consensus file.')
            return
        name = data.get('name') or os.path.splitext(os.path.basename(path))[0]
        seq  = data['consensus'].replace('-', '')
        self._custom_entries.append((name, seq))
        preview = seq[:30] + ('…' if len(seq) > 30 else '')
        self._cust_lb.insert(
            tk.END, f'{name}  ({len(seq)} aa/nt)  {preview}')
        self._cust_status.config(text=f'Added: {name}')

    def _custom_upload_other(self):
        path = filedialog.askopenfilename(
            title='Upload sequence file',
            filetypes=[('All Supported', '*.fa *.fasta *.txt *.seq *.gpsat *.vgat'),
                       ('All Files', '*.*')])
        if not path: return
        ext = os.path.splitext(path)[1].lower()
        if ext in ('.gpsat', '.vgat'):
            self._custom_upload_gpsat_path(path)
            return
        # Try as FASTA
        try:
            entries, _ = parse_fasta_entries(input_path=path)
            if entries:
                added = 0
                for e in entries:
                    name = e.get('protein_name') or e['id']
                    seq  = e['sequence']
                    self._custom_entries.append((name, seq))
                    preview = seq[:30] + ('…' if len(seq) > 30 else '')
                    self._cust_lb.insert(
                        tk.END, f'{name}  ({len(seq)} aa/nt)  {preview}')
                    added += 1
                self._cust_status.config(
                    text=f'{added} sequence(s) added.')
                return
        except Exception:
            pass
        # Try as plain sequence text
        try:
            with open(path) as fh:
                raw = fh.read().strip().upper().replace(' ', '').replace('\n', '')
            name = os.path.splitext(os.path.basename(path))[0]
            self._custom_entries.append((name, raw))
            preview = raw[:30] + ('…' if len(raw) > 30 else '')
            self._cust_lb.insert(
                tk.END, f'{name}  ({len(raw)} aa/nt)  {preview}')
            self._cust_status.config(text=f'Added: {name}')
        except Exception as e:
            messagebox.showerror('Error', f'Could not parse file:\n{e}')

    def _custom_upload_gpsat_path(self, path):
        try:
            with open(path) as fh:
                data = json.load(fh)
        except Exception as e:
            messagebox.showerror('Error', str(e)); return
        name = data.get('name') or os.path.splitext(os.path.basename(path))[0]
        seq  = data.get('consensus', '').replace('-', '')
        self._custom_entries.append((name, seq))
        preview = seq[:30] + ('…' if len(seq) > 30 else '')
        self._cust_lb.insert(
            tk.END, f'{name}  ({len(seq)} aa/nt)  {preview}')
        self._cust_status.config(text=f'Added: {name}')

    def _custom_remove(self):
        sel = list(self._cust_lb.curselection())
        for i in reversed(sel):
            self._cust_lb.delete(i)
            del self._custom_entries[i]

    def _custom_clear(self):
        if not self._custom_entries:
            return
        if messagebox.askyesno('Clear', 'Remove all sequences from the pool?'):
            self._custom_entries.clear()
            self._cust_lb.delete(0, tk.END)
            self._cust_status.config(text='')

    def _custom_align_export(self):
        if len(self._custom_entries) < 2:
            messagebox.showerror('Too few sequences',
                                 'Add at least 2 sequences before aligning.')
            return
        if self._mafft_busy:
            return
        path = filedialog.asksaveasfilename(
            title='Save custom alignment as .gpsat',
            defaultextension='.gpsat',
            filetypes=[('GPSAT file', '*.gpsat'), ('All Files', '*.*')])
        if not path: return

        self._mafft_busy = True
        self._cust_status.config(
            text=f'Aligning {len(self._custom_entries)} sequences via MAFFT…')
        seqs = list(self._custom_entries)
        file_name = os.path.splitext(os.path.basename(path))[0]

        def worker():
            try:
                aligned = _run_mafft(seqs)
                consensus, variability, counts_data = _votes_in_memory(aligned)
                a2, c2, v2, cd2 = _remove_gap_consensus_cols(
                    aligned, consensus, variability, counts_data)
                self.after(0, lambda: self._custom_export_done(
                    path, file_name, c2, v2, cd2))
            except Exception as exc:
                self.after(0, lambda m=str(exc): self._custom_export_error(m))

        threading.Thread(target=worker, daemon=True).start()

    def _custom_export_done(self, path, file_name, consensus, variability, counts_data):
        self._mafft_busy = False
        try:
            _write_gpsat(path, file_name, consensus, variability, counts_data)
            self._cust_status.config(text=f'Saved: {os.path.basename(path)}')
            messagebox.showinfo('Saved',
                                f'.gpsat file saved:\n{path}\n\n'
                                f'Consensus length: {len(consensus)} positions\n'
                                f'Load this file in the Compare view as a student '
                                f'or master sequence.')
        except Exception as e:
            messagebox.showerror('Error', str(e))

    def _custom_export_error(self, msg):
        self._mafft_busy = False
        self._cust_status.config(text='Error — see dialog')
        messagebox.showerror('MAFFT Error', msg)

    # ── Shared MAFFT error handler ────────────────────────────────────────────

    def _cov_mafft_error(self, msg):
        self._mafft_busy      = False
        self._uniprot_pending = False
        self._status.config(text=f'MAFFT error: {msg}')


# ─── LOGOMAKER HELPER ─────────────────────────────────────────────────────────

def _render_logo(seqs, consensus, max_cols=200, start_pos=None, end_pos=None):
    """Build a logomaker figure from aligned sequences.
    start_pos / end_pos are 1-indexed, inclusive (None = auto).
    Returns a matplotlib Figure.  Raises if logomaker unavailable."""
    if not _LOGO_AVAILABLE:
        raise RuntimeError('logomaker not installed')

    length = min(len(s) for _, s in seqs) if seqs else len(consensus)

    if start_pos is not None or end_pos is not None:
        # User-specified range (1-indexed, clamp to valid bounds)
        trim_start = max(0, (start_pos - 1) if start_pos is not None else 0)
        trim_end   = min(length, end_pos if end_pos is not None else length)
        if trim_start >= trim_end:
            trim_start, trim_end = 0, length
    elif length > max_cols:
        # Auto: pick the most variable window
        var = []
        for i in range(length):
            residues = [s[i] for _, s in seqs if i < len(s) and s[i] != '-']
            if residues:
                counts = {r: residues.count(r) for r in set(residues)}
                top = max(counts.values())
                var.append((top - len(residues)) / len(residues) * -1)
            else:
                var.append(0.0)
        best_start = 0
        best_score = float('-inf')
        for i in range(length - max_cols + 1):
            score = sum(var[i:i + max_cols])
            if score > best_score:
                best_score = score
                best_start = i
        trim_start = best_start
        trim_end   = best_start + max_cols
    else:
        trim_start = 0
        trim_end   = length

    # Determine sequence alphabet
    all_chars = set()
    for _, s in seqs:
        all_chars.update(s[trim_start:trim_end].replace('-', ''))
    nucl_chars = set('ATGCU')
    aa_chars   = set('ACDEFGHIKLMNPQRSTVWY')
    is_nucl    = len(all_chars - nucl_chars) == 0 and len(all_chars) > 0

    alphabet   = list(nucl_chars if is_nucl else aa_chars)

    # Build counts matrix
    n_pos = trim_end - trim_start
    rows  = []
    for col_i in range(n_pos):
        col = trim_start + col_i
        counts = {a: 0 for a in alphabet}
        for _, s in seqs:
            if col < len(s):
                ch = s[col]
                if ch in counts:
                    counts[ch] += 1
        rows.append(counts)

    df = pd.DataFrame(rows, columns=alphabet).fillna(0)
    # Normalise to information content
    df = df.div(df.sum(axis=1).replace(0, 1), axis=0)

    fig, ax = plt.subplots(figsize=(min(20, n_pos * 0.25 + 2), 2.5))
    try:
        logomaker.Logo(df, ax=ax,
                       color_scheme='chemistry' if not is_nucl else 'classic',
                       vpad=0.05, width=0.9)
    except Exception:
        logomaker.Logo(df, ax=ax)
    ax.set_xlabel('Position', fontsize=8)
    ax.set_ylabel('Frequency', fontsize=8)
    ax.tick_params(labelsize=7)
    # Set x ticks to actual positions
    tick_positions = list(range(0, n_pos, max(1, n_pos // 10)))
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([str(trim_start + p + 1) for p in tick_positions],
                        fontsize=7)
    fig.tight_layout()
    return fig
