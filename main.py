import tkinter as tk
from tkinter import ttk
import functions2o as f2

# ─── THEME SETUP ─────────────────────────────────────────────────────────────

BG       = f2.BG
PANEL    = f2.PANEL
ENTRY_BG = f2.ENTRY_BG
FG       = f2.FG
DIM_FG   = f2.DIM_FG
ACCENT   = f2.ACCENT
SEL_FG   = f2.SEL_FG
BORDER   = f2.BORDER


def setup_theme(root):
    sz  = f2._font_size
    sz9 = max(8, sz - 1)

    root.configure(bg=BG)
    style = ttk.Style(root)
    style.theme_use('clam')

    # ── Base elements ──────────────────────────────────────────────────────
    style.configure('.',
                    background=BG,
                    foreground=FG,
                    font=('Helvetica', sz),
                    borderwidth=0,
                    relief='flat')

    style.configure('TFrame', background=BG)
    style.configure('Panel.TFrame', background=PANEL)

    # ── Labels ────────────────────────────────────────────────────────────
    style.configure('TLabel',        background=BG,    foreground=FG)
    style.configure('Panel.TLabel',  background=PANEL, foreground=FG)
    style.configure('Section.TLabel',background=PANEL, foreground=FG,
                    font=('Helvetica', sz, 'bold'))
    style.configure('Dim.TLabel',    background=PANEL, foreground=DIM_FG,
                    font=('Helvetica', sz9))

    # ── Buttons ───────────────────────────────────────────────────────────
    style.configure('TButton',
                    background=ENTRY_BG,
                    foreground=FG,
                    padding=(10, 5),
                    borderwidth=1,
                    relief='flat')
    style.map('TButton',
              background=[('active', BORDER), ('pressed', BORDER)],
              foreground=[('active', FG)])

    style.configure('Accent.TButton',
                    background=ACCENT,
                    foreground=SEL_FG,
                    padding=(12, 6),
                    font=('Helvetica', sz, 'bold'))
    style.map('Accent.TButton',
              background=[('active', '#D96000'), ('pressed', '#A34A00')],
              foreground=[('active', SEL_FG)])

    # ── Entry ─────────────────────────────────────────────────────────────
    style.configure('TEntry',
                    fieldbackground=ENTRY_BG,
                    foreground=FG,
                    insertcolor=FG,
                    bordercolor=BORDER,
                    lightcolor=BORDER,
                    darkcolor=BORDER)

    # ── Scrollbar ─────────────────────────────────────────────────────────
    style.configure('TScrollbar',
                    background=ENTRY_BG,
                    troughcolor=BG,
                    arrowcolor=DIM_FG,
                    bordercolor=BG,
                    gripcount=0)
    style.map('TScrollbar',
              background=[('active', BORDER)])

    # ── Notebook (tabs) ───────────────────────────────────────────────────
    style.configure('TNotebook',
                    background=BG,
                    borderwidth=0,
                    tabmargins=[0, 0, 0, 0])
    style.configure('TNotebook.Tab',
                    background=ENTRY_BG,
                    foreground=DIM_FG,
                    padding=[16, 8],
                    borderwidth=0,
                    font=('Helvetica', sz))
    style.map('TNotebook.Tab',
              background=[('selected', PANEL), ('active', BORDER)],
              foreground=[('selected', ACCENT), ('active', FG)])

    # ── Checkbutton ───────────────────────────────────────────────────────
    style.configure('TCheckbutton',       background=PANEL, foreground=FG,
                    font=('Helvetica', sz9))
    style.configure('Panel.TCheckbutton', background=PANEL, foreground=FG,
                    font=('Helvetica', sz9))
    style.map('TCheckbutton',
              background=[('active', PANEL)], foreground=[('active', FG)])
    style.map('Panel.TCheckbutton',
              background=[('active', PANEL)], foreground=[('active', FG)])

    # ── Separator ─────────────────────────────────────────────────────────
    style.configure('TSeparator', background=BORDER)


# ─── MAIN APPLICATION ─────────────────────────────────────────────────────────

def main():
    root = tk.Tk()
    root.title('Genomic & Proteomic Sequence Analysis Toolkit')
    root.geometry('720x600')
    root.minsize(600, 480)

    setup_theme(root)

    # ── Header banner ─────────────────────────────────────────────────────
    banner = tk.Frame(root, bg=PANEL, pady=10)
    banner.pack(fill='x')

    title_lbl = tk.Label(banner,
                          text='Genomic & Proteomic Sequence Analysis Toolkit',
                          bg=PANEL, fg=ACCENT,
                          font=('Helvetica', 15, 'bold'))
    title_lbl.pack()
    subtitle_lbl = tk.Label(banner,
                             text='Prepare Data  ·  Alignment Viewer'
                                  '  (Alignment · Compare · Custom)',
                             bg=PANEL, fg=DIM_FG,
                             font=('Helvetica', 9))
    subtitle_lbl.pack()

    ttk.Separator(root, orient='horizontal').pack(fill='x')

    # ── Font size slider bar ───────────────────────────────────────────────
    font_bar = tk.Frame(root, bg=PANEL, pady=4)
    font_bar.pack(fill='x', padx=14)
    tk.Label(font_bar, text='Font Size:', bg=PANEL, fg=DIM_FG,
             font=('Helvetica', 9)).pack(side='left')
    _font_size_var = tk.IntVar(value=f2._font_size)
    _font_val_lbl  = tk.Label(font_bar, text=f'{f2._font_size} pt',
                               bg=PANEL, fg=FG,
                               font=('Helvetica', 9), width=5, anchor='w')
    _font_val_lbl.pack(side='left', padx=(6, 0))
    _font_scale_widget = ttk.Scale(font_bar, from_=7, to=18, orient='horizontal',
                                   variable=_font_size_var, length=160)
    _font_scale_widget.pack(side='left', padx=(6, 0))
    ttk.Separator(root, orient='horizontal').pack(fill='x')

    # ── Notebook ──────────────────────────────────────────────────────────
    nb = ttk.Notebook(root)
    nb.pack(fill='both', expand=True, padx=0, pady=0)

    tab1 = f2.ProteinFilterTab(nb)
    tab2 = f2.AlignmentViewerTab(nb)

    nb.add(tab1, text='  1 · Prepare Data  ')
    nb.add(tab2, text='  2 · Alignment Viewer  ')

    # ── Status bar ────────────────────────────────────────────────────────
    status_bar = tk.Frame(root, bg=ENTRY_BG, height=24)
    status_bar.pack(fill='x', side='bottom')
    sb_left = tk.Label(status_bar, text='Ready',
                        bg=ENTRY_BG, fg=DIM_FG,
                        font=('Helvetica', 9), anchor='w')
    sb_left.pack(side='left', padx=10)
    sb_right = tk.Label(status_bar, text='KeLab  ·  GPSAT',
                         bg=ENTRY_BG, fg=DIM_FG,
                         font=('Helvetica', 9), anchor='e')
    sb_right.pack(side='right', padx=10)

    # ── Font change callback (defined after all widget refs exist) ────────
    def _on_font_change(value):
        sz = int(float(value))
        f2._font_size = sz
        _font_val_lbl.config(text=f'{sz} pt')
        setup_theme(root)
        f2.apply_font_scale()
        # Update the fixed tk.Label widgets in the banner and status bar
        title_lbl.config(font=('Helvetica', max(11, sz + 5), 'bold'))
        subtitle_lbl.config(font=('Helvetica', max(7, sz - 1)))
        sb_left.config(font=('Helvetica', max(7, sz - 1)))
        sb_right.config(font=('Helvetica', max(7, sz - 1)))
        # Redraw active canvas view in the alignment tab
        tab2._redraw_current()

    _font_scale_widget.config(command=_on_font_change)

    root.mainloop()


if __name__ == '__main__':
    main()
