import tkinter as tk
from tkinter import ttk, filedialog
import functions as func
import sys
import io

selected_proteins = []
# Function to open a new GUI window for Picking Protein Types

def open_picking_protein_types():


    new_window = tk.Toplevel(root)
    new_window.title("Picking Protein Types")
    new_window.geometry("800x600")  # Increased width and height for better visibility
    
    label = tk.Label(new_window, text="Picking Protein Types Analysis", font=("Times New Roman", 14, "bold"))
    label.pack(pady=10)
    
    def select_ncbi_file():
        file_path = filedialog.askopenfilename(title="Select NCBI Input File")
        ncbi_entry.delete(0, tk.END)
        ncbi_entry.insert(0, file_path)
    
    def select_allpros_folder():
        folder_path = filedialog.askdirectory(title="Select Folder for allpros.txt")
        allpros_folder_entry.delete(0, tk.END)
        allpros_folder_entry.insert(0, folder_path)
    
    def select_keylengths_folder():
        folder_path = filedialog.askdirectory(title="Select Folder for KeyLengths")
        keylengths_folder_entry.delete(0, tk.END)
        keylengths_folder_entry.insert(0, folder_path)
    
    input_frame = tk.Frame(new_window)
    input_frame.pack(pady=10)
    
    ncbi_label = tk.Label(input_frame, text="NCBI Input File:")
    ncbi_label.grid(row=0, column=0)
    ncbi_entry = tk.Entry(input_frame, width=40)
    ncbi_entry.grid(row=0, column=1)
    ncbi_button = ttk.Button(input_frame, text="Browse", command=select_ncbi_file)
    ncbi_button.grid(row=0, column=2)
    
    allpros_label = tk.Label(input_frame, text="Select Folder for allpros.txt:")
    allpros_label.grid(row=1, column=0)
    allpros_folder_entry = tk.Entry(input_frame, width=40)
    allpros_folder_entry.grid(row=1, column=1)
    allpros_folder_button = ttk.Button(input_frame, text="Browse", command=select_allpros_folder)
    allpros_folder_button.grid(row=1, column=2)
    
    allpros_name_label = tk.Label(input_frame, text="Enter File Name for allpros.txt:")
    allpros_name_label.grid(row=2, column=0)
    allpros_name_entry = tk.Entry(input_frame, width=40)
    allpros_name_entry.grid(row=2, column=1, columnspan=2)
    
    keylengths_folder_label = tk.Label(input_frame, text="KeyLengths Folder:")
    keylengths_folder_label.grid(row=3, column=0)
    keylengths_folder_entry = tk.Entry(input_frame, width=40)
    keylengths_folder_entry.grid(row=3, column=1)
    keylengths_folder_button = ttk.Button(input_frame, text="Browse", command=select_keylengths_folder)
    keylengths_folder_button.grid(row=3, column=2)
    
    keylengths_name_label = tk.Label(input_frame, text="Key Lengths File Name:")
    keylengths_name_label.grid(row=4, column=0)
    keylengths_name_entry = tk.Entry(input_frame, width=40)
    keylengths_name_entry.grid(row=4, column=1, columnspan=2)
    
    display_frame = tk.Frame(new_window)
    display_frame.pack(pady=10)
    
    output_text = tk.Text(display_frame, height=10, width=40, state=tk.DISABLED)
    output_text.grid(row=0, column=0, padx=5)
    
    protein_listbox = tk.Listbox(display_frame, selectmode=tk.MULTIPLE, width=40, height=10)
    protein_listbox.grid(row=0, column=1, padx=5)
    
    def compute():
        file_name = ncbi_entry.get()
        all_pros_folder = allpros_folder_entry.get()
        all_pros_file_name = allpros_name_entry.get()
        all_pros_file_path = f"{all_pros_folder}/{all_pros_file_name}"
        
        output_text.config(state=tk.NORMAL)
        output_text.delete(1.0, tk.END)
        sys.stdout = io.StringIO()
        
        output = func.determining_protein_types(file_name, all_pros_file_path)
        print('By looking at the protein types above, decide which protein types you want to explore more.')
        
        output_text.insert(tk.END, sys.stdout.getvalue())
        output_text.config(state=tk.DISABLED)
        sys.stdout = sys.__stdout__
        
        protein_listbox.delete(0, tk.END)
        for item in output:
            protein_listbox.insert(tk.END, item)
    
    def select_proteins():
        selected_proteins = [protein_listbox.get(i) for i in protein_listbox.curselection()]
        
        keylengths_folder = keylengths_folder_entry.get()
        keylengths_file_name = keylengths_name_entry.get()
        output_file = f"{keylengths_folder}/{keylengths_file_name}"
        
        all_pros_file_path = f"{allpros_folder_entry.get()}/{allpros_name_entry.get()}"
        func.key_protein_lengths(selected_proteins, output_file, all_pros_file_path)
        print(' ')
        print('Your file has been created in your request pathway.')   
    
    button_frame = tk.Frame(new_window)
    button_frame.pack(pady=10)
    
    compute_button = ttk.Button(button_frame, text="Compute", command=compute)
    compute_button.grid(row=0, column=0, padx=5)
    
    select_button = ttk.Button(button_frame, text="Select Proteins", command=select_proteins)
    select_button.grid(row=0, column=1, padx=5)
    
    close_button = ttk.Button(button_frame, text="Close", command=new_window.destroy)
    close_button.grid(row=0, column=2, padx=5)

def open_longest_protein_window():
    new_window = tk.Toplevel(root)
    new_window.title("Determine Longest Protein")
    new_window.geometry("800x600")  # Increased width and height for better visibility
    
    input_frame = tk.Frame(new_window)
    input_frame.pack(side=tk.LEFT, padx=10, pady=10)
    
    output_frame = tk.Frame(new_window)
    output_frame.pack(side=tk.RIGHT, padx=10, pady=10)
    
    def select_file():
        file_path = filedialog.askopenfilename(title="Select Input File")
        file_entry.delete(0, tk.END)
        file_entry.insert(0, file_path)
    
    def select_output_folder():
        folder_path = filedialog.askdirectory(title="Select Output Folder")
        folder_entry.delete(0, tk.END)
        folder_entry.insert(0, folder_path)
    
    def select_consensus_input_file():
        file_path = filedialog.askopenfilename(title="Select Consensus Input File")
        consensus_input_entry.delete(0, tk.END)
        consensus_input_entry.insert(0, file_path)
    
    def select_consensus_output_folder():
        folder_path = filedialog.askdirectory(title="Select Consensus Output Folder")
        consensus_output_folder_entry.delete(0, tk.END)
        consensus_output_folder_entry.insert(0, folder_path)
    
    # Max Lengths Section
    tk.Label(input_frame, text="Select Input File:").pack(pady=5)
    file_entry = tk.Entry(input_frame, width=40)
    file_entry.pack()
    tk.Button(input_frame, text="Browse", command=select_file).pack()
    
    tk.Label(input_frame, text="Enter Genotypes (comma-separated):").pack(pady=5)
    genotype_entry = tk.Entry(input_frame, width=40)
    genotype_entry.pack()
    
    tk.Label(input_frame, text="Select Output Folder:").pack(pady=5)
    folder_entry = tk.Entry(input_frame, width=40)
    folder_entry.pack()
    tk.Button(input_frame, text="Browse", command=select_output_folder).pack()
    
    def run_analysis():
        last_file_name = file_entry.get()
        genotypes = genotype_entry.get().split(',')
        folder_name = folder_entry.get()
        func.finding_longest_protein(last_file_name, genotypes, folder_name)
    
    tk.Button(input_frame, text="Run Analysis", command=run_analysis).pack(pady=10)
    
    # Consensus Determination Section
    tk.Label(output_frame, text="Select Consensus Input File:").pack(pady=5)
    consensus_input_entry = tk.Entry(output_frame, width=40)
    consensus_input_entry.pack()
    tk.Button(output_frame, text="Browse", command=select_consensus_input_file).pack()
    
    tk.Label(output_frame, text="Select Consensus Output Folder:").pack(pady=5)
    consensus_output_folder_entry = tk.Entry(output_frame, width=40)
    consensus_output_folder_entry.pack()
    tk.Button(output_frame, text="Browse", command=select_consensus_output_folder).pack()
    
    tk.Label(output_frame, text="Enter Consensus Output File Name:").pack(pady=5)
    consensus_output_name_entry = tk.Entry(output_frame, width=40)
    consensus_output_name_entry.pack()
    
    def generate_consensus_file():
        input_file = consensus_input_entry.get()
        output_folder = consensus_output_folder_entry.get()
        output_file = f"{output_folder}/{consensus_output_name_entry.get()}"
        func.consensus_determination_file_creator(input_file, output_file)
    
    tk.Button(output_frame, text="Generate Consensus Determination File", command=generate_consensus_file).pack(pady=10)
    

def open_votes_calculation_window():
    new_window = tk.Toplevel(root)
    new_window.title("Votes Calculation")
    new_window.geometry("500x300")
    
    def select_input_file():
        file_path = filedialog.askopenfilename(title="Select Alignment Input File")
        input_entry.delete(0, tk.END)
        input_entry.insert(0, file_path)
    
    def select_output_folder():
        folder_path = filedialog.askdirectory(title="Select Output Folder")
        output_folder_entry.delete(0, tk.END)
        output_folder_entry.insert(0, folder_path)
    
    tk.Label(new_window, text="Input Alignment File:").pack()
    input_entry = tk.Entry(new_window, width=40)
    input_entry.pack()
    tk.Button(new_window, text="Browse", command=select_input_file).pack()
    
    tk.Label(new_window, text="Output Folder:").pack()
    output_folder_entry = tk.Entry(new_window, width=40)
    output_folder_entry.pack()
    tk.Button(new_window, text="Browse", command=select_output_folder).pack()
    
    tk.Label(new_window, text="Output File Name:").pack()
    output_file_entry = tk.Entry(new_window, width=40)
    output_file_entry.pack()
    
    tk.Button(new_window, text="Calculate Votes", command=lambda: func.protein_votes_file(input_entry.get(), f"{output_folder_entry.get()}/{output_file_entry.get()}")).pack()

def open_alignment_graph_window():
    new_window = tk.Toplevel(root)
    new_window.title("Create Alignment Graph")
    new_window.geometry("500x300")
    
    def select_votes_file():
        file_path = filedialog.askopenfilename(title="Select Votes File")
        votes_entry.delete(0, tk.END)
        votes_entry.insert(0, file_path)
    
    tk.Label(new_window, text="Select Votes File:").pack()
    votes_entry = tk.Entry(new_window, width=40)
    votes_entry.pack()
    tk.Button(new_window, text="Browse", command=select_votes_file).pack()
    
    tk.Label(new_window, text="Enter Graph Title:").pack()
    title_entry = tk.Entry(new_window, width=40)
    title_entry.pack()
    
    def generate_graph():
        votes_file = votes_entry.get()
        title = title_entry.get()
        func.create_variability_graph(votes_file, title)
    
    tk.Button(new_window, text="Generate Graph", command=generate_graph).pack(pady=10)

def open_consensus_printer_window():
    new_window = tk.Toplevel(root)
    new_window.title("Consensus Sequence Printer")
    new_window.geometry("500x300")
    
    def select_votes_file():
        file_path = filedialog.askopenfilename(title="Select Votes File")
        votes_entry.delete(0, tk.END)
        votes_entry.insert(0, file_path)
    
    def run_consensus():
        file_name = votes_entry.get()
        output_capture = io.StringIO()  # Create a StringIO buffer
        sys.stdout = output_capture  # Redirect stdout

        try:
            func.get_consensus(file_name)  # Call the function
            output_text.config(state=tk.NORMAL)
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, output_capture.getvalue())  # Insert the captured output
            output_text.config(state=tk.DISABLED)
        except Exception as e:
            output_text.config(state=tk.NORMAL)
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, f"Error: {str(e)}")
            output_text.config(state=tk.DISABLED)
        finally:
            sys.stdout = sys.__stdout__  # Restore original stdout



    
    tk.Label(new_window, text="Select Votes File:").pack()
    votes_entry = tk.Entry(new_window, width=40)
    votes_entry.pack()
    tk.Button(new_window, text="Browse", command=select_votes_file).pack()
    
    tk.Button(new_window, text="Run Consensus", command=run_consensus).pack(pady=10)
    
    output_text = tk.Text(new_window, height=10, width=50, state=tk.DISABLED)
    output_text.pack()    

# Initialize main window
root = tk.Tk()
root.title("Viral Genome Analysis Tool")

# Set the initial dimensions to fit all widgets
root.geometry("1560x600")  # Increased width and height for better visibility

# Configure root for dynamic resizing
root.rowconfigure(0, weight=1)
root.columnconfigure(0, weight=1)

# Main Frame to Hold Title, Description, and Options
main_frame = tk.Frame(root)
main_frame.grid(row=0, column=0, columnspan=3, padx=5, pady=5, sticky="nsew")

# Configure main_frame for resizing
main_frame.rowconfigure(0, weight=1)  # Title row
main_frame.rowconfigure(1, weight=3)  # Description row (increased weight for more space)
main_frame.columnconfigure(0, weight=2)  # Title and description column
main_frame.columnconfigure(1, weight=1)  # Options column

# Title Section
title_frame = tk.Frame(main_frame, relief="groove", borderwidth=3)
title_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
title_label = tk.Label(title_frame, text="Viral Genome Analysis Tool", font=("Times New Roman", 14, "bold"))
title_label.pack(padx=5, pady=5)

# Description Section
description_frame = tk.Frame(main_frame, relief="groove", borderwidth=3)
description_frame.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
description_label = tk.Label(
    description_frame,
    text=(
        "This software allows users to parse through and analyze specific viral proteins for consensus sequence analysis, "
        "ORF determination, protein translation as well as graphical analysis. Please refer to the documentation listed "
        "here to download a PDF version of documentation for a step-by-step process. For a self-guided experience, please "
        "refer to the links below for analysis."
    ),
    wraplength=600,  # Increased wraplength to fit the frame
    justify="left",  # Align text to the left
    font=("Times New Roman", 12)
)
description_label.pack(padx=10, pady=10, fill="both", expand=True)

# Options Section (Next to Title and Description)
options_frame = tk.Frame(main_frame, relief="groove", borderwidth=3)
options_frame.grid(row=0, column=1, rowspan=2, padx=10, pady=5, sticky="nsew")
options_label = tk.Label(options_frame, text="Options", font=("Times New Roman", 12))
options_label.grid(row=0, column=0, columnspan=4, pady=5)

option_buttons = {
    "Picking Protein Types": open_picking_protein_types,
    "Determining Longest Protein": lambda: open_analysis_window("Determining Longest Protein"),
    "Votes Calculation": open_votes_calculation_window,
    "Creating Alignment graph": lambda: open_analysis_window("Creating Alignment graph"),
    "ORF match file creator": lambda: open_analysis_window("ORF match file creator"),
    "ORF Analyzer": lambda: open_analysis_window("ORF Analyzer"),
    "Genome Modifier": lambda: open_analysis_window("Genome Modifier"),
    "Protein Translator": lambda: open_analysis_window("Protein Translator"),
    "Consensus Sequence Printer": lambda: open_analysis_window("Consensus Sequence Printer")
}

# Update "Determining Longest Protein" button to open the new window
option_buttons["Determining Longest Protein"] = open_longest_protein_window
option_buttons["Creating Alignment graph"] = open_alignment_graph_window
option_buttons["Consensus Sequence Printer"] = open_consensus_printer_window

# Configure options_frame for dynamic resizing
options_frame.rowconfigure(0, weight=1)  # Title row
for i in range(3):  # Configure rows for buttons
    options_frame.rowconfigure(i + 1, weight=1)
for j in range(3):  # Configure columns for buttons
    options_frame.columnconfigure(j, weight=1)

# Arrange buttons in a grid (3 rows, 3 columns) with explicit width
button_grid = [
    list(option_buttons.keys())[:3],  # First row (3 buttons)
    list(option_buttons.keys())[3:6],  # Second row (3 buttons)
    list(option_buttons.keys())[6:]   # Third row (remaining buttons)
]

for r, row in enumerate(button_grid):
    for c, option in enumerate(row):
        btn = ttk.Button(options_frame, text=option, command=option_buttons[option], width=25)  # Set explicit width
        btn.grid(row=r + 1, column=c, padx=10, pady=10, sticky="nsew")

# Run Application
root.mainloop()
