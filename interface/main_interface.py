import tkinter as tk
from tkinter import filedialog, messagebox

# Create the main window
root = tk.Tk()
root.title("Genomic Data Processor")
root.geometry("600x650")

# Default values
DEFAULTS = {
    "database": 'databases/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24.fasta',
    "matrix": 'databases/BLOSUM62.fasta',  # Default BLOSUM62 matrix file
    "parallel": True,
    "match": 3,
    "mismatch": -1,
    "gap_open": -10,
    "gap_extend": -4,
    "process": True
}

# Function to select a file
def select_file(entry_field):
    file_path = filedialog.askopenfilename()
    entry_field.delete(0, tk.END)
    entry_field.insert(0, file_path)

# Function to select a directory
def select_directory(entry_field):
    directory = filedialog.askdirectory()
    entry_field.delete(0, tk.END)
    entry_field.insert(0, directory)

# Function to process the input
def process_inputs():
    inputs = {
        "fasta_file": fasta_file_entry.get(),
        "output_path": output_path_entry.get(),
        "gene": gene_entry.get(),
        "database": database_entry.get(),
        "matrix": matrix_entry.get(),
        "parallel": parallel_var.get(),
        "match": int(match_entry.get()),
        "mismatch": int(mismatch_entry.get()),
        "gap_open": int(gap_open_entry.get()),
        "gap_extend": int(gap_extend_entry.get()),
        "process": process_var.get()
    }
    messagebox.showinfo("Inputs", f"Processed inputs: {inputs}")

# Create labels, entry fields, and buttons for each argument
fields = [
    ("Fasta File", "Path to the fasta file with the whole genome for the strain"),
    ("Output Path", "Path to where to save the output files"),
    ("Organism name", "Name of the organism to process"),
    ("Database", "Path to the chromoprotein database", DEFAULTS["database"]),
    ("BLOSUM62 Matrix", "Path to the BLOSUM62 matrix file", DEFAULTS["matrix"]),
    ("Match Score", "Score for a match", DEFAULTS["match"]),
    ("Mismatch Penalty", "Penalty for a mismatch", DEFAULTS["mismatch"]),
    ("Gap Open Penalty", "Gap opening penalty", DEFAULTS["gap_open"]),
    ("Gap Extend Penalty", "Gap extension penalty", DEFAULTS["gap_extend"])
]

entries = {}

for idx, (label_text, help_text, *default) in enumerate(fields):
    label = tk.Label(root, text=label_text, anchor='w')
    label.grid(row=idx, column=0, padx=10, pady=5, sticky='e')
    entry = tk.Entry(root, width=50, justify='center')
    if default:
        entry.insert(0, default[0])
    entry.grid(row=idx, column=1, padx=10, pady=5, sticky='ew')
    entries[label_text] = entry
    if label_text in ["Fasta File", "Output Path", "Organism name", "Database", "BLOSUM62 Matrix"]:
        button = tk.Button(root, text="Browse", command=lambda e=entry: select_file(e) if label_text != "Output Path" else select_directory(e))
        button.grid(row=idx, column=2, padx=10, pady=5, sticky='ew')

# Special entries
fasta_file_entry = entries["Fasta File"]
output_path_entry = entries["Output Path"]
gene_entry = entries["Organism name"]
database_entry = entries["Database"]
matrix_entry = entries["BLOSUM62 Matrix"]
match_entry = entries["Match Score"]
mismatch_entry = entries["Mismatch Penalty"]
gap_open_entry = entries["Gap Open Penalty"]
gap_extend_entry = entries["Gap Extend Penalty"]

# Create checkboxes for boolean arguments
process_var = tk.BooleanVar(value=DEFAULTS["process"])
process_checkbox = tk.Checkbutton(root, text="Enable DNA to Protein Processing", variable=process_var)
process_checkbox.grid(row=9, column=0, columnspan=3, padx=10, pady=5, sticky='w')

# Create the process button
process_button = tk.Button(root, text="Process", command=process_inputs)
process_button.grid(row=10, column=0, columnspan=3, padx=10, pady=20, sticky='ew')

# Configure column weights to center-align elements
for col in range(3):
    root.grid_columnconfigure(col, weight=1)

# Configure row weights to center-align elements
root.grid_rowconfigure(list(range(10)), weight=1)

# Start the Tkinter event loop
root.mainloop()
