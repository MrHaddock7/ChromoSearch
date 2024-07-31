if __name__ == "__main__":
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog, messagebox
    import logging
    import threading

    import sys
    import os

    root = tk.Tk()
    root.title("Chromoprotein Genome Processor")
    root.geometry("700x750")
    root.resizable(False, False)

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler("interface_project.log")]
    )

    logger = logging.getLogger(__name__)   

    # Add the directory containing the module to the system path
    sys.path.append(os.path.abspath('/Users/klonk/Desktop/Chromoprotein_Strikes_Back/ChromoSearch'))
    sys.path.append(os.path.abspath('/Users/klonk/Desktop/Chromoprotein_Strikes_Back/ChromoSearch/databases/chromoproteins_uniprot'))
    from chromosearch import main as ChromoSearch

    # Create the main window

    logger.info('Started a main window')

    # Default values
    DEFAULTS = {
        "fasta_file": '/Users/klonk/Desktop/genomes/k12.fasta',
        "output_path": '/Users/klonk/Desktop/genomes',
        "database": '/Users/klonk/Desktop/Chromoprotein_Strikes_Back/ChromoSearch/databases/chromoproteins_uniprot/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24',
        # "matrix": True,  # Default BLOSUM62 matrix file
        "parallel": True,
        "match": 3,
        "mismatch": -1,
        "gap_open": -10,
        "gap_extend": -4,
        "process": True,
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

    def execute_threading():
        thread = threading.Thread(target=process_inputs)
        thread.start()

    # Function to process the input
    def process_inputs():
        inputs = {
            "fasta_file": fasta_file_entry.get(),
            "output_path": output_path_entry.get(),
            "gene": gene_entry.get(),
            "database": database_entry,
            # "matrix": matrix_entry.get(),
            "match": int(match_entry.get()),
            "mismatch": int(mismatch_entry.get()),
            "gap_open": int(gap_open_entry.get()),
            "gap_extend": int(gap_extend_entry.get()),
            "process": process_var.get()
        }
        # messagebox.showinfo("Processing", 'Your genome is being processed')

        logger = logging.getLogger(__name__)

        ChromoSearch(fasta_path=inputs["fasta_file"], 
                    output_path=inputs["output_path"], 
                    gene=inputs["gene"],
                    database=inputs["database"],
                    blastpnsw=True,
                    save_intermediates=True,
                    process=inputs["process"]
                    )

    # Create labels, entry fields, and buttons for each argument
    fields = [
        ("Fasta File", "Path to the fasta file with the whole genome for the strain", DEFAULTS["fasta_file"]), 
        ("Output Path", "Path to where to save the output files", DEFAULTS["output_path"]),
        ("Organism name", "Name of the organism to process"),
        ("Database", "Path to the chromoprotein database", DEFAULTS["database"]),
        # ("BLOSUM62 Matrix", "Path to the BLOSUM62 matrix file", DEFAULTS["matrix"]),
        ("Match Score", "Score for a match", DEFAULTS["match"]),
        ("Mismatch Penalty", "Penalty for a mismatch", DEFAULTS["mismatch"]),
        ("Gap Open Penalty", "Gap opening penalty", DEFAULTS["gap_open"]),
        ("Gap Extend Penalty", "Gap extension penalty", DEFAULTS["gap_extend"])
    ]

    entries = {}

    for idx, (label_text, help_text, *default) in enumerate(fields):
        label = tk.Label(root, text=label_text, anchor='w')
        label.grid(row=idx, column=0, padx=10, pady=5, sticky='e')
        
        if label_text == "Database":
            def on_option_selected(value):
                print(f"Selected: {value}")

            selected_option = tk.StringVar()
            selected_option.set("Chromoproteins")  # Set default value
            options = ["Chromoproteins", "Pigment-pathway enzymes"]

            dropdown = tk.OptionMenu(root, selected_option, *options, command=on_option_selected)
            dropdown.grid(row=idx, column=1, padx=10, pady=5, sticky='ew')  # Grid instead of pack
            if selected_option.get() == "Chromoproteins":
                entries[label_text] = "databases/chromoproteins_uniprot/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24"  # Store the StringVar instead of the dropdown widget
            elif selected_option.get() == "Pigment-pathway enzymes":
                entries[label_text] = "databases/pigment_biosynthesis_chromoproteins/uniprotkb_go_manual_0046148_NOT_taxonom_2024_07_01"
            # else:
            #     raise ValueError("Invalid argument for the database")


        else:
            entry = tk.Entry(root, width=50, justify='center')
            if default:
                entry.insert(0, default[0])
            entry.grid(row=idx, column=1, padx=10, pady=5, sticky='ew')
            entries[label_text] = entry
            
            if label_text in ["Fasta File", "Output Path", "Organism name"]:
                button = tk.Button(root, text="Browse", command=lambda e=entry: select_file(e) if label_text != "Output Path" else select_directory(e))
                button.grid(row=idx, column=2, padx=10, pady=5, sticky='ew')

    class RedirectText(object):
        def __init__(self, text_widget):
            self.text_widget = text_widget

        def write(self, string):
            self.text_widget.insert(tk.END, string)
            self.text_widget.see(tk.END)  # Scroll to the end

        def flush(self):
            pass
    

    # Create a Text widget for displaying terminal output
    text_widget = tk.Text(root, wrap='word', height=10, width=80)
    text_widget.grid(row=len(fields), column=0, columnspan=3, padx=10, pady=5, sticky='nsew')

    # Redirect stdout to the Text widget
    sys.stdout = RedirectText(text_widget)

    # Special entries
    fasta_file_entry = entries["Fasta File"]
    output_path_entry = entries["Output Path"]
    gene_entry = entries["Organism name"]
    database_entry = entries["Database"]
    # matrix_entry = entries["BLOSUM62 Matrix"]
    match_entry = entries["Match Score"]
    mismatch_entry = entries["Mismatch Penalty"]
    gap_open_entry = entries["Gap Open Penalty"]
    gap_extend_entry = entries["Gap Extend Penalty"]

    # Create checkboxes for boolean arguments
    process_var = tk.BooleanVar(value=DEFAULTS["process"])
    process_checkbox = tk.Checkbutton(root, text="Enable DNA to Protein Processing", variable=process_var)
    process_checkbox.grid(row=9, column=0, columnspan=3, padx=10, pady=5, sticky='w')

    # Create the process button
    process_button = tk.Button(root, text="Process", command=execute_threading)
    process_button.grid(row=10, column=0, columnspan=3, padx=10, pady=20, sticky='ew')

    # Configure column weights to center-align elements
    for col in range(3):
        root.grid_columnconfigure(col, weight=1)

    # Configure row weights to center-align elements
    root.grid_rowconfigure(list(range(10)), weight=1)

    # Start the Tkinter event loop
    root.mainloop()
