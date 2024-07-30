import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import logging
import threading
import sys
import os
import threading

sys.path.append(os.path.abspath('/Users/klonk/Desktop/Chromoprotein_Strikes_Back/ChromoSearch'))
sys.path.append(os.path.abspath('/Users/klonk/Desktop/Chromoprotein_Strikes_Back/ChromoSearch/databases'))
from chromosearch import main as ChromoSearch

def select_file(self, entry_field):
    file_path = filedialog.askopenfilename()
    self.entry_field.delete(0, tk.END)
    self.entry_field.insert(0, file_path)

# Function to select a directory
def select_directory(self, entry_field):
    directory = filedialog.askdirectory()
    self.entry_field.delete(0, tk.END)
    entry_field.insert(0, directory)

class RedirectText(object):
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END)  # Scroll to the end

    def flush(self):
        pass

class Main_app:

    class Arguments:
        def __init__(self):
            self.fastafile = '/Users/klonk/Desktop/genomes/k12.fasta'
            self.outputfile = '/Users/klonk/Desktop/genomes'
            self.database = '/Users/klonk/Desktop/Chromoprotein_Strikes_Back/ChromoSearch/databases/chromoproteins_uniprot/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24'
            self.parallel = True
            self.match = 3
            self.mismatch = -1
            self.gap_open_entry = -10
            self.gap_extend = -4
            self.process = True

    def __init__(self, root):

        self.exe_threading = False
        root.title("Chromoprotein Genome Processor")
        root.geometry("700x750")
        # root.resizable(False, False)
        self.arguments = self.Arguments()

        ## Input arguments for chromosearch.py

        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler("interface_project.log")]
        )

        logger = logging.getLogger(__name__)  

        logger.info('Started a main window')

        self._fields = [
            ("Fasta File", "Path to the fasta file with the whole genome for the strain", self.arguments.fastafile), 
            ("Output Path", "Path to where to save the output files", self.arguments.outputfile),
            ("Organism name", "Name of the organism to process"),
            ("Database", "Path to the chromoprotein database", self.arguments.database),
            # ("BLOSUM62 Matrix", "Path to the BLOSUM62 matrix file", DEFAULTS["matrix"]),
            ("Match Score", "Score for a match", self.arguments.match),
            ("Mismatch Penalty", "Penalty for a mismatch", self.arguments.mismatch),
            ("Gap Open Penalty", "Gap opening penalty", self.arguments.gap_open_entry),
            ("Gap Extend Penalty", "Gap extension penalty", self.arguments.gap_extend)
        ]

        self._entries = {}

        for idx, (label_text, help_text, *default) in enumerate(self._fields):
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
                    self._entries[label_text] = "databases/chromoproteins_uniprot/uniprotkb_chromophore_keyword_KW_0157_AND_reviewed_2024_06_24"  # Store the StringVar instead of the dropdown widget
                elif selected_option.get() == "Pigment-pathway enzymes":
                    self._entries[label_text] = "databases/pigment_biosynthesis_chromoproteins/uniprotkb_go_manual_0046148_NOT_taxonom_2024_07_01"
                # else:
                #     raise ValueError("Invalid argument for the database")


            else:
                entry = tk.Entry(root, width=50, justify='center')
                if default:
                    entry.insert(0, default[0])
                entry.grid(row=idx, column=1, padx=10, pady=5, sticky='ew')
                self._entries[label_text] = entry
                
                if label_text in ["Fasta File", "Output Path", "Organism name"]:
                    button = tk.Button(root, text="Browse", command=lambda e=entry: select_file(e) if label_text != "Output Path" else select_directory(e))
                    button.grid(row=idx, column=2, padx=10, pady=5, sticky='ew')

        text_widget = tk.Text(root, wrap='word', height=10, width=80)
        text_widget.grid(row=len(self._fields), column=0, columnspan=3, padx=10, pady=5, sticky='nsew')

        process_var = tk.BooleanVar(value=self.arguments.process)
        process_checkbox = tk.Checkbutton(root, text="Enable DNA to Protein Processing")
        process_checkbox.grid(row=9, column=0, columnspan=3, padx=10, pady=5, sticky='w')

        # Create the process button
        process_button = tk.Button(root, text="Process", command=self.execute_threading)
        process_button.grid(row=10, column=0, columnspan=3, padx=10, pady=20, sticky='ew')

        sys.stdout = RedirectText(text_widget)
        
    def execute_chromosearch(self):
        print('exec')
        ChromoSearch(fasta_path=self._entries["Fasta File"].get(), 
                    output_path=self._entries["Output Path"].get(), 
                    gene=self._entries["Organism name"].get(),
                    database=self._entries["Database"],
                    blastpnsw=True,
                    save_intermediates=True,
                    process=self.arguments.process
                    )

    def execute_threading(self):
        thread = threading.Thread(target=self.execute_chromosearch())
        thread.start()
    
if __name__ == '__main__':
    root = tk.Tk()
    app = Main_app(root)
    tk.mainloop()

        




