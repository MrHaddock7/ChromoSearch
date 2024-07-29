import tkinter as tk
from tkinter import ttk

def on_option_selected(value):
    print(f"Selected: {value}")

# Create the main window
root = tk.Tk()
root.title("Dropdown Menu Example")

# Create a StringVar to hold the value of the selected option
selected_option = tk.StringVar()
selected_option.set("Option 1")  # Set default value

# Define the list of options
options = ["Option 1", "Option 2", "Option 3", "Option 4"]

# Create the OptionMenu
dropdown = ttk.OptionMenu(root, selected_option, options[0], *options, command=on_option_selected)

# Pack the OptionMenu into the window
dropdown.pack(pady=10)

# Run the application
root.mainloop()
