import tkinter as tk

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("My App")

        # Create the button and set its command to the say_hello method
        self.process_button = tk.Button(root, text="Process", command=self.say_hello)
        self.process_button.grid(row=0, column=0, padx=10, pady=10)

    def say_hello(self):
        print("hej")

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()