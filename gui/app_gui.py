
import tkinter as tk
from tkinter import filedialog, messagebox
from utils.file_utils import write_dna_file, read_dna_file
from dna.compressor import compress, decompress
from dna.translator import get_complementary_sequence

def process_dna(input_data):
    # Save input to sample_dna.txt
    write_dna_file("data/sample_dna.txt", input_data)

    # Process input
    complementary_dna = get_complementary_sequence(input_data)
    compressed = compress(complementary_dna)
    decompressed = decompress(compressed, len(input_data))

    return complementary_dna, compressed, decompressed

def browse_file(entry):
    filepath = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
    if filepath:
        content = read_dna_file(filepath)
        entry.delete(0, tk.END)
        entry.insert(0, content)

def run_gui():
    def on_process():
        dna_input = input_entry.get().strip().upper()

        if not all(base in 'ATCG' for base in dna_input):
            messagebox.showerror("Invalid Input", "DNA sequence must contain only A, T, C, G.")
            return

        complementary, compressed, decompressed = process_dna(dna_input)

        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, f"Complementary: {complementary}\n\n")
        result_text.insert(tk.END, f"Compressed: {list(compressed)}\n\n")
        result_text.insert(tk.END, f"Decompressed: {decompressed}\n")

    window = tk.Tk()
    window.title("DNA Mutation Analyzer")
    window.geometry("700x500")

    tk.Label(window, text="Enter DNA Sequence or Upload File:", font=("Arial", 12)).pack(pady=10)

    input_frame = tk.Frame(window)
    input_frame.pack()

    input_entry = tk.Entry(input_frame, width=60, font=("Arial", 12))
    input_entry.pack(side=tk.LEFT, padx=5)

    browse_btn = tk.Button(input_frame, text="📂 Browse File", command=lambda: browse_file(input_entry))
    browse_btn.pack(side=tk.LEFT)

    process_btn = tk.Button(window, text="🧬 Process DNA", font=("Arial", 12), command=on_process)
    process_btn.pack(pady=20)

    result_text = tk.Text(window, height=15, wrap=tk.WORD, font=("Consolas", 11))
    result_text.pack(padx=10, fill=tk.BOTH, expand=True)

    window.mainloop()
