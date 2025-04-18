import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dna.compressor import compress, decompress
from dna.translator import get_complementary_sequence
from dna.mutation_detector import get_result
from dna.mutation_detector_2 import get_result2
from utils.file_utils import read_dna_file, write_dna_file

class DNASequencingApp(tk.Tk):
    def __init__(self):
        super().__init__()
        
        self.title("DNA Sequencing Toolkit")
        self.geometry("900x650")
        self.minsize(900, 650)
        
        # Configure style
        self.style = ttk.Style()
        self.style.configure("TFrame", background="#f0f0f0")
        self.style.configure("TButton", font=("Arial", 11))
        self.style.configure("TLabel", font=("Arial", 11), background="#f0f0f0")
        self.style.configure("Header.TLabel", font=("Arial", 14, "bold"), background="#f0f0f0")
        
        # Main container
        self.main_container = ttk.Frame(self)
        self.main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create the menu frame and content frame
        self.create_menu_frame()
        self.create_content_frame()
        
        # Initialize with home content
        self.show_home()
        
        # For storing current mutation detection type
        self.current_mutation_type = tk.StringVar(value="tp53")
    
    def create_menu_frame(self):
        """Create the menu sidebar"""
        self.menu_frame = ttk.Frame(self.main_container, width=200, style="TFrame")
        self.menu_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # App title
        ttk.Label(self.menu_frame, text="DNA Sequencing", 
                  font=("Arial", 16, "bold"), style="Header.TLabel").pack(pady=(0, 20))
        
        # Menu buttons
        menu_options = [
            ("Home", self.show_home),
            ("Compress DNA", self.show_compress),
            ("Decompress DNA", self.show_decompress),
            ("Complementary DNA", self.show_complementary),
            ("Species Matcher", self.show_species_matcher),
            ("Mutation Detector", self.show_mutation_selector)  # Changed to show selector first
        ]
        
        self.menu_buttons = []
        for text, command in menu_options:
            btn = ttk.Button(self.menu_frame, text=text, command=command, width=25)
            btn.pack(pady=5, padx=5, fill=tk.X)
            self.menu_buttons.append(btn)
    
    def create_content_frame(self):
        """Create the main content area"""
        self.content_frame = ttk.Frame(self.main_container, style="TFrame")
        self.content_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
    
    def clear_content(self):
        """Clear the content frame for new content"""
        for widget in self.content_frame.winfo_children():
            widget.destroy()
    
    def show_home(self):
        """Display the home/welcome screen"""
        self.clear_content()
        
        # Welcome message
        ttk.Label(self.content_frame, text="Welcome to DNA Sequencing Toolkit", 
                 font=("Arial", 18, "bold"), style="Header.TLabel").pack(pady=(20, 10))
        
        ttk.Label(self.content_frame, text="Select an option from the menu to get started.", 
                 font=("Arial", 12), style="TLabel").pack(pady=(0, 20))
        
        # Feature frames
        features = [
            ("Compress DNA", "Compress DNA sequences to save storage space"),
            ("Decompress DNA", "Decompress previously compressed DNA sequences"),
            ("Complementary DNA", "Generate the complementary sequence of a DNA strand"),
            ("Species Matcher", "Compare DNA sequences between species"),
            ("Mutation Detector", "Detect mutations in DNA sequences (TP53 or HBB)")
        ]
        
        for i, (title, desc) in enumerate(features):
            frame = ttk.Frame(self.content_frame, style="TFrame")
            frame.pack(fill=tk.X, pady=10, padx=20)
            
            ttk.Label(frame, text=f"{i+1}. {title}", font=("Arial", 12, "bold"), 
                     style="TLabel").pack(anchor=tk.W)
            ttk.Label(frame, text=desc, style="TLabel").pack(anchor=tk.W, padx=(20, 0))
    
    # New method to show mutation type selection
    def show_mutation_selector(self):
        """Show mutation type selection screen"""
        self.clear_content()
        
        ttk.Label(self.content_frame, text="Select Mutation Detection Type", 
                 font=("Arial", 16, "bold"), style="Header.TLabel").pack(pady=(20, 30))
        
        # Description
        ttk.Label(self.content_frame, text="Choose which type of mutation you want to detect:", 
                 font=("Arial", 12), style="TLabel").pack(pady=(0, 20))
        
        # Create a frame for the selection buttons
        btn_frame = ttk.Frame(self.content_frame)
        btn_frame.pack(fill=tk.X, padx=50, pady=10)
        
        # TP53 Button with description
        tp53_frame = ttk.Frame(btn_frame)
        tp53_frame.pack(fill=tk.X, pady=10)
        
        tp53_btn = ttk.Button(tp53_frame, text="TP53 Mutation Detection", 
                             command=lambda: self.show_mutation_detector("tp53"),
                             width=30)
        tp53_btn.pack(pady=5)
        
        ttk.Label(tp53_frame, text="Detect mutations in the TP53 gene(responsible for cell division and cell repairing)", 
                 wraplength=500).pack(pady=(0, 10))
        
        # HBB Button with description
        hbb_frame = ttk.Frame(btn_frame)
        hbb_frame.pack(fill=tk.X, pady=10)
        
        hbb_btn = ttk.Button(hbb_frame, text="HBB Mutation Detection", 
                            command=lambda: self.show_mutation_detector("hbb"),
                            width=30)
        hbb_btn.pack(pady=5)
        
        ttk.Label(hbb_frame, text="Detect mutations in the HBB gene(responsible for hemoglobin production)", 
                 wraplength=500).pack(pady=(0, 10))
        
        # Additional info
        info_frame = ttk.LabelFrame(self.content_frame, text="Information")
        info_frame.pack(fill=tk.X, padx=20, pady=20)
        
        ttk.Label(info_frame, text="TP53 Mutations: Changes in the TP53 gene can lead to loss of tumor suppression, "
                                 "which is associated with many cancer types.\n\n"
                                 "HBB Mutations: Mutations in the HBB gene can lead to sickle cell anemia and other "
                                 "hemoglobin disorders.", 
                 wraplength=600, justify=tk.LEFT).pack(padx=10, pady=10)
    
    def show_mutation_detector(self, mutation_type):
        """Show mutation detector interface for the specified mutation type"""
        self.clear_content()
        self.current_mutation_type = mutation_type
        
        if mutation_type == "tp53":
            title = "TP53 Mutation Detector(responsible for cell division and cell repairing)"
        else:  # hbb
            title = "HBB Mutation Detector(responsible for hemoglobin production)"
        
        ttk.Label(self.content_frame, text=title, 
                 font=("Arial", 16, "bold"), style="Header.TLabel").pack(pady=(20, 30))
        
        # Input frame
        input_frame = ttk.Frame(self.content_frame)
        input_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Label(input_frame, text="Enter DNA sequence to check for mutations:").pack(anchor=tk.W, pady=(0, 5))
        
        self.mutation_dna_text = scrolledtext.ScrolledText(input_frame, height=6, width=50, wrap=tk.WORD)
        self.mutation_dna_text.pack(fill=tk.X)
        
        # Back button and Detect button in the same row
        btn_frame = ttk.Frame(self.content_frame)
        btn_frame.pack(fill=tk.X, padx=20, pady=20)
        
        back_btn = ttk.Button(btn_frame, text="Back to Selection", 
                             command=self.show_mutation_selector)
        back_btn.pack(side=tk.LEFT, padx=(0, 10))
        
        detect_btn = ttk.Button(btn_frame, text="Detect Mutations", 
                              command=self.detect_mutations)
        detect_btn.pack(side=tk.LEFT)
        
        # Results frame
        result_frame = ttk.LabelFrame(self.content_frame, text="Mutation Analysis")
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        self.mutation_result = scrolledtext.ScrolledText(result_frame, height=10, wrap=tk.WORD)
        self.mutation_result.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Status message
        self.mutation_status_var = tk.StringVar()
        ttk.Label(self.content_frame, textvariable=self.mutation_status_var, 
                 font=("Arial", 11, "italic")).pack(pady=10)
    
    def detect_mutations(self):
        try:
            # Get DNA input
            input_dna = self.mutation_dna_text.get("1.0", tk.END).strip()
            if not input_dna:
                messagebox.showerror("Error", "Please enter a DNA sequence.")
                return
            
            # Validate DNA characters
            input_dna = ''.join(input_dna.split()).upper()
            for char in input_dna:
                if char not in ['A', 'T', 'C', 'G']:
                    messagebox.showerror("Error", f"Invalid DNA character found: '{char}'")
                    return
            
            # Get mutation type and run appropriate detector
            if self.current_mutation_type == "tp53":
                # Use existing TP53 detector
                results = get_result2(input_dna)
            else:
                # This would call your HBB detector once implemented
                # results = detect_hbb_mutations(input_dna)
                # For now, show a placeholder message
                results = get_result(input_dna)
            
            # Display results
            #self.mutation_result.delete("1.0", tk.END)
            
            if results:
                if isinstance(results, list):
                    for result in results:
                        self.mutation_result.insert(tk.END, f"{result}\n\n")
                else:
                    self.mutation_result.insert(tk.END, results)
                
                if "✅ Congratulations!" in str(results):
                    self.mutation_status_var.set("No mutations detected!")
                else:
                    self.mutation_status_var.set("Mutations detected or information provided!")
            else:
                self.mutation_result.insert(tk.END, "No results or error in analysis.")
                self.mutation_status_var.set("Analysis complete, but no clear results.")
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.mutation_status_var.set(f"Error: {str(e)}")
    
    # Rest of the class methods (show_compress, show_decompress, etc.) remain unchanged
    def show_compress(self):
        """Show DNA compression interface"""
        self.clear_content()
        
        ttk.Label(self.content_frame, text="Compress DNA Sequence", 
                 font=("Arial", 16, "bold"), style="Header.TLabel").pack(pady=(20, 30))
        
        # Input frame
        input_frame = ttk.Frame(self.content_frame)
        input_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Radio buttons for input type
        self.input_method = tk.StringVar(value="text")
        
        rb_frame = ttk.Frame(input_frame)
        rb_frame.pack(fill=tk.X, pady=10)
        
        ttk.Radiobutton(rb_frame, text="Enter DNA sequence", value="text", 
                       variable=self.input_method, command=self.toggle_compress_input).pack(side=tk.LEFT, padx=(0, 20))
        ttk.Radiobutton(rb_frame, text="Load from file", value="file", 
                       variable=self.input_method, command=self.toggle_compress_input).pack(side=tk.LEFT)
        
        # Text input
        self.dna_text_frame = ttk.Frame(input_frame)
        self.dna_text_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(self.dna_text_frame, text="Enter DNA sequence:").pack(anchor=tk.W, pady=(0, 5))
        self.dna_text = scrolledtext.ScrolledText(self.dna_text_frame, height=10, width=50, wrap=tk.WORD)
        self.dna_text.pack(fill=tk.X)
        
        # File input
        self.file_frame = ttk.Frame(input_frame)
        
        ttk.Label(self.file_frame, text="Select DNA file:").pack(anchor=tk.W, pady=(0, 5))
        
        file_select_frame = ttk.Frame(self.file_frame)
        file_select_frame.pack(fill=tk.X)
        
        self.file_path = tk.StringVar()
        ttk.Entry(file_select_frame, textvariable=self.file_path, width=50).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(file_select_frame, text="Browse", command=self.browse_dna_file).pack(side=tk.RIGHT, padx=(10, 0))
        
        # Output options
        output_frame = ttk.Frame(self.content_frame)
        output_frame.pack(fill=tk.X, padx=20, pady=20)
        
        ttk.Label(output_frame, text="Output file name:").pack(anchor=tk.W, pady=(0, 5))
        
        self.output_name = tk.StringVar(value="compressed_dna.bin")
        ttk.Entry(output_frame, textvariable=self.output_name, width=50).pack(anchor=tk.W)
        
        # Compress button
        ttk.Button(self.content_frame, text="Compress DNA", 
                  command=self.compress_dna).pack(pady=20)
        
        # Status message
        self.status_var = tk.StringVar()
        ttk.Label(self.content_frame, textvariable=self.status_var, 
                 font=("Arial", 11, "italic")).pack(pady=10)
    
    def toggle_compress_input(self):
        if self.input_method.get() == "text":
            self.dna_text_frame.pack(fill=tk.X, pady=10)
            self.file_frame.pack_forget()
        else:
            self.dna_text_frame.pack_forget()
            self.file_frame.pack(fill=tk.X, pady=10)
    
    def browse_dna_file(self):
        file_path = filedialog.askopenfilename(
            title="Select DNA file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if file_path:
            self.file_path.set(file_path)
    
    def compress_dna(self):
        try:
            # Get DNA input
            if self.input_method.get() == "text":
                input_dna = self.dna_text.get("1.0", tk.END).strip()
                if not input_dna:
                    messagebox.showerror("Error", "Please enter a DNA sequence.")
                    return
                
                # Validate DNA characters
                input_dna = ''.join(input_dna.split()).upper()
                for char in input_dna:
                    if char not in ['A', 'T', 'C', 'G']:
                        messagebox.showerror("Error", f"Invalid DNA character found: '{char}'")
                        return
            else:
                # Read from file
                path = self.file_path.get()
                if not path:
                    messagebox.showerror("Error", "Please select a file.")
                    return
                try:
                    input_dna = read_dna_file(path)
                    # Fix the tuple issue from your original read_dna_file function
                    if isinstance(input_dna, tuple):
                        input_dna = input_dna[0]
                except Exception as e:
                    messagebox.showerror("Error", str(e))
                    return
            
            # Compress DNA
            compressed_dna = compress(input_dna)
            
            # Save to file
            output_path = os.path.join("compressed_files", self.output_name.get())
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            with open(output_path, "wb") as f:
                f.write(len(input_dna).to_bytes(4, byteorder='big'))
                f.write(compressed_dna)
            
            self.status_var.set(f"Compressed DNA saved to {output_path}")
            messagebox.showinfo("Success", f"DNA compressed successfully!\nSaved to {output_path}")
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_var.set(f"Error: {str(e)}")
    
    def show_decompress(self):
        """Show DNA decompression interface"""
        self.clear_content()
        
        ttk.Label(self.content_frame, text="Decompress DNA Sequence", 
                 font=("Arial", 16, "bold"), style="Header.TLabel").pack(pady=(20, 30))
        
        # File input
        file_frame = ttk.Frame(self.content_frame)
        file_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Label(file_frame, text="Select compressed DNA file (.bin):").pack(anchor=tk.W, pady=(0, 5))
        
        file_select_frame = ttk.Frame(file_frame)
        file_select_frame.pack(fill=tk.X)
        
        self.compressed_file_path = tk.StringVar()
        ttk.Entry(file_select_frame, textvariable=self.compressed_file_path, width=50).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(file_select_frame, text="Browse", command=self.browse_compressed_file).pack(side=tk.RIGHT, padx=(10, 0))
        
        # Output options
        output_frame = ttk.Frame(self.content_frame)
        output_frame.pack(fill=tk.X, padx=20, pady=20)
        
        ttk.Label(output_frame, text="Output file name:").pack(anchor=tk.W, pady=(0, 5))
        
        self.decomp_output_name = tk.StringVar(value="decompressed_dna.txt")
        ttk.Entry(output_frame, textvariable=self.decomp_output_name, width=50).pack(anchor=tk.W)
        
        # Decompress button
        ttk.Button(self.content_frame, text="Decompress DNA", 
                  command=self.decompress_dna).pack(pady=20)
        
        # Results frame
        result_frame = ttk.LabelFrame(self.content_frame, text="Decompressed DNA")
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        self.decomp_result = scrolledtext.ScrolledText(result_frame, height=10, wrap=tk.WORD)
        self.decomp_result.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Status message
        self.decomp_status_var = tk.StringVar()
        ttk.Label(self.content_frame, textvariable=self.decomp_status_var, 
                 font=("Arial", 11, "italic")).pack(pady=10)
    
    def browse_compressed_file(self):
        file_path = filedialog.askopenfilename(
            title="Select compressed DNA file",
            filetypes=[("Binary files", "*.bin"), ("All files", "*.*")]
        )
        if file_path:
            self.compressed_file_path.set(file_path)
    
    def decompress_dna(self):
        try:
            # Get compressed file path
            bin_path = self.compressed_file_path.get()
            if not bin_path:
                messagebox.showerror("Error", "Please select a compressed file.")
                return
            
            # Decompress DNA
            with open(bin_path, "rb") as f:
                original_length = int.from_bytes(f.read(4), byteorder='big')
                compressed_data = bytearray(f.read())
                decompressed_dna = decompress(compressed_data, original_length)
            
            # Save to file
            output_path = self.decomp_output_name.get()
            with open(output_path, "w") as f:
                f.write(decompressed_dna)
            
            # Display result
            self.decomp_result.delete("1.0", tk.END)
            if len(decompressed_dna) > 1000:
                preview = decompressed_dna[:1000] + "... (truncated, see saved file for complete sequence)"
                self.decomp_result.insert(tk.END, preview)
            else:
                self.decomp_result.insert(tk.END, decompressed_dna)
            
            self.decomp_status_var.set(f"Decompressed DNA saved to {output_path}")
            messagebox.showinfo("Success", f"DNA decompressed successfully!\nSaved to {output_path}")
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.decomp_status_var.set(f"Error: {str(e)}")
    
    def show_complementary(self):
        """Show complementary DNA interface"""
        self.clear_content()
        
        ttk.Label(self.content_frame, text="Generate Complementary DNA", 
                 font=("Arial", 16, "bold"), style="Header.TLabel").pack(pady=(20, 30))
        
        # Input frame
        input_frame = ttk.Frame(self.content_frame)
        input_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Label(input_frame, text="Enter DNA sequence:").pack(anchor=tk.W, pady=(0, 5))
        
        self.comp_dna_text = scrolledtext.ScrolledText(input_frame, height=6, width=50, wrap=tk.WORD)
        self.comp_dna_text.pack(fill=tk.X)
        
        # Generate button
        ttk.Button(self.content_frame, text="Generate Complementary Sequence", 
                  command=self.generate_complementary).pack(pady=20)
        
        # Results frame
        result_frame = ttk.LabelFrame(self.content_frame, text="Complementary DNA")
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        self.comp_result = scrolledtext.ScrolledText(result_frame, height=6, wrap=tk.WORD)
        self.comp_result.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Status message
        self.comp_status_var = tk.StringVar()
        ttk.Label(self.content_frame, textvariable=self.comp_status_var, 
                 font=("Arial", 11, "italic")).pack(pady=10)
    
    def generate_complementary(self):
        try:
            # Get DNA input
            input_dna = self.comp_dna_text.get("1.0", tk.END).strip()
            if not input_dna:
                messagebox.showerror("Error", "Please enter a DNA sequence.")
                return
            
            # Validate DNA characters
            input_dna = ''.join(input_dna.split()).upper()
            for char in input_dna:
                if char not in ['A', 'T', 'C', 'G']:
                    messagebox.showerror("Error", f"Invalid DNA character found: '{char}'")
                    return
            
            # Generate complementary sequence
            complementary_dna = get_complementary_sequence(input_dna)
            
            # Display result
            self.comp_result.delete("1.0", tk.END)
            self.comp_result.insert(tk.END, complementary_dna)
            
            self.comp_status_var.set("Complementary DNA generated successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.comp_status_var.set(f"Error: {str(e)}")
    
    def show_species_matcher(self):
        """Show species matcher interface"""
        self.clear_content()
        
        ttk.Label(self.content_frame, text="DNA Species Matcher", 
                 font=("Arial", 16, "bold"), style="Header.TLabel").pack(pady=(20, 30))
        
        # Species name input
        name_frame = ttk.Frame(self.content_frame)
        name_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Label(name_frame, text="Enter species name:").pack(anchor=tk.W, pady=(0, 5))
        
        self.species_name_var = tk.StringVar()
        ttk.Entry(name_frame, textvariable=self.species_name_var, width=40).pack(anchor=tk.W)
        
        # DNA sequence input
        seq_frame = ttk.Frame(self.content_frame)
        seq_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Label(seq_frame, text="Enter DNA sequence for the species:").pack(anchor=tk.W, pady=(0, 5))
        
        self.species_dna_text = scrolledtext.ScrolledText(seq_frame, height=6, width=50, wrap=tk.WORD)
        self.species_dna_text.pack(fill=tk.X)
        
        # Compare button
        ttk.Button(self.content_frame, text="Compare with Human DNA", 
                  command=self.compare_species).pack(pady=20)
        
        # Results frame
        result_frame = ttk.LabelFrame(self.content_frame, text="Comparison Results")
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        self.species_result = scrolledtext.ScrolledText(result_frame, height=6, wrap=tk.WORD)
        self.species_result.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Graph frame
        self.graph_frame = ttk.Frame(self.content_frame)
        self.graph_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # Status message
        self.species_status_var = tk.StringVar()
        ttk.Label(self.content_frame, textvariable=self.species_status_var, 
                 font=("Arial", 11, "italic")).pack(pady=10)
    
    def compare_species(self):
        try:
            # Get inputs
            species_name = self.species_name_var.get().strip()
            if not species_name:
                messagebox.showerror("Error", "Please enter a species name.")
                return
            
            species_seq = self.species_dna_text.get("1.0", tk.END).strip()
            if not species_seq:
                messagebox.showerror("Error", "Please enter a DNA sequence.")
                return
            
            # Clean sequence
            species_seq = ''.join(s.upper() for s in species_seq if s.upper() in 'ATGC')
            
            if not species_seq:
                messagebox.showerror("Error", "The sequence provided is empty or contains no valid DNA bases (A, T, G, C).")
                return
            
            # Create thread for analysis to prevent UI freeze
            self.species_status_var.set("Analyzing DNA sequences... Please wait.")
            self.update()
            
            threading.Thread(target=self.run_species_analysis, 
                           args=(species_name, species_seq)).start()
            
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.species_status_var.set(f"Error: {str(e)}")
    
    def run_species_analysis(self, species_name, species_seq):
        try:
            # Human COI reference sequence (from your species_matcher.py)
            human_coi = """
            GTCCTACTATCCATGCAGGTATCTTCTATCTTTGGGGCATGAGCGGGCATAGTAGGCACAGCCCTAAGCCTCCTCATTCG
            AGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTCGTAA
            TAATCTTCTTCATAGTAATACCAATAATAATCGGAGGCTTCGGAAACTGACTAGTCCCCCTTATAATCGGTGCCCCCGAC
            ATAGCATTCCCACGAATAAATAACATAAGCTTCTGACTCCTCCCTCCATCCTTTCTCCTTCTTCTCGCATCCTCCGGAGT
            AGAAGCTGGCGCAGGAACAGGCTGAACAGTCTACCCTCCCCTAGCAGGAAACTACTCCCACCCTGGAGCCTCCGTAGACC
            TGGCAATCTTCTCCCTCCACCTGGCAGGTATTTCCTCCATCCTAGGAGCAATTAACTTCATCACCACAGCTATCAACATA
            AAACCCCCTGCAATATCCCAGTATCAAACTCCCCTATTCGTCTGATCAGTCCTAATTACCGCCGTCCTACTCCTCCTGTC
            CCTGCCCGTCCTCGCTGCAGGAATCACAATGCTGCTCACAGACCGCAACCTTAACACCACCTTCTTCGACCCGGCAGGAG
            GAGGAGACCCAGTCCTGTACCAACACCTATTCT
            """
            human_coi = ''.join(s.upper() for s in human_coi if s.upper() in 'ATGC')
            
            # Calculate overall similarity
            from analysis.species_matcher import calculate_similarity, find_diagnostic_sites, sliding_window_similarity
            
            similarity = calculate_similarity(human_coi, species_seq)
            
            # Get sliding window similarities for graph
            window_similarities = sliding_window_similarity(human_coi, species_seq)
            
            # Find diagnostic sites
            min_len = min(len(human_coi), len(species_seq))
            differences = []
            for pos in range(min_len):
                if human_coi[pos] != species_seq[pos]:
                    differences.append(f"Position {pos+1}: Human={human_coi[pos]}, {species_name}={species_seq[pos]}")
                    if len(differences) >= 10:
                        differences.append(f"... and more differences")
                        break
            
            # Update the UI in the main thread
            self.after(0, lambda: self.update_species_results(
                species_name, similarity, differences, window_similarities))
            
        except Exception as e:
            self.after(0, lambda: messagebox.showerror("Error"))
            self.after(0, lambda: self.species_status_var.set("Error:kuch to gadbad h"))
    
    def update_species_results(self, species_name, similarity, differences, window_similarities):
        # Update text results
        self.species_result.delete("1.0", tk.END)
        self.species_result.insert(tk.END, f"Overall DNA Similarity: {similarity:.2f}%\n\n")
        self.species_result.insert(tk.END, f"Diagnostic Sites where {species_name} differs from human:\n")
        self.species_result.insert(tk.END, "="*50 + "\n")
        
        for diff in differences:
            self.species_result.insert(tk.END, f"  {diff}\n")
        
        # Clear previous graph if any
        for widget in self.graph_frame.winfo_children():
            widget.destroy()
        
        # Create new graph
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(window_similarities)
        ax.set_title(f"DNA Similarity: Human vs {species_name}")
        ax.set_xlabel("Position (50 bp window)")
        ax.set_ylabel("Percent Similarity")
        ax.set_ylim(0, 100)
        ax.axhline(y=90, color='g', linestyle='--', alpha=0.3, label="High conservation (90%)")
        ax.axhline(y=50, color='r', linestyle='--', alpha=0.3, label="Low conservation (50%)")
        ax.legend()
        
        canvas = FigureCanvasTkAgg(fig, master=self.graph_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Update status
        self.species_status_var.set("Analysis complete!")
        
        # Save the plot to file
        try:
            plt.savefig("dna_comparison.png")
        except:
            pass

if __name__ == "__main__":
    app = DNASequencingApp()
    app.mainloop()
