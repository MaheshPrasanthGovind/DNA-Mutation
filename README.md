#GenoMutant Pro 🧬

Advanced DNA Mutation Analysis & Visualization Platform

GenoMutant Pro is a web-based application built with Streamlit and Biopython that enables users to explore, simulate, and analyze DNA mutations and their impact on protein sequences. Designed for educational, research, and exploratory purposes, it provides both interactive visualization and detailed genetic analysis.

Features
1. DNA Sequence Input

Custom DNA entry: Users can input any valid DNA sequence.

Random sequence generation: Generate synthetic DNA sequences of variable lengths.

Example sequences: Preloaded biologically relevant sequences, e.g., Hemoglobin Beta, Insulin, p53 Tumor Suppressor.

2. Mutation Simulation

Point mutations: Substitute a single nucleotide.

Insertions: Add new nucleotides at a specified position.

Deletions: Remove nucleotides from the sequence.

Automatic detection of frameshift vs. in-frame mutations.

3. Analysis

Protein translation: Converts DNA sequences to amino acid sequences using the standard genetic code.

Mutation impact assessment: Categorizes mutations as silent, missense, nonsense, frameshift, or in-frame.

Calculates changes in DNA and protein length.

Codon-level analysis for point mutations.

4. Visualization

Base composition comparison: Interactive Plotly charts showing nucleotide distributions before and after mutation.

Optional color coding to highlight mutation impacts.

5. Export & Sharing

Download mutation analysis results as a CSV file.

Easily shareable and reproducible for research or educational purposes.

6. Educational Resources

In-app guide on mutation types, genetic code properties, and clinical relevance.

Designed to enhance understanding of molecular biology concepts.

Tech Stack

Frontend: Streamlit

Backend / Bioinformatics: Biopython, Plotly, Pandas

Python Version: 3.9+

Deployment: Can be hosted on Streamlit Cloud, Heroku, or any Python-compatible web server.
