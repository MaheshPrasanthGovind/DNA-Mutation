## 🧬🧬Try out the app using the code : [GenoMutant Pro](https://genomutantpro-mutationanalysis-gssfefbmyurpqnjjzwwf32.streamlit.app/)
# GenoMutant Pro 🧬

Advanced DNA Mutation Analysis & Visualization Platform

GenoMutant Pro is a web-based application built with Streamlit and Biopython that enables users to explore, simulate, and analyze DNA mutations and their impact on protein sequences. Designed for educational, research, and exploratory purposes, it provides both interactive visualization and detailed genetic analysis.

### Features:

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

## Tech Stack

Frontend: Streamlit
Backend / Bioinformatics: Biopython, Plotly, Pandas
Python Version: 3.9+
Deployment: Can be hosted on Streamlit Cloud, Heroku, or any Python-compatible web server.


## Educational Value

GenoMutant Pro is designed to bridge the gap between theory and practice in genetics:
Visualize the effects of mutations on both DNA and protein sequences.
Explore frameshift vs. in-frame mutations interactively.
Understand the molecular basis of diseases related to point mutations or indels.
This makes it an ideal tool for students, educators, and early-stage researchers.

## Future Improvements

Highlight mutated amino acids directly in the protein sequence.
Allow multiple simultaneous mutations.
Integrate 3D protein structure visualization.
Enable batch mutation simulations with statistical analysis.

## What I Learned

Developing GenoMutant Pro was an invaluable learning experience that strengthened both my technical and scientific skills:

Bioinformatics & Genetics: Gained hands-on experience with DNA and protein sequences, codon translation, and mutation impacts, reinforcing my understanding of molecular biology concepts such as frameshifts, missense, and nonsense mutations.
Python & Libraries: Deepened my proficiency in Python, Streamlit, Biopython, Plotly, and Pandas for building interactive, data-driven web applications.
Data Visualization: Learned to design intuitive, informative visualizations to communicate complex biological data effectively.
Software Design & User Experience: Practiced designing a user-friendly interface that balances educational content, analysis tools, and visual appeal.
Problem-Solving & Debugging: Navigated challenges with sequence validation, translation errors, and mutation logic, which taught me systematic debugging and careful attention to edge cases.
Project Management: Improved skills in structuring a project from concept to deployment, including modular coding, caching for efficiency, and documenting workflows for others.
This project helped me combine my passion for biology with computational skills, and I now feel more confident tackling real-world bioinformatics problems in research or academic settings.

## License
This project is for educational and research purposes only. Commercial use is prohibited without prior permission.
