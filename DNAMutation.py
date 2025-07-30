import streamlit as st
import random
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio.Seq import Seq
from collections import Counter
import base64

# Page configuration
st.set_page_config(
    page_title="GenoMutant Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for styling
st.markdown("""
<style>
    .main-header {
        text-align: center;
        padding: 2rem 0;
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        color: white;
        border-radius: 10px;
        margin-bottom: 2rem;
    }
    .subtitle {
        text-align: center;
        color: #666;
        font-size: 1.2em;
        margin-bottom: 2rem;
    }
    .mutation-highlight {
        background-color: #ffebee;
        padding: 1rem;
        border-radius: 5px;
        border-left: 4px solid #f44336;
    }
    .success-box {
        background-color: #e8f5e8;
        padding: 1rem;
        border-radius: 5px;
        border-left: 4px solid #4caf50;
    }
    .footer {
        text-align: center;
        padding: 2rem 0;
        color: #888;
        border-top: 1px solid #eee;
        margin-top: 3rem;
    }
    .metric-card {
        background: white;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        text-align: center;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("""
<div class="main-header">
    <h1>🧬 GenoMutant Pro</h1>
    <h3>Advanced DNA Mutation Analysis & Visualization Platform</h3>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="subtitle">
    Explore the fascinating world of genetic mutations and their impact on protein sequences. 
    Designed for researchers, educators, and students in molecular biology and genetics.
</div>
""", unsafe_allow_html=True)

# Helper Functions
@st.cache_data
def validate_dna(sequence):
    """Validate DNA sequence contains only valid nucleotides."""
    return all(base.upper() in 'ATGC' for base in sequence.replace(' ', ''))

@st.cache_data
def get_codon_and_aa(sequence, nucleotide_pos):
    """Extract codon and amino acid at specific position."""
    if not (0 <= nucleotide_pos < len(sequence)):
        return None, None
    
    codon_start_pos = (nucleotide_pos // 3) * 3
    codon_end_pos = codon_start_pos + 3
    codon = sequence[codon_start_pos:codon_end_pos]
    
    if len(codon) == 3:
        amino_acid = str(Seq(codon).translate(table=1))
        return codon, amino_acid
    return codon, None

@st.cache_data
def generate_random_dna(length):
    """Generate random DNA sequence."""
    return ''.join(random.choice('ATGC') for _ in range(length))

def apply_point_mutation(dna_seq, pos, new_base):
    """Apply point mutation."""
    if not (0 <= pos < len(dna_seq)):
        return None, None
    original_base = dna_seq[pos]
    mutated_seq = dna_seq[:pos] + new_base + dna_seq[pos+1:]
    return mutated_seq, original_base

def apply_insertion_mutation(dna_seq, pos, inserted_seq):
    """Apply insertion mutation."""
    if not (0 <= pos <= len(dna_seq)) or not inserted_seq:
        return None, None
    mutated_seq = dna_seq[:pos] + inserted_seq + dna_seq[pos:]
    return mutated_seq, inserted_seq

def apply_deletion_mutation(dna_seq, pos, delete_length):
    """Apply deletion mutation."""
    if delete_length <= 0 or not (0 <= pos < len(dna_seq)) or (pos + delete_length) > len(dna_seq):
        return None, None
    deleted_seq = dna_seq[pos:pos + delete_length]
    mutated_seq = dna_seq[:pos] + dna_seq[pos + delete_length:]
    return mutated_seq, deleted_seq

def create_sequence_comparison_chart(original_seq, mutated_seq, mutation_pos, mutation_type):
    """Create visualization comparing original and mutated sequences."""
    fig = go.Figure()
    
    # Base composition analysis
    orig_composition = Counter(original_seq)
    mut_composition = Counter(mutated_seq)
    
    bases = ['A', 'T', 'G', 'C']
    orig_counts = [orig_composition.get(base, 0) for base in bases]
    mut_counts = [mut_composition.get(base, 0) for base in bases]
    
    fig.add_trace(go.Bar(
        name='Original',
        x=bases,
        y=orig_counts,
        marker_color='lightblue'
    ))
    
    fig.add_trace(go.Bar(
        name='Mutated',
        x=bases,
        y=mut_counts,
        marker_color='lightcoral'
    ))
    
    fig.update_layout(
        title=f'Base Composition Comparison ({mutation_type.title()} Mutation)',
        xaxis_title='Nucleotide Base',
        yaxis_title='Count',
        barmode='group',
        height=400
    )
    
    return fig

def analyze_mutation_impact(original_dna, mutated_dna, mutation_type, mutation_pos, change_length=1):
    """Comprehensive mutation impact analysis."""
    try:
        original_protein = str(Seq(original_dna).translate(table=1, to_stop=True))
        mutated_protein = str(Seq(mutated_dna).translate(table=1, to_stop=True))
    except:
        original_protein = "Translation error"
        mutated_protein = "Translation error"
    
    analysis = {
        'mutation_type': mutation_type,
        'position': mutation_pos,
        'original_dna': original_dna,
        'mutated_dna': mutated_dna,
        'original_protein': original_protein,
        'mutated_protein': mutated_protein,
        'dna_length_change': len(mutated_dna) - len(original_dna),
        'protein_length_change': len(mutated_protein) - len(original_protein)
    }
    
    # Determine mutation impact
    if mutation_type == 'point':
        orig_codon, orig_aa = get_codon_and_aa(original_dna, mutation_pos)
        mut_codon, mut_aa = get_codon_and_aa(mutated_dna, mutation_pos)
        
        analysis['original_codon'] = orig_codon
        analysis['mutated_codon'] = mut_codon
        analysis['original_aa'] = orig_aa
        analysis['mutated_aa'] = mut_aa
        
        if orig_aa == mut_aa:
            analysis['impact'] = 'Silent Mutation'
            analysis['description'] = 'No change in amino acid sequence'
            analysis['severity'] = 'Low'
        elif mut_aa == '*':
            analysis['impact'] = 'Nonsense Mutation'
            analysis['description'] = 'Premature stop codon introduced'
            analysis['severity'] = 'High'
        else:
            analysis['impact'] = 'Missense Mutation'
            analysis['description'] = f'Amino acid changed from {orig_aa} to {mut_aa}'
            analysis['severity'] = 'Medium'
    
    elif mutation_type in ['insertion', 'deletion']:
        is_frameshift = (change_length % 3 != 0)
        
        if is_frameshift:
            analysis['impact'] = 'Frameshift Mutation'
            analysis['description'] = 'Reading frame altered, affecting all downstream amino acids'
            analysis['severity'] = 'High'
        else:
            analysis['impact'] = f'In-frame {mutation_type.title()}'
            analysis['description'] = f'Amino acids {"added" if mutation_type == "insertion" else "removed"} without frameshift'
            analysis['severity'] = 'Medium'
    
    return analysis

# Sidebar
st.sidebar.header("🔬 Mutation Parameters")

# DNA Input Section
input_method = st.sidebar.selectbox(
    "Choose DNA input method:",
    ["Enter custom sequence", "Generate random sequence", "Use example sequence"]
)

if input_method == "Enter custom sequence":
    dna_input = st.sidebar.text_area(
        "Enter DNA sequence:",
        placeholder="ATGCGTGACTGACTGACGTA",
        help="Enter only A, T, G, C nucleotides"
    ).upper().replace(' ', '')
elif input_method == "Generate random sequence":
    seq_length = st.sidebar.slider("Sequence length:", 20, 500, 100)
    if st.sidebar.button("Generate Random DNA"):
        st.session_state.random_dna = generate_random_dna(seq_length)
    dna_input = st.session_state.get('random_dna', generate_random_dna(seq_length))
else:
    examples = {
        "Hemoglobin Beta": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGGAGAAGTAACAGGGTACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAGATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATTTCTGCATATAAATATTTCTGCATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGA",
        "Insulin": "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCGCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC",
        "p53 Tumor Suppressor": "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCG"
    }
    selected_example = st.sidebar.selectbox("Choose example:", list(examples.keys()))
    dna_input = examples[selected_example]

# Validation
if dna_input and not validate_dna(dna_input):
    st.sidebar.error("⚠️ Invalid DNA sequence! Please use only A, T, G, C nucleotides.")
    st.stop()

# Main content
if dna_input:
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("📊 Original Sequence Analysis")
        
        # Sequence metrics
        metrics_col1, metrics_col2, metrics_col3, metrics_col4 = st.columns(4)
        
        with metrics_col1:
            st.metric("Length", f"{len(dna_input)} bp")
        
        with metrics_col2:
            gc_content = (dna_input.count('G') + dna_input.count('C')) / len(dna_input) * 100
            st.metric("GC Content", f"{gc_content:.1f}%")
        
        with metrics_col3:
            try:
                protein_seq = str(Seq(dna_input).translate(table=1, to_stop=True))
                st.metric("Protein Length", f"{len(protein_seq)} aa")
            except:
                st.metric("Protein Length", "Error")
        
        with metrics_col4:
            st.metric("Codons", f"{len(dna_input)//3}")
        
        # Display sequences in expandable sections
        with st.expander("🧬 View DNA Sequence", expanded=False):
            # Format DNA sequence in blocks of 10 for readability
            formatted_dna = ' '.join([dna_input[i:i+10] for i in range(0, len(dna_input), 10)])
            st.code(formatted_dna, language=None)
        
        try:
            protein_seq = str(Seq(dna_input).translate(table=1, to_stop=True))
            with st.expander("🧪 View Protein Sequence", expanded=False):
                formatted_protein = ' '.join([protein_seq[i:i+10] for i in range(0, len(protein_seq), 10)])
                st.code(formatted_protein, language=None)
        except:
            st.warning("Unable to translate protein sequence")
    
    with col2:
        st.subheader("⚙️ Mutation Settings")
        
        mutation_type = st.selectbox(
            "Mutation type:",
            ["point", "insertion", "deletion"],
            help="Choose the type of mutation to apply"
        )
        
        if mutation_type == "point":
            position = st.number_input(
                "Position (0-indexed):",
                min_value=0,
                max_value=len(dna_input)-1,
                value=len(dna_input)//2
            )
            
            current_base = dna_input[position]
            st.info(f"Current base at position {position}: **{current_base}**")
            
            available_bases = [b for b in ['A', 'T', 'G', 'C'] if b != current_base]
            new_base = st.selectbox("New base:", available_bases)
            
        elif mutation_type == "insertion":
            position = st.number_input(
                "Insertion position:",
                min_value=0,
                max_value=len(dna_input),
                value=len(dna_input)//2
            )
            
            inserted_sequence = st.text_input(
                "Sequence to insert:",
                value="ATG",
                help="Enter the DNA sequence to insert"
            ).upper()
            
        else:  # deletion
            position = st.number_input(
                "Deletion start position:",
                min_value=0,
                max_value=len(dna_input)-1,
                value=len(dna_input)//2
            )
            
            max_delete = len(dna_input) - position
            delete_length = st.number_input(
                "Number of bases to delete:",
                min_value=1,
                max_value=max_delete,
                value=min(3, max_delete)
            )

# Apply Mutation Button
if st.button("🧬 Apply Mutation", type="primary", use_container_width=True):
    if not dna_input:
        st.error("Please provide a DNA sequence first!")
    else:
        # Apply mutation based on type
        if mutation_type == "point":
            result = apply_point_mutation(dna_input, position, new_base)
            if result[0]:
                mutated_dna = result[0]
                change_length = 1
                st.session_state.mutation_applied = True
                st.session_state.mutation_data = {
                    'original': dna_input,
                    'mutated': mutated_dna,
                    'type': mutation_type,
                    'position': position,
                    'change_length': change_length
                }
            else:
                st.error("Failed to apply point mutation!")
                
        elif mutation_type == "insertion":
            if validate_dna(inserted_sequence):
                result = apply_insertion_mutation(dna_input, position, inserted_sequence)
                if result[0]:
                    mutated_dna = result[0]
                    change_length = len(inserted_sequence)
                    st.session_state.mutation_applied = True
                    st.session_state.mutation_data = {
                        'original': dna_input,
                        'mutated': mutated_dna,
                        'type': mutation_type,
                        'position': position,
                        'change_length': change_length
                    }
                else:
                    st.error("Failed to apply insertion!")
            else:
                st.error("Invalid insertion sequence! Use only A, T, G, C.")
                
        else:  # deletion
            result = apply_deletion_mutation(dna_input, position, delete_length)
            if result[0]:
                mutated_dna = result[0]
                change_length = delete_length
                st.session_state.mutation_applied = True
                st.session_state.mutation_data = {
                    'original': dna_input,
                    'mutated': mutated_dna,
                    'type': mutation_type,
                    'position': position,
                    'change_length': change_length
                }
            else:
                st.error("Failed to apply deletion!")

# Display Results
if st.session_state.get('mutation_applied', False):
    mutation_data = st.session_state.mutation_data
    
    st.markdown("---")
    st.header("🔬 Mutation Analysis Results")
    
    # Analyze mutation impact
    analysis = analyze_mutation_impact(
        mutation_data['original'],
        mutation_data['mutated'],
        mutation_data['type'],
        mutation_data['position'],
        mutation_data['change_length']
    )
    
    # Results overview
    col1, col2, col3 = st.columns(3)
    
    with col1:
        severity_color = {"Low": "🟢", "Medium": "🟡", "High": "🔴"}
        st.metric(
            "Mutation Impact",
            analysis['impact'],
            help=f"Severity: {severity_color.get(analysis['severity'], '⚪')} {analysis['severity']}"
        )
    
    with col2:
        st.metric(
            "DNA Length Change",
            f"{analysis['dna_length_change']:+d} bp"
        )
    
    with col3:
        st.metric(
            "Protein Length Change",
            f"{analysis['protein_length_change']:+d} aa"
        )
    
    # Detailed analysis
    st.subheader("📋 Detailed Analysis")
    
    impact_color = {
        "Silent Mutation": "success",
        "Missense Mutation": "warning", 
        "Nonsense Mutation": "error",
        "Frameshift Mutation": "error",
        "In-frame Insertion": "warning",
        "In-frame Deletion": "warning"
    }
    
    alert_type = impact_color.get(analysis['impact'], "info")
    
    if alert_type == "success":
        st.success(f"**{analysis['impact']}**: {analysis['description']}")
    elif alert_type == "warning":
        st.warning(f"**{analysis['impact']}**: {analysis['description']}")
    elif alert_type == "error":
        st.error(f"**{analysis['impact']}**: {analysis['description']}")
    else:
        st.info(f"**{analysis['impact']}**: {analysis['description']}")
    
    # Sequence comparison
    st.subheader("🔍 Sequence Comparison")
    
    comp_col1, comp_col2 = st.columns(2)
    
    with comp_col1:
        st.write("**Original DNA:**")
        orig_formatted = ' '.join([mutation_data['original'][i:i+10] for i in range(0, len(mutation_data['original']), 10)])
        st.code(orig_formatted, language=None)
        
        st.write("**Original Protein:**")
        if analysis['original_protein'] != "Translation error":
            orig_prot_formatted = ' '.join([analysis['original_protein'][i:i+10] for i in range(0, len(analysis['original_protein']), 10)])
            st.code(orig_prot_formatted, language=None)
        else:
            st.code("Translation error", language=None)
    
    with comp_col2:
        st.write("**Mutated DNA:**")
        mut_formatted = ' '.join([mutation_data['mutated'][i:i+10] for i in range(0, len(mutation_data['mutated']), 10)])
        st.code(mut_formatted, language=None)
        
        st.write("**Mutated Protein:**")
        if analysis['mutated_protein'] != "Translation error":
            mut_prot_formatted = ' '.join([analysis['mutated_protein'][i:i+10] for i in range(0, len(analysis['mutated_protein']), 10)])
            st.code(mut_prot_formatted, language=None)
        else:
            st.code("Translation error", language=None)
    
    # Codon analysis for point mutations
    if mutation_data['type'] == 'point' and 'original_codon' in analysis:
        st.subheader("🧪 Codon Analysis")
        
        codon_col1, codon_col2 = st.columns(2)
        
        with codon_col1:
            st.write(f"**Original Codon:** `{analysis['original_codon']}` → `{analysis['original_aa']}`")
        
        with codon_col2:
            st.write(f"**Mutated Codon:** `{analysis['mutated_codon']}` → `{analysis['mutated_aa']}`")
    
    # Visualization
    st.subheader("📊 Base Composition Analysis")
    fig = create_sequence_comparison_chart(
        mutation_data['original'],
        mutation_data['mutated'],
        mutation_data['position'],
        mutation_data['type']
    )
    st.plotly_chart(fig, use_container_width=True)
    
    # Download results
    st.subheader("💾 Export Results")
    
    results_data = {
        'Parameter': ['Original DNA', 'Mutated DNA', 'Original Protein', 'Mutated Protein', 
                     'Mutation Type', 'Position', 'Impact', 'Severity', 'Description'],
        'Value': [mutation_data['original'], mutation_data['mutated'], 
                 analysis['original_protein'], analysis['mutated_protein'],
                 analysis['mutation_type'], analysis['position'], 
                 analysis['impact'], analysis['severity'], analysis['description']]
    }
    
    results_df = pd.DataFrame(results_data)
    csv = results_df.to_csv(index=False)
    
    st.download_button(
        label="📥 Download Results as CSV",
        data=csv,
        file_name=f"mutation_analysis_{analysis['mutation_type']}_{analysis['position']}.csv",
        mime="text/csv"
    )

# Educational content
with st.expander("📚 Learn About DNA Mutations", expanded=False):
    st.markdown("""
    ### Types of DNA Mutations
    
    **Point Mutations (Substitutions):**
    - **Silent**: No change in amino acid sequence
    - **Missense**: Different amino acid produced
    - **Nonsense**: Premature stop codon created
    
    **Insertion/Deletion Mutations:**
    - **Frameshift**: Not divisible by 3, shifts reading frame
    - **In-frame**: Divisible by 3, maintains reading frame
    
    ### Genetic Code Properties
    - **Redundancy**: Multiple codons can code for the same amino acid
    - **Reading Frame**: DNA is read in triplets (codons)
    - **Start/Stop Codons**: ATG starts translation, TAA/TAG/TGA stop it
    
    ### Clinical Relevance
    Different mutations have varying impacts on protein function and organism health.
    This tool helps understand these relationships for educational and research purposes.
    """)

# Footer
st.markdown("""
<div class="footer">
    <p><strong>GenoMutant Pro</strong> | Advanced DNA Mutation Analysis Platform</p>
    <p>Built with Streamlit & BioPython | For Educational & Research Use</p>
    <p>💡 <em>Advancing understanding of genetic variation and molecular biology</em></p>
</div>
""", unsafe_allow_html=True)
