import streamlit as st
import random
import re

# --- Set Page Config and Custom CSS for Dark Theme ---
st.set_page_config(page_title="🧬 DNA Mutation Lab", page_icon="🔬", layout="wide")

st.markdown(
    """
    <style>
    body {
        background-color: #000000;
        color: white;
    }
    .output-box {
        background-color: #444444;
        padding: 10px;
        border-radius: 8px;
        font-family: monospace;
        white-space: pre-wrap;
        word-break: break-word;
        border: 1px solid #6272a4; /* Added a subtle border */
    }
    .highlight-red {
        color: #FF5555; /* A slightly softer red */
        font-weight: bold;
        background-color: #330000; /* Subtle background for highlight */
        padding: 2px 0px; /* Adjust padding for better look */
        border-radius: 3px;
    }
    .stButton>button {
        background-color: #6272a4; /* Darker background for buttons */
        color: white;
        border-radius: 8px;
        border: none;
        padding: 10px 20px;
        font-size: 16px;
        transition: background-color 0.3s ease;
    }
    .stButton>button:hover {
        background-color: #8be9fd; /* Lighter on hover */
        color: black;
    }
    .stSelectbox, .stNumberInput, .stTextInput, .stTextArea {
        background-color: #282a36; /* Darker input fields */
        border-radius: 8px;
        padding: 10px;
    }
    .stTextInput>div>div>input, .stTextArea>div>div>textarea {
        background-color: #282a36;
        color: #f8f8f2;
        border: 1px solid #44475a;
        border-radius: 8px;
        padding: 8px;
    }
    .stSelectbox>div>div>div>span {
        color: #f8f8f2;
    }
    .stNumberInput>div>div>input {
        background-color: #282a36;
        color: #f8f8f2;
        border: 1px solid #44475a;
        border-radius: 8px;
        padding: 8px;
    }
    .biological-impact-box {
        background-color: #282a36;
        padding: 15px;
        border-radius: 8px;
        border: 1px solid #bd93f9; /* Purple border for biological info */
        margin-top: 15px;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# --- Genetic Code Dictionary ---
GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',  # * for Stop
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',  # AUG is Start and Met
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


# --- Utility Functions ---

def validate_dna(seq):
    """Validate that the sequence contains only A,T,G,C (case insensitive)."""
    pattern = re.compile("^[ATGCatgc]+$")
    return bool(pattern.fullmatch(seq))

def generate_random_dna(length=50):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ATGC') for _ in range(length))

def apply_point_mutation(seq, position, new_base):
    """Apply a point mutation at the given 0-based position."""
    original_base = seq[position]
    if new_base == original_base:
        return seq, original_base # Return original sequence, but indicate no change
    mutated_seq = seq[:position] + new_base + seq[position+1:]
    return mutated_seq, original_base

def apply_insertion_mutation(seq, position, insert_seq):
    """Insert a sequence at the given 0-based position."""
    mutated_seq = seq[:position] + insert_seq + seq[position:]
    return mutated_seq, insert_seq

def apply_deletion_mutation(seq, position, del_length):
    """Delete a substring of given length starting at 0-based position."""
    deleted_seq = seq[position:position+del_length]
    mutated_seq = seq[:position] + seq[position+del_length:]
    return mutated_seq, deleted_seq

def transcribe_dna_to_mrna(dna_seq):
    """Transcribes a DNA sequence (template strand) to an mRNA sequence."""
    mrna = dna_seq.replace('T', 'U')
    return mrna

def translate_mrna_to_protein(mrna_seq):
    """Translates an mRNA sequence into an amino acid sequence."""
    protein = []
    # Start translation from the first 'AUG' (Methionine) if present, otherwise from the beginning
    start_codon_idx = mrna_seq.find('AUG')
    if start_codon_idx == -1: # No start codon found, translate from beginning (non-standard)
        coding_sequence = mrna_seq
    else:
        coding_sequence = mrna_seq[start_codon_idx:]

    for i in range(0, len(coding_sequence) - len(coding_sequence) % 3, 3):
        codon = coding_sequence[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, 'X') # 'X' for unknown codon
        protein.append(amino_acid)
        if amino_acid == '*': # Stop codon
            break
    return "".join(protein)


def highlight_mutation(original_seq, mutated_seq, position, mutation_type, length=1):
    """
    Highlight the mutated parts in red.
    Returns (highlighted original, highlighted mutated) sequences as HTML strings.
    Positions are 0-based for internal function use.
    """
    def wrap_red(text):
        return f"<span class='highlight-red'>{text}</span>"

    orig_highlight = original_seq
    mut_highlight = mutated_seq

    if mutation_type == "point":
        if 0 <= position < len(original_seq) and 0 <= position < len(mutated_seq):
            orig_highlight = original_seq[:position] + wrap_red(original_seq[position]) + original_seq[position+1:]
            mut_highlight = mutated_seq[:position] + wrap_red(mutated_seq[position]) + mutated_seq[position+1:]
    elif mutation_type == "insertion":
        if 0 <= position <= len(mutated_seq) - length:
            mut_highlight = mutated_seq[:position] + wrap_red(mutated_seq[position:position+length]) + mutated_seq[position+length:]
    elif mutation_type == "deletion":
        if 0 <= position <= len(original_seq) - length:
            orig_highlight = original_seq[:position] + wrap_red(original_seq[position:position+length]) + original_seq[position+length:]

    return orig_highlight, mut_highlight

# --- Session State Initialization ---
# Initialize session state variables if they don't exist
if "dna_input" not in st.session_state:
    st.session_state.dna_input = ""
if "mutation_type" not in st.session_state:
    st.session_state.mutation_type = "point"
if "position" not in st.session_state:
    st.session_state.position = 1
if "new_base" not in st.session_state:
    st.session_state.new_base = "A"
if "insert_seq" not in st.session_state:
    st.session_state.insert_seq = ""
if "del_length" not in st.session_state:
    st.session_state.del_length = 1

# --- Streamlit UI Layout ---

st.markdown("<h1 style='color:#50FA7B;'>🧬 DNA Mutation Lab 🧪🔬</h1>", unsafe_allow_html=True)
st.markdown("<hr style='border: 1px solid #6272a4;'>", unsafe_allow_html=True)

# DNA sequence input
dna_input_area = st.text_area(
    "🧬 Input DNA Sequence (A, T, G, C only):",
    value=st.session_state.dna_input,
    height=100,
    max_chars=1000,
    help="✍️ Type or paste your DNA sequence here. It will be automatically converted to uppercase."
)
st.session_state.dna_input = dna_input_area.upper() # Always store uppercase

# Random DNA generator button
if st.button("✨ Generate Random DNA Sequence"):
    random_seq = generate_random_dna(50)
    st.session_state.dna_input = random_seq
    st.session_state.position = 1 # Reset position for new sequence
    st.rerun() # Use st.rerun() for immediate UI update

# Create two columns for mutation type and position
col1, col2 = st.columns(2)

with col1:
    mutation_type = st.selectbox(
        "⚙️ Select Mutation Type:",
        options=["point", "insertion", "deletion"],
        index=["point", "insertion", "deletion"].index(st.session_state.mutation_type),
        key="mutation_type_select",
        help="🎯 Point: Replace a single base. ➕ Insertion: Add new bases. ➖ Deletion: Remove existing bases."
    )
    st.session_state.mutation_type = mutation_type

with col2:
    max_pos = len(st.session_state.dna_input) if st.session_state.dna_input else 1
    # For insertion, position can be after the last base (len(seq) + 1)
    # For deletion and point, position cannot be beyond the sequence length
    if mutation_type == "insertion":
        current_max_pos = max_pos + 1
    else:
        current_max_pos = max_pos

    position = st.number_input(
        "📍 Mutation Position (1-based index):",
        min_value=1,
        max_value=current_max_pos,
        value=min(st.session_state.position, current_max_pos), # Adjust value if max changes
        step=1,
        key="mutation_position_input",
        help="🔢 The 1-based position where the mutation will occur."
    )
    st.session_state.position = position

# Dynamic inputs based on mutation type
if mutation_type == "point":
    new_base = st.selectbox(
        "➡️ New Base (Substitution):",
        options=["A", "T", "G", "C"],
        index=["A", "T", "G", "C"].index(st.session_state.new_base),
        key="new_base_select",
        help="Choose the nucleotide (A, T, G, or C) to replace the existing one."
    )
    st.session_state.new_base = new_base
elif mutation_type == "insertion":
    insert_seq = st.text_input(
        "➕ Sequence to Insert:",
        value=st.session_state.insert_seq,
        max_chars=100,
        help="✍️ Enter the sequence of bases (e.g., GGTCA) to insert.",
        key="insert_seq_input"
    ).upper()
    st.session_state.insert_seq = insert_seq
elif mutation_type == "deletion":
    current_dna_length = len(st.session_state.dna_input)
    # Max deletion length depends on the position chosen
    max_del_length = current_dna_length - (position - 1) if current_dna_length >= (position - 1) else 1

    del_length = st.number_input(
        "✂️ Length of Deletion:",
        min_value=1,
        max_value=max_del_length,
        value=min(st.session_state.del_length, max_del_length),
        step=1,
        key="del_length_input",
        help="📏 Specify how many bases to remove from the DNA sequence."
    )
    st.session_state.del_length = del_length

# --- Apply Mutation Button and Display Results ---

if st.button("🚀 Apply Mutation and Simulate"):
    seq = st.session_state.dna_input
    pos_idx = st.session_state.position - 1  # Convert to 0-based index for internal functions

    # Validation
    if not seq:
        st.error("🚫 Oops! Please enter a DNA sequence first to simulate mutations. 🧬")
    elif not validate_dna(seq):
        st.error("❌ Invalid DNA sequence detected! Please use only A, T, G, and C. 🧐")
    elif mutation_type == "point" and (pos_idx < 0 or pos_idx >= len(seq)):
        st.error(f"⚠️ Position {st.session_state.position} is out of bounds for a point mutation in a sequence of length {len(seq)}. Try again! 📏")
    elif mutation_type == "insertion" and (pos_idx < 0 or pos_idx > len(seq)): # Insertion can be at len(seq)
        st.error(f"⚠️ Position {st.session_state.position} is out of range for insertion in a sequence of length {len(seq)}. 🚧")
    elif mutation_type == "deletion" and (st.session_state.del_length < 1 or pos_idx < 0 or pos_idx + st.session_state.del_length > len(seq)):
        st.error(f"⚠️ Invalid deletion length or range! Cannot delete {st.session_state.del_length} bases starting at position {st.session_state.position} from a sequence of length {len(seq)}. 🤔")
    elif mutation_type == "point" and not validate_dna(st.session_state.new_base):
        st.error("❌ Invalid new base for point mutation. Please use A, T, G, or C. 💡")
    elif mutation_type == "insertion" and (not st.session_state.insert_seq or not validate_dna(st.session_state.insert_seq)):
        st.error("❌ Invalid sequence to insert. It must contain only A, T, G, and C. 🧬")
    else:
        mutated_seq = ""
        consequence_summary = ""
        biological_impact_details = "" # This will be populated in the next step

        mutation_length_for_highlight = 1 # Default for point, will be overridden

        if mutation_type == "point":
            mutated_seq, original_base = apply_point_mutation(seq, pos_idx, st.session_state.new_base)
            if mutated_seq == seq:
                consequence_summary = f"No actual change occurred! 🔄 The base '{original_base}' was already at position {st.session_state.position}. Sequence remains identical. ✅"
            else:
                consequence_summary = f"Point mutation: Substituted '{original_base}' with '{st.session_state.new_base}' at position {st.session_state.position}. 🎯"
            mutation_length_for_highlight = 1
        elif mutation_type == "insertion":
            mutated_seq, inserted = apply_insertion_mutation(seq, pos_idx, st.session_state.insert_seq)
            consequence_summary = f"Insertion: Inserted '{inserted}' (length {len(inserted)}) at position {st.session_state.position}. ➕"
            mutation_length_for_highlight = len(st.session_state.insert_seq)
        elif mutation_type == "deletion":
            mutated_seq, deleted = apply_deletion_mutation(seq, pos_idx, st.session_state.del_length)
            consequence_summary = f"Deletion: Removed '{deleted}' (length {st.session_state.del_length}) from position {st.session_state.position} to {st.session_state.position + st.session_state.del_length - 1}. ➖"
            mutation_length_for_highlight = st.session_state.del_length

        # --- Perform Transcription and Translation ---
        original_mrna = transcribe_dna_to_mrna(seq)
        original_protein = translate_mrna_to_protein(original_mrna)

        mutated_mrna = transcribe_dna_to_mrna(mutated_seq)
        mutated_protein = translate_mrna_to_protein(mutated_mrna)


        # --- Highlight and display results ---
        orig_highlight_dna, mut_highlight_dna = highlight_mutation(
            seq,
            mutated_seq,
            pos_idx, # Pass 0-based index to highlight function
            mutation_type,
            length=mutation_length_for_highlight
        )

        st.markdown("<hr style='border: 1px dashed #44475a;'>", unsafe_allow_html=True) # Dotted separator

        # Display DNA sequences
        st.markdown("### 🧬 Original DNA Sequence:")
        st.markdown(
            f"<div class='output-box'>{orig_highlight_dna}</div>",
            unsafe_allow_html=True
        )

        st.markdown("### 🧬 Mutated DNA Sequence:")
        st.markdown(
            f"<div class='output-box'>{mut_highlight_dna}</div>",
            unsafe_allow_html=True
        )

        st.markdown("---") # Separator for RNA/Protein

        # Display mRNA sequences
        st.markdown("### 📊 mRNA Sequences (Transcribed from DNA):")
        st.markdown(f"**Original mRNA:** `{original_mrna}`")
        st.markdown(f"**Mutated mRNA:** `{mutated_mrna}`")


        # Display Protein sequences
        st.markdown("### 🧪 Protein Sequences (Translated from mRNA):")
        st.markdown(f"**Original Protein:** `{original_protein}`")
        st.markdown(f"**Mutated Protein:** `{mutated_protein}`")

        st.markdown("<hr style='border: 1px dashed #44475a;'>", unsafe_allow_html=True) # Dotted separator

        st.markdown("### 🔬 Summary of Mutation:")
        st.markdown(
            f"<div style='background-color:#222;color:#ff5555;padding:15px;border-radius:10px;font-weight:bold; border: 1px solid #ff5555;'>{consequence_summary}</div>",
            unsafe_allow_html=True
        )

        # The biological impact details will be populated in the next step based on mutation classification
        st.markdown("### 💡 Potential Biological Impact:")
        with st.expander("Click to learn more about the biological impact of this mutation type"):
            # This section will be dynamic based on classification in the next step
            st.markdown(
                f"<div class='biological-impact-box'>{biological_impact_details if biological_impact_details else 'Details will appear here after mutation classification is implemented!'}</div>",
                unsafe_allow_html=True
            )
