import random
import streamlit as st
from Bio.Seq import Seq
if "random_dna" not in st.session_state:
    st.session_state.random_dna = ""

if "mutation_result" not in st.session_state:
    st.session_state.mutation_result = ""
if "dna_input" not in st.session_state:
    st.session_state.dna_input = ""

# Add all others you plan to use...


# Set page config and custom styles for background and output boxes
st.set_page_config(page_title="DNA Mutation Simulator", layout="wide")

# Custom CSS for black background and gray output boxes with readable text
st.markdown(
    """
    <style>
    .main {
        background-color: #000000;
        color: white;
        font-family: 'Courier New', Courier, monospace;
    }
    .output-box {
        background-color: #2f2f2f;
        padding: 10px;
        border-radius: 8px;
        color: #f0f0f0;
        font-family: monospace;
        white-space: pre-wrap;
        overflow-x: auto;
    }
    .highlight {
        color: #ff5555;
        font-weight: bold;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

def validate_dna(sequence):
    return all(base.upper() in 'ATGC' for base in sequence)

def get_codon_and_aa(sequence, nucleotide_pos):
    if not (0 <= nucleotide_pos < len(sequence)):
        return None, None

    codon_start_pos = (nucleotide_pos // 3) * 3
    codon_end_pos = codon_start_pos + 3

    codon = sequence[codon_start_pos:codon_end_pos]

    amino_acid = str(Seq(codon).translate(table=1))

    return codon, amino_acid

def highlight_mutation(original_seq, mutated_seq, pos, mutation_type, change_len=1):
    # Highlight mutated parts in red using HTML span with 'highlight' class
    def highlight_span(text):
        return f'<span class="highlight">{text}</span>'

    if mutation_type == 'point':
        highlighted_orig = (
            original_seq[:pos]
            + highlight_span(original_seq[pos])
            + original_seq[pos + 1 :]
        )
        highlighted_mut = (
            mutated_seq[:pos]
            + highlight_span(mutated_seq[pos])
            + mutated_seq[pos + 1 :]
        )
    elif mutation_type == 'insertion':
        # Highlight inserted sequence in mutated sequence
        highlighted_orig = original_seq
        highlighted_mut = (
            mutated_seq[:pos]
            + highlight_span(mutated_seq[pos : pos + change_len])
            + mutated_seq[pos + change_len :]
        )
    elif mutation_type == 'deletion':
        # Highlight deleted sequence in original sequence
        highlighted_orig = (
            original_seq[:pos]
            + highlight_span(original_seq[pos : pos + change_len])
            + original_seq[pos + change_len :]
        )
        highlighted_mut = mutated_seq
    else:
        highlighted_orig = original_seq
        highlighted_mut = mutated_seq

    return highlighted_orig, highlighted_mut

def apply_point_mutation(dna_seq_str, pos, new_base):
    length = len(dna_seq_str)
    if not (0 <= pos < length):
        return None, None

    original_base = dna_seq_str[pos]

    mutated_seq = dna_seq_str[:pos] + new_base + dna_seq_str[pos + 1 :]
    return mutated_seq, original_base

def apply_insertion_mutation(dna_seq_str, pos, inserted_seq):
    if not inserted_seq:
        return None, None

    length = len(dna_seq_str)
    if not (0 <= pos <= length):
        return None, None

    mutated_seq = dna_seq_str[:pos] + inserted_seq + dna_seq_str[pos:]
    return mutated_seq, inserted_seq

def apply_deletion_mutation(dna_seq_str, pos, delete_length):
    if delete_length <= 0:
        return None, None

    length = len(dna_seq_str)
    if length < delete_length:
        return None, None

    if not (0 <= pos < length and (pos + delete_length) <= length):
        return None, None

    deleted_seq = dna_seq_str[pos : pos + delete_length]

    mutated_seq = dna_seq_str[:pos] + dna_seq_str[pos + delete_length :]
    return mutated_seq, deleted_seq

# Start of Streamlit app UI
st.title("🚀 DNA Mutation Simulator!")
st.markdown(
    """
    This tool helps you understand how tiny changes in DNA can lead to big differences in proteins.
    Let's begin by providing a DNA sequence to work with!
    """
)

# DNA input: manual or random
dna_input_mode = st.radio(
    "Choose DNA input method:",
    options=["Manual Input", "Random Sequence"],
    index=0,
)

original_dna_str = ""

if dna_input_mode == "Manual Input":
    dna_input = st.text_area(
        "🚀 Enter your DNA sequence (A, T, G, C only):",
        height=80,
        max_chars=500,
        placeholder="Example: ATGCGTGACTGACTGACGTA",
    ).upper().strip()

    if dna_input:
        if validate_dna(dna_input):
            original_dna_str = dna_input
        else:
            st.error("❌ Invalid DNA sequence! Only A, T, G, C allowed.")
else:
    length = st.slider(
        "🔢 Select length of random DNA sequence:",
        min_value=10,
        max_value=500,
        value=50,
        step=1,
    )
    if st.button("Generate Random DNA Sequence"):
        original_dna_str = ''.join(random.choice('ATGC') for _ in range(length))
        st.success(f"✅ Random DNA sequence of length {length} generated!")

if original_dna_str:
    st.markdown(f"🧬 **Original DNA Sequence (Length: {len(original_dna_str)} bases):**")
    st.markdown(f'<div class="output-box">{original_dna_str}</div>', unsafe_allow_html=True)

    # Translate to protein if possible
    if len(original_dna_str) >= 3:
        try:
            original_protein = str(Seq(original_dna_str).translate(table=1, to_stop=True))
            st.markdown(f"➡️ **Original Protein Sequence (Length: {len(original_protein)} amino acids):**")
            st.markdown(f'<div class="output-box">{original_protein}</div>', unsafe_allow_html=True)
        except Exception:
            st.warning("⚠️ Warning: DNA might be too short or irregular for full protein translation.")
    else:
        st.warning("⚠️ Warning: DNA too short (<3 bases) for protein translation.")
# --- Mutation Section ---

st.header("🧪 Apply Mutations to DNA Sequence")

mutation_type = st.selectbox(
    "Choose Mutation Type:",
    options=["Point Mutation 🧬", "Insertion ➕", "Deletion ➖"],
)

mutated_dna_str = ""
mutation_position = None
mutation_info = ""
highlighted_original = ""
highlighted_mutated = ""

if original_dna_str:
    seq_len = len(original_dna_str)

    if mutation_type == "Point Mutation 🧬":
        mutation_position = st.number_input(
            "Enter mutation position (0-based index):",
            min_value=0,
            max_value=seq_len - 1,
            value=0,
            step=1,
        )
        original_base = original_dna_str[mutation_position]
        new_base = st.selectbox(
            f"Replace base '{original_base}' at position {mutation_position} with:",
            options=[b for b in 'ATGC' if b != original_base],
        )

        if st.button("Apply Point Mutation"):
            mutated_seq, orig_base = apply_point_mutation(
                original_dna_str, mutation_position, new_base
            )
            if mutated_seq:
                mutated_dna_str = mutated_seq
                highlighted_original, highlighted_mutated = highlight_mutation(
                    original_dna_str, mutated_dna_str, mutation_position, 'point'
                )
                mutation_info = f"Point mutation: position {mutation_position}, {orig_base} → {new_base}."

    elif mutation_type == "Insertion ➕":
        mutation_position = st.number_input(
            "Enter insertion position (0-based index):",
            min_value=0,
            max_value=seq_len,
            value=0,
            step=1,
        )
        inserted_seq = st.text_input(
            "Enter DNA sequence to insert (A, T, G, C only):",
            max_chars=50,
        ).upper().strip()

        if st.button("Apply Insertion"):
            if validate_dna(inserted_seq) and inserted_seq:
                mutated_seq, inserted = apply_insertion_mutation(
                    original_dna_str, mutation_position, inserted_seq
                )
                if mutated_seq:
                    mutated_dna_str = mutated_seq
                    highlighted_original, highlighted_mutated = highlight_mutation(
                        original_dna_str, mutated_dna_str, mutation_position, 'insertion', len(inserted_seq)
                    )
                    mutation_info = f"Insertion: position {mutation_position}, inserted '{inserted_seq}'."
            else:
                st.error("Invalid insertion sequence! Use only A, T, G, C.")

    elif mutation_type == "Deletion ➖":
        mutation_position = st.number_input(
            "Enter deletion start position (0-based index):",
            min_value=0,
            max_value=seq_len - 1,
            value=0,
            step=1,
        )
        max_del_len = seq_len - mutation_position
        delete_length = st.number_input(
            "Enter number of bases to delete:",
            min_value=1,
            max_value=max_del_len,
            value=1,
            step=1,
        )

        if st.button("Apply Deletion"):
            mutated_seq, deleted_seq = apply_deletion_mutation(
                original_dna_str, mutation_position, delete_length
            )
            if mutated_seq:
                mutated_dna_str = mutated_seq
                highlighted_original, highlighted_mutated = highlight_mutation(
                    original_dna_str, mutated_dna_str, mutation_position, 'deletion', delete_length
                )
                mutation_info = f"Deletion: position {mutation_position}, deleted '{deleted_seq}'."

    # Display results after mutation
    if mutated_dna_str:
        st.markdown("### 🧬 DNA Sequences with Mutation Highlighted")
        st.markdown(
            "**Original Sequence:**",
        )
        st.markdown(
            f'<div class="output-box" style="color:white;">{highlighted_original}</div>',
            unsafe_allow_html=True,
        )
        st.markdown(
            "**Mutated Sequence:**",
        )
        st.markdown(
            f'<div class="output-box" style="color:white;">{highlighted_mutated}</div>',
            unsafe_allow_html=True,
        )

        # Translate mutated DNA
        try:
            mutated_protein = str(Seq(mutated_dna_str).translate(table=1, to_stop=True))
            original_protein = str(Seq(original_dna_str).translate(table=1, to_stop=True))
            st.markdown("### 🧫 Protein Sequences")
            st.markdown(f"**Original Protein:**")
            st.markdown(f'<div class="output-box">{original_protein}</div>', unsafe_allow_html=True)
            st.markdown(f"**Mutated Protein:**")
            st.markdown(f'<div class="output-box">{mutated_protein}</div>', unsafe_allow_html=True)

            # Show mutation info and consequence
            st.markdown(f"### ℹ️ Mutation Info")
            st.info(mutation_info)

            # Check protein change
            if mutated_protein == original_protein:
                st.success("🟢 Mutation is synonymous: no change in protein sequence.")
            else:
                st.warning("🔴 Mutation is non-synonymous: protein sequence changed!")

        except Exception as e:
            st.error(f"⚠️ Error translating DNA sequences: {str(e)}")
# --- Utility Functions ---

def validate_dna(seq):
    return all(base in "ATGC" for base in seq)

def apply_point_mutation(seq, position, new_base):
    if 0 <= position < len(seq) and new_base in "ATGC":
        original_base = seq[position]
        mutated_seq = seq[:position] + new_base + seq[position + 1:]
        return mutated_seq, original_base
    return None, None

def apply_insertion_mutation(seq, position, insert_seq):
    if 0 <= position <= len(seq) and validate_dna(insert_seq):
        mutated_seq = seq[:position] + insert_seq + seq[position:]
        return mutated_seq, insert_seq
    return None, None

def apply_deletion_mutation(seq, position, length):
    if 0 <= position < len(seq) and (position + length) <= len(seq):
        deleted_seq = seq[position:position+length]
        mutated_seq = seq[:position] + seq[position+length:]
        return mutated_seq, deleted_seq
    return None, None

def highlight_mutation(original, mutated, position, mutation_type, length=1):
    if mutation_type == "point":
        original_highlighted = (
            original[:position]
            + f"<span style='color:red;font-weight:bold'>{original[position]}</span>"
            + original[position+1:]
        )
        mutated_highlighted = (
            mutated[:position]
            + f"<span style='color:red;font-weight:bold'>{mutated[position]}</span>"
            + mutated[position+1:]
        )
    elif mutation_type == "insertion":
        original_highlighted = original
        mutated_highlighted = (
            mutated[:position]
            + f"<span style='color:red;font-weight:bold'>{mutated[position:position+length]}</span>"
            + mutated[position+length:]
        )
    elif mutation_type == "deletion":
        original_highlighted = (
            original[:position]
            + f"<span style='color:red;font-weight:bold'>{original[position:position+length]}</span>"
            + original[position+length:]
        )
        mutated_highlighted = mutated
    else:
        original_highlighted = original
        mutated_highlighted = mutated

    return original_highlighted, mutated_highlighted

# --- Random DNA Generator ---

def generate_random_dna(length=60):
    import random
    return ''.join(random.choices("ATGC", k=length))

st.sidebar.markdown("### 🧬 Random DNA Generator")
random_dna_len = st.sidebar.slider("Random DNA Length", 10, 300, 60)
if st.sidebar.button("🎲 Generate Random DNA"):
    st.session_state.dna_input = generate_random_dna(random_dna_len)
# --- Mutation Logic and Display ---

if st.session_state.dna_input:
    seq = st.session_state.dna_input.upper()
    st.markdown("<h3 style='color:white;'>🔬 Mutation Result</h3>", unsafe_allow_html=True)

    if not validate_dna(seq):
        st.error("Invalid DNA sequence. Please enter a valid sequence using only A, T, G, and C.")
    else:
        mutated_seq = seq
        consequence = ""
        if mutation_type == "point":
            mutated_seq, original_base = apply_point_mutation(seq, position, new_base)
            consequence = f"Substituted {original_base} with {new_base} at position {position}."
        elif mutation_type == "insertion":
            mutated_seq, inserted = apply_insertion_mutation(seq, position, insert_seq)
            consequence = f"Inserted {inserted} at position {position}."
        elif mutation_type == "deletion":
            mutated_seq, deleted = apply_deletion_mutation(seq, position, del_length)
            consequence = f"Deleted {deleted} from position {position} to {position + del_length - 1}."

        if mutated_seq:
            orig_highlight, mut_highlight = highlight_mutation(
                seq, mutated_seq, position, mutation_type,
                length=(len(insert_seq) if mutation_type == "insertion" else del_length)
            )

            st.markdown("#### 🧬 Original Sequence")
            st.markdown(f"<div style='background-color:#444;padding:10px;border-radius:8px;'>{orig_highlight}</div>", unsafe_allow_html=True)

            st.markdown("#### 🧬 Mutated Sequence")
            st.markdown(f"<div style='background-color:#444;padding:10px;border-radius:8px;'>{mut_highlight}</div>", unsafe_allow_html=True)

            st.markdown(f"#### ⚠️ Biological Consequence")
            st.markdown(f"<div style='background-color:#333;color:white;padding:10px;border-radius:8px;'>{consequence}</div>", unsafe_allow_html=True)
        else:
            st.error("⚠️ Mutation failed due to invalid parameters.")
