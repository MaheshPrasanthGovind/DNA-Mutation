import streamlit as st
import random
import re

# Set page config for dark theme
st.set_page_config(page_title="🧬 DNA Mutation Simulator", page_icon="🧬", layout="wide")

# Custom CSS for background and output box colors
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
    }
    .highlight-red {
        color: red;
        font-weight: bold;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# --- Utility Functions ---

def validate_dna(seq):
    """Validate that the sequence contains only A,T,G,C (case insensitive)."""
    pattern = re.compile("^[ATGCatgc]+$")
    return bool(pattern.fullmatch(seq))

def generate_random_dna(length=50):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ATGC') for _ in range(length))

def apply_point_mutation(seq, position, new_base):
    """Apply a point mutation at the given 1-based position."""
    # Positions are 1-based; convert to 0-based for string manipulation
    index = position - 1
    if index < 0 or index >= len(seq):
        return None, None  # invalid position
    original_base = seq[index]
    mutated_seq = seq[:index] + new_base + seq[index+1:]
    return mutated_seq, original_base

def apply_insertion_mutation(seq, position, insert_seq):
    """Insert a sequence at the given 1-based position."""
    index = position - 1
    if index < 0 or index > len(seq):
        return None, None
    mutated_seq = seq[:index] + insert_seq + seq[index:]
    return mutated_seq, insert_seq

def apply_deletion_mutation(seq, position, del_length):
    """Delete a substring of given length starting at 1-based position."""
    index = position - 1
    if index < 0 or index + del_length > len(seq):
        return None, None
    deleted_seq = seq[index:index+del_length]
    mutated_seq = seq[:index] + seq[index+del_length:]
    return mutated_seq, deleted_seq

def highlight_mutation(orig_seq, mutated_seq, position, mutation_type, length=1):
    """
    Highlight the mutated parts in red.
    Return (highlighted original, highlighted mutated) sequences as HTML strings.
    """
    def wrap_red(text):
        return f"<span class='highlight-red'>{text}</span>"

    pos = position - 1  # zero-based

    if mutation_type == "point":
        orig_highlight = (
            orig_seq[:pos] + wrap_red(orig_seq[pos]) + orig_seq[pos+1:]
        )
        mut_highlight = (
            mutated_seq[:pos] + wrap_red(mutated_seq[pos]) + mutated_seq[pos+1:]
        )
    elif mutation_type == "insertion":
        orig_highlight = orig_seq
        mut_highlight = (
            mutated_seq[:pos] + wrap_red(mutated_seq[pos:pos+length]) + mutated_seq[pos+length:]
        )
    elif mutation_type == "deletion":
        orig_highlight = (
            orig_seq[:pos] + wrap_red(orig_seq[pos:pos+length]) + orig_seq[pos+length:]
        )
        mut_highlight = mutated_seq
    else:
        orig_highlight = orig_seq
        mut_highlight = mutated_seq

    return orig_highlight, mut_highlight

# --- Initialize session state variables ---

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

# --- Sidebar Input ---

st.sidebar.header("🧬 DNA Mutation Simulator")

# DNA sequence input
dna_input = st.sidebar.text_area(
    "Enter DNA sequence (A, T, G, C only):",
    value=st.session_state.dna_input,
    height=120,
    max_chars=1000,
    help="Enter your DNA sequence here."
)

if dna_input:
    st.session_state.dna_input = dna_input.upper()

# Random DNA generator button
if st.sidebar.button("🔀 Generate Random DNA"):
    random_seq = generate_random_dna(50)
    st.session_state.dna_input = random_seq
    st.experimental_rerun()

# Mutation type selection
mutation_type = st.sidebar.selectbox(
    "Select Mutation Type:",
    options=["point", "insertion", "deletion"],
    index=["point", "insertion", "deletion"].index(st.session_state.mutation_type),
)
st.session_state.mutation_type = mutation_type

# Position input (1-based)
position = st.sidebar.number_input(
    "Mutation Position (1-based index):",
    min_value=1,
    max_value=10000,
    value=st.session_state.position,
    step=1,
    help="Position in the DNA sequence to mutate."
)
st.session_state.position = position

# Dynamic inputs based on mutation type
if mutation_type == "point":
    new_base = st.sidebar.selectbox(
        "New Base (Substitution):",
        options=["A", "T", "G", "C"],
        index=["A", "T", "G", "C"].index(st.session_state.new_base),
    )
    st.session_state.new_base = new_base

elif mutation_type == "insertion":
    insert_seq = st.sidebar.text_input(
        "Sequence to Insert:",
        value=st.session_state.insert_seq,
        max_chars=100,
        help="Sequence of bases to insert at position."
    ).upper()
    st.session_state.insert_seq = insert_seq

elif mutation_type == "deletion":
    del_length = st.sidebar.number_input(
        "Length of Deletion:",
        min_value=1,
        max_value=1000,
        value=st.session_state.del_length,
        step=1,
        help="Number of bases to delete starting at position."
    )
    st.session_state.del_length = del_length

# --- Chunk 2 ---

import random
import streamlit as st

# Function to generate a random DNA sequence
def generate_random_dna(length=50):
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(length))

# Function to validate the DNA sequence
def validate_dna(seq):
    valid_bases = {'A', 'T', 'G', 'C'}
    return all(base in valid_bases for base in seq)

# Function to apply point mutation
def apply_point_mutation(seq, position, new_base):
    original_base = seq[position]
    mutated_seq = seq[:position] + new_base + seq[position+1:]
    return mutated_seq, original_base

# Function to apply insertion mutation
def apply_insertion_mutation(seq, position, insert_seq):
    mutated_seq = seq[:position] + insert_seq + seq[position:]
    return mutated_seq, insert_seq

# Function to apply deletion mutation
def apply_deletion_mutation(seq, position, length):
    deleted_seq = seq[position:position+length]
    mutated_seq = seq[:position] + seq[position+length:]
    return mutated_seq, deleted_seq

# Function to highlight mutations in the sequence
def highlight_mutation(original_seq, mutated_seq, position, mutation_type, length=1):
    # Helper to wrap mutated parts with red span
    def red(text):
        return f"<span style='color:red;font-weight:bold;'>{text}</span>"

    if mutation_type == "point":
        orig_highlight = (
            original_seq[:position] + red(original_seq[position]) + original_seq[position+1:]
        )
        mut_highlight = (
            mutated_seq[:position] + red(mutated_seq[position]) + mutated_seq[position+1:]
        )
    elif mutation_type == "insertion":
        orig_highlight = original_seq
        mut_highlight = (
            mutated_seq[:position] + red(mutated_seq[position:position+length]) + mutated_seq[position+length:]
        )
    elif mutation_type == "deletion":
        orig_highlight = (
            original_seq[:position] + red(original_seq[position:position+length]) + original_seq[position+length:]
        )
        mut_highlight = mutated_seq
    else:
        orig_highlight = original_seq
        mut_highlight = mutated_seq

    return orig_highlight, mut_highlight

# Session state initialization for mutation controls
if 'dna_input' not in st.session_state:
    st.session_state.dna_input = ''
if 'mutation_type' not in st.session_state:
    st.session_state.mutation_type = 'point'
if 'position' not in st.session_state:
    st.session_state.position = 0
if 'new_base' not in st.session_state:
    st.session_state.new_base = 'A'
if 'insert_seq' not in st.session_state:
    st.session_state.insert_seq = ''
if 'del_length' not in st.session_state:
    st.session_state.del_length = 1


# ------------------- Chunk 3 -------------------

# Function to highlight mutated DNA parts in red 🔴
def highlight_mutation(original_seq, mutated_seq, position, mutation_type, length=1):
    """
    Highlights mutations in the sequences:
    - original_seq: Original DNA sequence (string)
    - mutated_seq: Mutated DNA sequence (string)
    - position: mutation position (int, 1-based index)
    - mutation_type: 'point', 'insertion', or 'deletion'
    - length: length of insertion/deletion
    Returns highlighted HTML strings for original and mutated sequences
    """
    start = position - 1  # Convert to 0-based index

    if mutation_type == "point":
        orig_highlight = (
            original_seq[:start]
            + f"<span style='color:red;font-weight:bold;'>{original_seq[start]}</span>"
            + original_seq[start + 1 :]
        )
        mut_highlight = (
            mutated_seq[:start]
            + f"<span style='color:red;font-weight:bold;'>{mutated_seq[start]}</span>"
            + mutated_seq[start + 1 :]
        )
    elif mutation_type == "insertion":
        # Highlight insertion region in mutated sequence
        orig_highlight = original_seq
        mut_highlight = (
            mutated_seq[:start]
            + f"<span style='color:red;font-weight:bold;'>{mutated_seq[start:start+length]}</span>"
            + mutated_seq[start + length :]
        )
    elif mutation_type == "deletion":
        # Highlight deleted region in original sequence
        orig_highlight = (
            original_seq[:start]
            + f"<span style='color:red;font-weight:bold;'>{original_seq[start:start+length]}</span>"
            + original_seq[start + length :]
        )
        mut_highlight = mutated_seq
    else:
        orig_highlight = original_seq
        mut_highlight = mutated_seq

    return orig_highlight, mut_highlight


# Function to validate DNA sequence 🧬
def validate_dna(seq):
    """Validate that sequence contains only A, T, G, C (case-insensitive)"""
    valid_nucleotides = set("ATGC")
    return all(nucleotide in valid_nucleotides for nucleotide in seq.upper())


# Function to generate random DNA sequence 🔀
def generate_random_dna(length):
    """Generate a random DNA sequence of specified length"""
    nucleotides = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(nucleotides) for _ in range(length))


# Mutation Functions (from chunk 2 but here for reference)
def apply_point_mutation(seq, position, new_base):
    original_base = seq[position - 1]
    if new_base == original_base:
        return None, None  # No mutation needed
    mutated_seq = seq[:position - 1] + new_base + seq[position:]
    return mutated_seq, original_base


def apply_insertion_mutation(seq, position, insert_seq):
    mutated_seq = seq[:position] + insert_seq + seq[position:]
    return mutated_seq, insert_seq


def apply_deletion_mutation(seq, position, length):
    deleted_seq = seq[position - 1: position - 1 + length]
    mutated_seq = seq[:position - 1] + seq[position - 1 + length:]
    return mutated_seq, deleted_seq


# Streamlit UI logic for mutation visualization and inputs 👇
if 'dna_input' not in st.session_state:
    st.session_state.dna_input = ""

if 'mutation_type' not in st.session_state:
    st.session_state.mutation_type = "point"

if 'position' not in st.session_state:
    st.session_state.position = 1

if 'new_base' not in st.session_state:
    st.session_state.new_base = "A"

if 'insert_seq' not in st.session_state:
    st.session_state.insert_seq = ""

if 'del_length' not in st.session_state:
    st.session_state.del_length = 1


st.markdown("<h2 style='color:#ff79c6;'>🧬 DNA Mutation Simulator 🔬</h2>", unsafe_allow_html=True)

# DNA sequence input or random generation
dna_input = st.text_area("Enter DNA Sequence (A, T, G, C only):", st.session_state.dna_input, height=80, max_chars=1000)
st.session_state.dna_input = dna_input.upper()

col1, col2 = st.columns(2)
with col1:
    mutation_type = st.selectbox(
        "Choose Mutation Type 🧪:",
        ("point", "insertion", "deletion"),
        index=["point", "insertion", "deletion"].index(st.session_state.mutation_type),
    )
    st.session_state.mutation_type = mutation_type
with col2:
    position = st.number_input(
        "Position (1-based) 🧮:",
        min_value=1,
        max_value=len(st.session_state.dna_input) if st.session_state.dna_input else 1,
        value=st.session_state.position,
        step=1,
    )
    st.session_state.position = position

if mutation_type == "point":
    new_base = st.selectbox(
        "New Base (for point mutation) 🔄:",
        ("A", "T", "G", "C"),
        index=["A", "T", "G", "C"].index(st.session_state.new_base),
    )
    st.session_state.new_base = new_base
elif mutation_type == "insertion":
    insert_seq = st.text_input(
        "Sequence to Insert ➕:",
        value=st.session_state.insert_seq,
        max_chars=50,
    ).upper()
    st.session_state.insert_seq = insert_seq
elif mutation_type == "deletion":
    del_length = st.number_input(
        "Deletion Length ➖:",
        min_value=1,
        max_value=len(st.session_state.dna_input) - position + 1 if st.session_state.dna_input else 1,
        value=st.session_state.del_length,
        step=1,
    )
    st.session_state.del_length = del_length

# Process mutation when button pressed
if st.button("Apply Mutation 🧬"):
    seq = st.session_state.dna_input
    if not seq:
        st.error("⚠️ Please enter a DNA sequence first!")
    elif not validate_dna(seq):
        st.error("🚫 Invalid DNA sequence. Use only A, T, G, and C.")
    elif position < 1 or position > len(seq):
        st.error("⚠️ Position out of range!")
    elif mutation_type == "point" and new_base not in ["A", "T", "G", "C"]:
        st.error("⚠️ Invalid base for point mutation.")
    elif mutation_type == "insertion" and not insert_seq:
        st.error("⚠️ Please enter sequence to insert.")
    elif mutation_type == "deletion" and (del_length < 1 or del_length > len(seq) - position + 1):
        st.error("⚠️ Invalid deletion length.")
    else:
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
                seq,
                mutated_seq,
                position,
                mutation_type,
                length=len(insert_seq) if mutation_type == "insertion" else del_length,
            )
            st.markdown("#### 🧬 Original Sequence")
            st.markdown(f"<div style='background-color:#444;padding:10px;border-radius:8px;color:white;font-family:monospace;'>{orig_highlight}</div>", unsafe_allow_html=True)

            st.markdown("#### 🧬 Mutated Sequence")
            st.markdown(f"<div style='background-color:#444;padding:10px;border-radius:8px;color:white;font-family:monospace;'>{mut_highlight}</div>", unsafe_allow_html=True)

            st.markdown(f"#### ⚠️ Biological Consequence")
            st.markdown(f"<div style='background-color:#333;color:#ff5555;padding:10px;border-radius:8px;font-family:monospace;'>{consequence}</div>", unsafe_allow_html=True)
        else:
            st.error("⚠️ Mutation failed due to invalid parameters.")
# ------------------------- Chunk 4: Output display, mutation handling, and UI polish -------------------------

import streamlit as st

def highlight_mutation(original_seq, mutated_seq, position, mutation_type, length=1):
    """
    Highlights mutated parts in red within the DNA sequences for display.
    """
    def highlight_seq(seq, start, length):
        return (
            seq[:start] +
            f"<span style='color:red;font-weight:bold;'>{seq[start:start+length]}</span>" +
            seq[start+length:]
        )

    if mutation_type == "point":
        orig_highlight = highlight_seq(original_seq, position, 1)
        mut_highlight = highlight_seq(mutated_seq, position, 1)
    elif mutation_type == "insertion":
        orig_highlight = original_seq  # No highlight in original for insertion
        mut_highlight = highlight_seq(mutated_seq, position, length)
    elif mutation_type == "deletion":
        orig_highlight = highlight_seq(original_seq, position, length)
        mut_highlight = mutated_seq  # No highlight in mutated for deletion
    else:
        orig_highlight = original_seq
        mut_highlight = mutated_seq

    return orig_highlight, mut_highlight

# UI & Mutation Processing

st.markdown("<h1 style='color:#50FA7B;'>🧬 DNA Mutation Simulator 🧬</h1>", unsafe_allow_html=True)
st.markdown("<hr style='border: 1px solid #6272a4;'>", unsafe_allow_html=True)

# Input DNA sequence box with gray background
dna_input = st.text_area(
    "🧬 Enter DNA Sequence (A, T, G, C only):",
    value=st.session_state.get("dna_input", ""),
    height=100,
    max_chars=1000,
    key="dna_input",
    help="Example: ATGCGTACGTA",
)

# Mutation type selector with emojis
mutation_type = st.selectbox(
    "⚙️ Select Mutation Type:",
    options=["point", "insertion", "deletion"],
    index=st.session_state.get("mutation_type", 0),
    key="mutation_type",
    help="Point: Single base substitution; Insertion: Add bases; Deletion: Remove bases."
)

# Position input
position = st.number_input(
    "📍 Mutation Position (1-based index):",
    min_value=1,
    max_value=len(dna_input) if dna_input else 1,
    value=st.session_state.get("position", 1),
    step=1,
    key="position"
)

# Additional inputs depending on mutation type
if mutation_type == "point":
    new_base = st.selectbox(
        "🔄 New Base:",
        options=["A", "T", "G", "C"],
        index=st.session_state.get("new_base_index", 0),
        key="new_base"
    )
elif mutation_type == "insertion":
    insert_seq = st.text_input(
        "➕ Insert Sequence:",
        value=st.session_state.get("insert_seq", ""),
        max_chars=50,
        key="insert_seq"
    )
elif mutation_type == "deletion":
    max_del_length = len(dna_input) - position + 1 if dna_input else 1
    del_length = st.number_input(
        "➖ Deletion Length:",
        min_value=1,
        max_value=max_del_length,
        value=st.session_state.get("del_length", 1),
        step=1,
        key="del_length"
    )

# Store input states back to session to preserve UI state on rerun
st.session_state.dna_input = dna_input
st.session_state.mutation_type = mutation_type
st.session_state.position = position
if mutation_type == "point":
    st.session_state.new_base_index = ["A", "T", "G", "C"].index(new_base)
elif mutation_type == "insertion":
    st.session_state.insert_seq = insert_seq
elif mutation_type == "deletion":
    st.session_state.del_length = del_length

# Perform mutation and display result when button is pressed
if st.button("🚀 Generate Mutation"):

    seq = dna_input.upper()
    pos_idx = position - 1  # Convert to 0-based index

    # Validation
    valid_bases = {"A", "T", "G", "C"}
    if not seq or any(base not in valid_bases for base in seq):
        st.error("❌ Invalid DNA sequence. Use only A, T, G, C.")
    elif pos_idx < 0 or pos_idx >= len(seq):
        st.error("❌ Position is out of range.")
    elif mutation_type == "point" and new_base not in valid_bases:
        st.error("❌ Invalid new base for point mutation.")
    elif mutation_type == "insertion" and (not insert_seq or any(b not in valid_bases for b in insert_seq.upper())):
        st.error("❌ Invalid insert sequence. Use only A, T, G, C.")
    elif mutation_type == "deletion" and (del_length < 1 or pos_idx + del_length > len(seq)):
        st.error("❌ Invalid deletion length or range.")
    else:
        # Apply mutations
        if mutation_type == "point":
            mutated_seq, original_base = apply_point_mutation(seq, pos_idx, new_base)
            bio_consequence = f"⚠️ Substituted {original_base} with {new_base} at position {position}."
        elif mutation_type == "insertion":
            mutated_seq, inserted = apply_insertion_mutation(seq, pos_idx, insert_seq.upper())
            bio_consequence = f"⚠️ Inserted '{inserted}' at position {position}."
        else:  # deletion
            mutated_seq, deleted = apply_deletion_mutation(seq, pos_idx, del_length)
            bio_consequence = f"⚠️ Deleted '{deleted}' from position {position} to {position + del_length - 1}."

        # Highlight mutated parts
        orig_highlight, mut_highlight = highlight_mutation(
            seq, mutated_seq, pos_idx, mutation_type,
            length=len(insert_seq) if mutation_type == "insertion" else (del_length if mutation_type == "deletion" else 1)
        )

        # Display sequences with highlights
        st.markdown("#### 🧬 Original Sequence")
        st.markdown(
            f"<div style='background-color:#444;color:white;padding:15px;border-radius:10px;overflow-wrap: break-word;font-family: monospace;'>{orig_highlight}</div>",
            unsafe_allow_html=True
        )

        st.markdown("#### 🧬 Mutated Sequence")
        st.markdown(
            f"<div style='background-color:#444;color:white;padding:15px;border-radius:10px;overflow-wrap: break-word;font-family: monospace;'>{mut_highlight}</div>",
            unsafe_allow_html=True
        )

        # Display biological consequence in a distinct gray box
        st.markdown("#### ⚠️ Biological Consequence")
        st.markdown(
            f"<div style='background-color:#222;color:#ff5555;padding:15px;border-radius:10px;font-weight:bold;'>{bio_consequence}</div>",
            unsafe_allow_html=True
        )
