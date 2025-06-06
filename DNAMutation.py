import streamlit as st
import random
import re

# Set page config and background styling
st.set_page_config(page_title="DNA Mutation Simulator 🧬", layout="wide")
page_bg_css = """
<style>
body {
    background-color: black;
    color: white;
}
.stMarkdown div {
    color: white;
}
.output-box {
    background-color: #555555;
    padding: 10px;
    border-radius: 8px;
    font-family: monospace;
    white-space: pre-wrap;
    word-wrap: break-word;
}
</style>
"""
st.markdown(page_bg_css, unsafe_allow_html=True)

# DNA validation function
def validate_dna(seq):
    """Validate that the sequence contains only A, T, C, G"""
    return bool(re.fullmatch(r"[ATCGatcg]*", seq))

# Generate a random DNA sequence
def generate_random_dna(length=50):
    return ''.join(random.choice('ATCG') for _ in range(length))

# Mutation functions
def apply_point_mutation(seq, position, new_base):
    original_base = seq[position]
    if original_base.upper() == new_base.upper():
        return None, None
    mutated_seq = seq[:position] + new_base.upper() + seq[position+1:]
    return mutated_seq, original_base.upper()

def apply_insertion_mutation(seq, position, insert_seq):
    mutated_seq = seq[:position] + insert_seq.upper() + seq[position:]
    return mutated_seq, insert_seq.upper()

def apply_deletion_mutation(seq, position, length):
    if position + length > len(seq):
        return None, None
    deleted_seq = seq[position:position+length]
    mutated_seq = seq[:position] + seq[position+length:]
    return mutated_seq, deleted_seq

# Highlight mutation in sequences
def highlight_mutation(original_seq, mutated_seq, position, mutation_type, length=1):
    """
    Returns two HTML strings with mutations highlighted in red.
    """
    def highlight_range(seq, start, end):
        return seq[:start] + \
               f"<span style='color:red;font-weight:bold;'>{seq[start:end]}</span>" + \
               seq[end:]
    
    # Highlight original sequence
    if mutation_type == "point":
        orig_highlight = highlight_range(original_seq, position, position + 1)
        mut_highlight = highlight_range(mutated_seq, position, position + 1)
    elif mutation_type == "insertion":
        orig_highlight = original_seq
        mut_highlight = highlight_range(mutated_seq, position, position + length)
    elif mutation_type == "deletion":
        orig_highlight = highlight_range(original_seq, position, position + length)
        mut_highlight = mutated_seq
    else:
        orig_highlight = original_seq
        mut_highlight = mutated_seq
    
    return orig_highlight, mut_highlight

# UI Header
st.title("DNA Mutation Simulator 🧬🧬")
st.markdown("Welcome! Simulate point mutations, insertions, and deletions in DNA sequences. 🚀")

# Input sequence or generate random
input_option = st.radio("Choose input option:", ["Enter DNA sequence 🧬", "Generate random sequence 🎲"])

if input_option == "Enter DNA sequence 🧬":
    dna_input = st.text_area("Enter DNA sequence (only A, T, G, C):", height=100, max_chars=1000)
else:
    length = st.slider("Select length for random sequence:", min_value=10, max_value=500, value=50)
    dna_input = generate_random_dna(length)
    st.markdown(f"**Random DNA Sequence ({length} bases):**")
    st.markdown(f"<div class='output-box'>{dna_input}</div>", unsafe_allow_html=True)

# Mutation type selection
mutation_type = st.selectbox("Select mutation type:", ["point", "insertion", "deletion"])

# Mutation parameters based on type
if mutation_type == "point":
    position = st.number_input("Mutation position (0-indexed):", min_value=0, max_value=len(dna_input)-1 if dna_input else 0, step=1)
    new_base = st.selectbox("New base:", ['A', 'T', 'C', 'G'])
elif mutation_type == "insertion":
    position = st.number_input("Insertion position (0-indexed):", min_value=0, max_value=len(dna_input), step=1)
    insert_seq = st.text_input("Sequence to insert (A, T, G, C only):")
elif mutation_type == "deletion":
    position = st.number_input("Deletion start position (0-indexed):", min_value=0, max_value=len(dna_input)-1 if dna_input else 0, step=1)
    del_length = st.number_input("Number of bases to delete:", min_value=1, max_value=len(dna_input)-position if dna_input else 1, step=1)

# Button to apply mutation
apply_mutation = st.button("Apply Mutation 🛠️")

# Mutation result output will come in next chunk for better clarity and separation

# Chunk 2

# Function to highlight mutated parts in red
def highlight_mutation(original_seq, mutated_seq, position, mutation_type, length=1):
    def color_text(seq, start, end):
        return (
            seq[:start] +
            f"<span style='color:red; font-weight:bold;'>{seq[start:end]}</span>" +
            seq[end:]
        )

    if mutation_type == "point":
        orig_highlight = color_text(original_seq, position, position + 1)
        mut_highlight = color_text(mutated_seq, position, position + 1)
    elif mutation_type == "insertion":
        orig_highlight = original_seq
        mut_highlight = color_text(mutated_seq, position, position + length)
    elif mutation_type == "deletion":
        orig_highlight = color_text(original_seq, position, position + length)
        mut_highlight = mutated_seq
    else:
        orig_highlight = original_seq
        mut_highlight = mutated_seq

    return orig_highlight, mut_highlight


# Helper functions for mutations
def apply_point_mutation(seq, position, new_base):
    if position < 0 or position >= len(seq):
        return None, None
    original_base = seq[position]
    mutated_seq = seq[:position] + new_base + seq[position+1:]
    return mutated_seq, original_base


def apply_insertion_mutation(seq, position, insert_seq):
    if position < 0 or position > len(seq):
        return None, None
    mutated_seq = seq[:position] + insert_seq + seq[position:]
    return mutated_seq, insert_seq


def apply_deletion_mutation(seq, position, del_length):
    if position < 0 or position + del_length > len(seq):
        return None, None
    deleted = seq[position:position+del_length]
    mutated_seq = seq[:position] + seq[position+del_length:]
    return mutated_seq, deleted


# Mutation consequence description functions
def consequence_point(original_base, new_base):
    return f"Substituted {original_base} with {new_base} at position."


def consequence_insertion(inserted_seq):
    return f"Inserted sequence '{inserted_seq}'."


def consequence_deletion(deleted_seq):
    return f"Deleted sequence '{deleted_seq}'."


# Random DNA generator
def generate_random_dna(length=50):
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(length))


# UI handling for random DNA generation
if st.button("🎲 Generate Random DNA Sequence"):
    rand_dna = generate_random_dna()
    st.session_state.dna_input = rand_dna
    st.experimental_rerun()


# Display mutation result when user clicks "Apply Mutation"
if st.button("🧪 Apply Mutation"):

    seq = st.session_state.dna_input.upper() if "dna_input" in st.session_state else ""

    if not validate_dna(seq):
        st.error("❌ Invalid DNA sequence! Use only A, T, G, and C.")
    else:
        mutation_type = st.session_state.mutation_type if "mutation_type" in st.session_state else ""
        position = st.session_state.position if "position" in st.session_state else None

        if mutation_type == "point":
            new_base = st.session_state.new_base if "new_base" in st.session_state else None
            mutated_seq, original_base = apply_point_mutation(seq, position, new_base)
            consequence = consequence_point(original_base, new_base)
            length = 1

        elif mutation_type == "insertion":
            insert_seq = st.session_state.insert_seq if "insert_seq" in st.session_state else None
            mutated_seq, inserted = apply_insertion_mutation(seq, position, insert_seq)
            consequence = consequence_insertion(inserted)
            length = len(insert_seq)

        elif mutation_type == "deletion":
            del_length = st.session_state.del_length if "del_length" in st.session_state else None
            mutated_seq, deleted = apply_deletion_mutation(seq, position, del_length)
            consequence = consequence_deletion(deleted)
            length = del_length

        else:
            mutated_seq = None
            consequence = ""
            length = 1

        if mutated_seq:
            orig_highlight, mut_highlight = highlight_mutation(
                seq, mutated_seq, position, mutation_type, length=length
            )

            st.markdown("#### 🧬 Original Sequence")
            st.markdown(f"<div style='background-color:#555; padding:10px; border-radius:8px;'>{orig_highlight}</div>", unsafe_allow_html=True)

            st.markdown("#### 🧬 Mutated Sequence")
            st.markdown(f"<div style='background-color:#555; padding:10px; border-radius:8px;'>{mut_highlight}</div>", unsafe_allow_html=True)

            st.markdown(f"#### ⚠️ Biological Consequence")
            st.markdown(f"<div style='background-color:#333; color:white; padding:10px; border-radius:8px;'>{consequence}</div>", unsafe_allow_html=True)

        else:
            st.error("⚠️ Mutation failed: Check your position and input values.")
# Chunk 3

# Function to validate DNA sequence - only A, T, G, C allowed
def validate_dna(seq):
    return all(base in "ATGC" for base in seq)


# Streamlit input widgets for mutation parameters

st.markdown("### ⚙️ Mutation Settings")

# Mutation type selector with emojis
mutation_type = st.selectbox(
    "Select Mutation Type 🧬",
    options=["point", "insertion", "deletion"],
    index=0,
    help="Choose the type of mutation to apply"
)

# Save mutation type in session state for usage elsewhere
st.session_state.mutation_type = mutation_type

# Input DNA sequence box with gray background
dna_input = st.text_area(
    "Enter DNA Sequence 🔡",
    value=st.session_state.get("dna_input", ""),
    height=100,
    max_chars=1000,
    help="Input your DNA sequence here (A, T, G, C only).",
    key="dna_input"
)

# Position input (0-indexed) with validation
position = st.number_input(
    "Mutation Position (0-indexed) 🔢",
    min_value=0,
    max_value=max(len(dna_input) - 1, 0),
    value=0,
    step=1,
    help="Position in the sequence where mutation happens.",
    key="position"
)
st.session_state.position = position

# Conditional inputs based on mutation type
if mutation_type == "point":
    new_base = st.selectbox(
        "New Base 🔄",
        options=["A", "T", "G", "C"],
        index=0,
        help="Choose the base to substitute at the position.",
        key="new_base"
    )

elif mutation_type == "insertion":
    insert_seq = st.text_input(
        "Sequence to Insert ➕",
        value=st.session_state.get("insert_seq", ""),
        max_chars=100,
        help="Enter sequence to insert at the position.",
        key="insert_seq"
    )

elif mutation_type == "deletion":
    del_length = st.number_input(
        "Number of Bases to Delete ➖",
        min_value=1,
        max_value=max(len(dna_input) - position, 1),
        value=1,
        step=1,
        help="How many bases to delete starting from position.",
        key="del_length"
    )
# Chunk 4

# Function to apply point mutation
def apply_point_mutation(seq, pos, new_base):
    if pos >= len(seq) or new_base not in "ATGC":
        return None, None
    original_base = seq[pos]
    mutated_seq = seq[:pos] + new_base + seq[pos + 1:]
    return mutated_seq, original_base

# Function to apply insertion mutation
def apply_insertion_mutation(seq, pos, insert_seq):
    if pos > len(seq):
        return None, None
    if not all(base in "ATGC" for base in insert_seq):
        return None, None
    mutated_seq = seq[:pos] + insert_seq + seq[pos:]
    return mutated_seq, insert_seq

# Function to apply deletion mutation
def apply_deletion_mutation(seq, pos, length):
    if pos + length > len(seq):
        return None, None
    deleted_seq = seq[pos:pos + length]
    mutated_seq = seq[:pos] + seq[pos + length:]
    return mutated_seq, deleted_seq

# Function to highlight mutations in red
def highlight_mutation(original, mutated, pos, mutation_type, length=1):
    def colorize(seq, start, end):
        return (
            seq[:start] +
            f"<span style='color:red; font-weight:bold;'>{seq[start:end]}</span>" +
            seq[end:]
        )

    orig_highlight = original
    mut_highlight = mutated
    if mutation_type == "point":
        orig_highlight = colorize(original, pos, pos + 1)
        mut_highlight = colorize(mutated, pos, pos + 1)
    elif mutation_type == "insertion":
        mut_highlight = colorize(mutated, pos, pos + length)
    elif mutation_type == "deletion":
        orig_highlight = colorize(original, pos, pos + length)
    return orig_highlight, mut_highlight


# Process mutation and display results

if dna_input:
    seq = dna_input.upper()
    st.markdown("<h3 style='color:white;'>🔬 Mutation Result</h3>", unsafe_allow_html=True)

    if not validate_dna(seq):
        st.error("❌ Invalid DNA sequence. Use only A, T, G, and C.")
    else:
        mutated_seq = None
        consequence = ""

        if mutation_type == "point":
            mutated_seq, original_base = apply_point_mutation(seq, position, new_base)
            if mutated_seq:
                consequence = f"Substituted {original_base} with {new_base} at position {position}."
            else:
                st.error("⚠️ Invalid point mutation parameters.")
        elif mutation_type == "insertion":
            mutated_seq, inserted = apply_insertion_mutation(seq, position, insert_seq)
            if mutated_seq:
                consequence = f"Inserted {inserted} at position {position}."
            else:
                st.error("⚠️ Invalid insertion parameters.")
        elif mutation_type == "deletion":
            mutated_seq, deleted = apply_deletion_mutation(seq, position, del_length)
            if mutated_seq:
                consequence = f"Deleted {deleted} from position {position} to {position + del_length - 1}."
            else:
                st.error("⚠️ Invalid deletion parameters.")

        if mutated_seq:
            orig_highlight, mut_highlight = highlight_mutation(
                seq, mutated_seq, position, mutation_type,
                length=(len(insert_seq) if mutation_type == "insertion" else del_length if mutation_type == "deletion" else 1)
            )

            st.markdown("#### 🧬 Original Sequence")
            st.markdown(f"<div style='background-color:#444;padding:10px;border-radius:8px;'>{orig_highlight}</div>", unsafe_allow_html=True)

            st.markdown("#### 🧬 Mutated Sequence")
            st.markdown(f"<div style='background-color:#444;padding:10px;border-radius:8px;'>{mut_highlight}</div>", unsafe_allow_html=True)

            st.markdown(f"#### ⚠️ Biological Consequence")
            st.markdown(f"<div style='background-color:#333;color:white;padding:10px;border-radius:8px;'>{consequence}</div>", unsafe_allow_html=True)
