import streamlit as st
import random

# --- Styling for Dark Theme ---
st.set_page_config(page_title="DNA Mutation Simulator", layout="wide")
st.markdown(
    """
    <style>
    body {
        background-color: #000;
    }
    .main {
        background-color: #000;
        color: white;
    }
    div.stTextInput > label, div.stNumberInput > label, div.stSelectbox > label {
        color: white;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# --- Utility Functions ---
def validate_dna(seq):
    return all(base in "ATGC" for base in seq.upper())

def random_dna(length=50):
    return ''.join(random.choices("ATGC", k=length))

def apply_point_mutation(seq, pos, new_base):
    if pos < 0 or pos >= len(seq) or new_base not in "ATGC":
        return None, None
    original_base = seq[pos]
    mutated_seq = seq[:pos] + new_base + seq[pos+1:]
    return mutated_seq, original_base

def apply_insertion_mutation(seq, pos, insert_seq):
    if pos < 0 or pos > len(seq) or not insert_seq or not validate_dna(insert_seq):
        return None, None
    mutated_seq = seq[:pos] + insert_seq + seq[pos:]
    return mutated_seq, insert_seq

def apply_deletion_mutation(seq, pos, del_len):
    if pos < 0 or pos + del_len > len(seq):
        return None, None
    deleted = seq[pos:pos+del_len]
    mutated_seq = seq[:pos] + seq[pos+del_len:]
    return mutated_seq, deleted

def highlight_mutation(original, mutated, pos, mutation_type, length=1):
    if mutation_type == "point":
        orig = (
            original[:pos] + f"<span style='color:red;'>{original[pos]}</span>" + original[pos+1:]
        )
        mut = (
            mutated[:pos] + f"<span style='color:red;'>{mutated[pos]}</span>" + mutated[pos+1:]
        )
    elif mutation_type == "insertion":
        orig = (
            original[:pos] + f"<span style='color:red;'></span>" + original[pos:]
        )
        mut = (
            mutated[:pos] + f"<span style='color:red;'>{mutated[pos:pos+length]}</span>" + mutated[pos+length:]
        )
    elif mutation_type == "deletion":
        orig = (
            original[:pos] + f"<span style='color:red;'>{original[pos:pos+length]}</span>" + original[pos+length:]
        )
        mut = (
            mutated[:pos] + f"<span style='color:red;'></span>" + mutated[pos:]
        )
    else:
        orig, mut = original, mutated
    return orig, mut
# Sidebar – Mutation Options 🛠️
st.sidebar.markdown("## 🛠️ Mutation Options")

mutation_type = st.sidebar.selectbox(
    "Select Mutation Type",
    options=["point", "insertion", "deletion"]
)

# Set defaults to avoid key errors
if "dna_input" not in st.session_state:
    st.session_state.dna_input = ""
if "position" not in st.session_state:
    st.session_state.position = 1
if "new_base" not in st.session_state:
    st.session_state.new_base = "A"
if "insert_seq" not in st.session_state:
    st.session_state.insert_seq = ""
if "del_length" not in st.session_state:
    st.session_state.del_length = 1

# DNA Input from User
dna_input = st.text_area("🧬 Enter DNA Sequence (A, T, G, C only)", height=150)
if dna_input:
    st.session_state.dna_input = dna_input.upper()

# Controls based on mutation type
if mutation_type == "point":
    position = st.sidebar.number_input("Position for Point Mutation", min_value=1, value=st.session_state.position)
    new_base = st.sidebar.selectbox("New Base", options=["A", "T", "G", "C"], index=0)
    st.session_state.position = position
    st.session_state.new_base = new_base

elif mutation_type == "insertion":
    position = st.sidebar.number_input("Position to Insert", min_value=1, value=st.session_state.position)
    insert_seq = st.sidebar.text_input("Sequence to Insert", value=st.session_state.insert_seq)
    st.session_state.position = position
    st.session_state.insert_seq = insert_seq.upper()

elif mutation_type == "deletion":
    position = st.sidebar.number_input("Position to Delete From", min_value=1, value=st.session_state.position)
    del_length = st.sidebar.number_input("Number of Bases to Delete", min_value=1, value=st.session_state.del_length)
    st.session_state.position = position
    st.session_state.del_length = del_length

# 🎲 Random Generation Button
if st.sidebar.button("🎲 Generate Random DNA"):
    random_seq = generate_random_dna(50)
    st.session_state.dna_input = random_seq
    st.experimental_rerun()
# ✅ Once we have a DNA input, process the mutation
if st.session_state.dna_input:
    seq = st.session_state.dna_input.upper()
    st.markdown("<h3 style='color:white;'>🔬 Mutation Result</h3>", unsafe_allow_html=True)

    if not validate_dna(seq):
        st.error("🚫 Invalid DNA sequence. Please use only A, T, G, and C.")
    else:
        mutated_seq = seq
        consequence = ""
        if mutation_type == "point":
            mutated_seq, original_base = apply_point_mutation(seq, st.session_state.position, st.session_state.new_base)
            consequence = f"🔁 Substituted **{original_base}** ➡️ **{st.session_state.new_base}** at position {st.session_state.position}."

        elif mutation_type == "insertion":
            mutated_seq, inserted = apply_insertion_mutation(seq, st.session_state.position, st.session_state.insert_seq)
            consequence = f"➕ Inserted **{inserted}** at position {st.session_state.position}."

        elif mutation_type == "deletion":
            mutated_seq, deleted = apply_deletion_mutation(seq, st.session_state.position, st.session_state.del_length)
            consequence = f"➖ Deleted **{deleted}** from position {st.session_state.position} to {st.session_state.position + st.session_state.del_length - 1}."

        if mutated_seq:
            orig_highlight, mut_highlight = highlight_mutation(
                seq, mutated_seq,
                st.session_state.position,
                mutation_type,
                length=(len(st.session_state.insert_seq) if mutation_type == "insertion" else st.session_state.del_length)
            )

            # 🧬 Display original and mutated sequences
            st.markdown("#### 🧬 Original Sequence")
            st.markdown(
                f"<div style='background-color:#444;padding:10px;border-radius:8px;color:white;'>{orig_highlight}</div>",
                unsafe_allow_html=True
            )

            st.markdown("#### 🧬 Mutated Sequence")
            st.markdown(
                f"<div style='background-color:#444;padding:10px;border-radius:8px;color:white;'>{mut_highlight}</div>",
                unsafe_allow_html=True
            )

            # ⚠️ Display consequence
            st.markdown("#### ⚠️ Biological Consequence")
            st.markdown(
                f"<div style='background-color:#333;color:white;padding:10px;border-radius:8px;'>{consequence}</div>",
                unsafe_allow_html=True
            )
        else:
            st.error("⚠️ Mutation failed due to invalid parameters.")
# ✅ --- Mutation Functions ---
def apply_point_mutation(seq, position, new_base):
    if 0 <= position < len(seq):
        original_base = seq[position]
        mutated_seq = seq[:position] + new_base + seq[position + 1:]
        return mutated_seq, original_base
    return None, ""

def apply_insertion_mutation(seq, position, insert_seq):
    if 0 <= position <= len(seq):
        mutated_seq = seq[:position] + insert_seq + seq[position:]
        return mutated_seq, insert_seq
    return None, ""

def apply_deletion_mutation(seq, position, length):
    if 0 <= position < len(seq) and position + length <= len(seq):
        deleted_seq = seq[position:position + length]
        mutated_seq = seq[:position] + seq[position + length:]
        return mutated_seq, deleted_seq
    return None, ""

# ✅ --- Highlighting Function ---
def highlight_mutation(original, mutated, pos, m_type, length=1):
    def wrap(seq, start, length):
        return seq[:start] + f"<span style='color:red;font-weight:bold'>" + seq[start:start+length] + "</span>" + seq[start+length:]

    if m_type == "point":
        return wrap(original, pos, 1), wrap(mutated, pos, 1)
    elif m_type == "insertion":
        return original[:pos] + "<span style='color:yellow'>&#x25B6;</span>" + original[pos:], wrap(mutated, pos, length)
    elif m_type == "deletion":
        return wrap(original, pos, length), mutated[:pos] + "<span style='color:yellow'>&#x25C0;</span>" + mutated[pos:]
    else:
        return original, mutated

# ✅ --- DNA Validation Function ---
def validate_dna(seq):
    return all(base in "ATGC" for base in seq.upper())

# ✅ --- Streamlit Footer (Optional Aesthetic) ---
st.markdown("""
<br><hr>
<center style='color:white;'>🧬 Built with ❤️ using Streamlit | Final Project © 2025</center>
""", unsafe_allow_html=True)
