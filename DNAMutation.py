import streamlit as st
import random
import re

# --- Set Page Config and Custom CSS for Dark Theme ---
st.set_page_config(page_title="🧬 DNA Mutation Lab", page_icon="🔬", layout="wide")

st.markdown(
    """
    <style>
    /* Force full page background to black */
    body {
        background-color: #000000 !important;
        color: #CCCCCC; /* Light grey for general text for better readability on black */
    }
    .stApp {
        background-color: #000000 !important;
    }
    .main {
        background-color: #000000 !important;
    }
    .block-container {
        background-color: #000000 !important;
    }
    header.st-emotion-cache-s1q2x3v { /* Targeting Streamlit's header */
        background-color: #000000 !important;
    }
    div[data-testid="stSidebar"] { /* Targeting Streamlit's sidebar if present */
        background-color: #000000 !important;
    }


    /* Output boxes (DNA, mRNA, Protein sequences) */
    .output-box {
        background-color: #1A1A1A; /* Very dark grey, subtly distinct from black */
        color: #CCCCCC; /* Text inside output boxes is light grey */
        padding: 10px;
        border-radius: 8px;
        font-family: monospace;
        white-space: pre-wrap;
        word-break: break-word;
        border: 1px solid #6272a4; /* Border color remains for structure */
    }

    /* Highlighting for DNA mutation */
    .highlight-red {
        color: #FF7777; /* Slightly brighter red for clarity */
        font-weight: bold;
        background-color: #660000; /* Darker red background */
        padding: 2px 0px;
        border-radius: 3px;
    }

    /* Highlighting for mRNA/Protein differences */
    .highlight-blue {
        color: #99EEFF; /* Brighter cyan/blue */
        font-weight: bold;
        background-color: #004455; /* Darker blue background */
        padding: 2px 0px;
        border-radius: 3px;
    }

    /* Buttons */
    .stButton>button {
        background-color: #6272a4; /* Dracula theme purple */
        color: white;
        border-radius: 8px;
        border: none;
        padding: 10px 20px;
        font-size: 16px;
        transition: background-color 0.3s ease;
    }
    .stButton>button:hover {
        background-color: #8be9fd; /* Dracula theme cyan on hover */
        color: black;
    }

    /* Input fields (text areas, number inputs, select boxes) */
    .stSelectbox, .stNumberInput, .stTextInput, .stTextArea {
        background-color: #222222; /* A slightly lighter dark grey for inputs */
        border-radius: 8px;
        padding: 10px;
    }
    .stTextInput>div>div>input, .stTextArea>div>div>textarea {
        background-color: #222222; /* Consistent background for actual input elements */
        color: #CCCCCC; /* Text inside inputs is light grey */
        border: 1px solid #44475a;
        border-radius: 8px;
        padding: 8px;
    }
    .stSelectbox>div>div>div>span {
        color: #CCCCCC; /* Text inside selectbox is light grey */
    }
    .stNumberInput>div>div>input {
        background-color: #222222; /* Consistent background for number input */
        color: #CCCCCC; /* Text inside number input is light grey */
        border: 1px solid #44475a;
        border-radius: 8px;
        padding: 8px;
    }

    /* Biological Impact Box (the expander content) */
    .biological-impact-box {
        background-color: #1A1A1A; /* Consistent with output boxes */
        color: #CCCCCC; /* Text inside this box is light grey */
        padding: 15px;
        border-radius: 8px;
        border: 1px solid #BD93F9; /* Dracula theme purple border */
        margin-top: 15px;
    }

    /* Specific styling for the summary box (the one showing "Point mutation: Substituted...") */
    /* This targets Streamlit's internal structure for markdown elements to apply specific style */
    div[data-testid="stMarkdownContainer"] div {
        background-color: #1A1A1A !important; /* Force background for summary box, consistent */
        color: #FF7777 !important; /* Red text for summary, slightly brighter */
        border: 1px solid #FF7777 !important; /* Match red highlight border */
        padding: 15px;
        border-radius: 10px;
        font-weight: bold;
    }

    /* Specific style for st.expander header text */
    .st-emotion-cache-1ftrzg7 p { /* Targeting the <p> tag inside the expander header */
        color: #CCCCCC !important; /* Ensure expander title is light grey */
    }
    .st-emotion-cache-1ftrzg7 .st-emotion-cache-j4d57l { /* Targeting the expander icon */
        color: #CCCCCC !important; /* Ensure expander icon is light grey */
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
    start_codon_idx = mrna_seq.find('AUG')
    
    if start_codon_idx == -1 or start_codon_idx % 3 != 0:
        coding_sequence = mrna_seq
    else:
        coding_sequence = mrna_seq[start_codon_idx:]

    for i in range(0, len(coding_sequence) - len(coding_sequence) % 3, 3):
        codon = coding_sequence[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, 'X')
        protein.append(amino_acid)
        if amino_acid == '*':
            break
    return "".join(protein)

def get_codon_at_position(sequence, pos_idx_0based):
    """
    Given a sequence and a 0-based position, returns the 3-base codon
    that contains that position and its 0-based start index.
    Assumes translation starts from index 0 for this specific lookup.
    Returns (codon, start_idx) or (None, None) if position is out of range or
    cannot form a full codon.
    """
    if pos_idx_0based < 0 or pos_idx_0based >= len(sequence):
        return None, None
    
    codon_start_idx = (pos_idx_0based // 3) * 3
    if codon_start_idx + 3 > len(sequence):
        return None, None
    
    codon = sequence[codon_start_idx : codon_start_idx + 3]
    return codon, codon_start_idx


def highlight_mutation_dna(original_seq, mutated_seq, position, mutation_type, length=1):
    """
    Highlight the mutated parts in red for DNA sequences.
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


def highlight_mrna_protein_diff(original_seq, mutated_seq, dna_mutation_pos_idx, mutation_type, mutation_len):
    """
    Highlights differences between original and mutated mRNA/Protein sequences in blue.
    For point mutations, highlights the specific affected codon/amino acid.
    For frameshifts, highlights from the first point of divergence.
    """
    def wrap_blue(text):
        return f"<span class='highlight-blue'>{text}</span>"

    orig_highlighted = list(original_seq)
    mut_highlighted = list(mutated_seq)

    mrna_start_highlight_idx = (dna_mutation_pos_idx // 3) * 3
    protein_start_highlight_idx = dna_mutation_pos_idx // 3

    if mutation_type == "point":
        orig_aa = original_seq[protein_start_highlight_idx] if 0 <= protein_start_highlight_idx < len(original_seq) else None
        mut_aa = mutated_seq[protein_start_highlight_idx] if 0 <= protein_start_highlight_idx < len(mutated_seq) else None

        if orig_aa != mut_aa:
            if 0 <= mrna_start_highlight_idx < len(orig_highlighted) - 2:
                for i in range(3):
                    orig_highlighted[mrna_start_highlight_idx + i] = wrap_blue(orig_highlighted[mrna_start_highlight_idx + i])
                    if mrna_start_highlight_idx + i < len(mut_highlighted):
                        mut_highlighted[mrna_start_highlight_idx + i] = wrap_blue(mut_highlighted[mrna_start_highlight_idx + i])
            
            if 0 <= protein_start_highlight_idx < len(orig_highlighted):
                orig_highlighted_aa = orig_highlighted[protein_start_highlight_idx]
                mut_highlighted_aa = mut_highlighted[protein_start_highlight_idx]
                
                if isinstance(orig_highlighted_aa, str) and orig_highlighted_aa == orig_aa:
                    orig_highlighted[protein_start_highlight_idx] = wrap_blue(orig_aa)
                if isinstance(mut_highlighted_aa, str) and mut_highlighted_aa == mut_aa:
                    mut_highlighted[protein_start_highlight_idx] = wrap_blue(mut_aa)


    elif mutation_type in ["insertion", "deletion"]:
        first_diff_idx = -1
        min_len = min(len(original_seq), len(mutated_seq))
        for i in range(min_len):
            if original_seq[i] != mutated_seq[i]:
                first_diff_idx = i
                break
        
        if first_diff_idx == -1 and len(original_seq) != len(mutated_seq):
            first_diff_idx = min_len 

        if first_diff_idx != -1:
            for i in range(first_diff_idx, len(orig_highlighted)):
                orig_highlighted[i] = wrap_blue(orig_highlighted[i])
            for i in range(first_diff_idx, len(mut_highlighted)):
                mut_highlighted[i] = wrap_blue(mut_highlighted[i])
        
        if first_diff_idx == -1 and len(original_seq) > len(mutated_seq):
            for i in range(len(mutated_seq), len(orig_highlighted)):
                orig_highlighted[i] = wrap_blue(orig_highlighted[i])
        elif first_diff_idx == -1 and len(mutated_seq) > len(original_seq):
            for i in range(len(original_seq), len(mut_highlighted)):
                mut_highlighted[i] = wrap_blue(mut_highlighted[i])


    return "".join(orig_highlighted), "".join(mut_highlighted)


# --- Session State Initialization ---
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

# New: Introduction for beginners
with st.expander("What are DNA Mutations? 🤔 (For Beginners)"):
    st.markdown("""
    **DNA** (Deoxyribonucleic Acid) is the genetic blueprint of life! It's made of a sequence of building blocks called **bases**: Adenine (A), Thymine (T), Guanine (G), and Cytosine (C). Think of it like a long string of letters that contains instructions for building and operating an organism.

    A **mutation** is a change in this DNA sequence. These changes can happen naturally or be caused by external factors. While some mutations have no effect, others can alter the instructions, leading to changes in proteins and potentially affecting an organism's traits or health.

    This lab simulates three main types of mutations:

    * **🎯 Point Mutation (Substitution):** This is when a **single base** in the DNA sequence is replaced by another base.
        * *Example:* `ATGC` becomes `ATTC` (G is replaced by T).

    * **➕ Insertion Mutation:** This is when **one or more extra bases** are added into the DNA sequence.
        * *Example:* `ATGC` becomes `ATG**GTC**C` (GTC is inserted).

    * **➖ Deletion Mutation:** This is when **one or more bases are removed** from the DNA sequence.
        * *Example:* `ATGCGG` becomes `ATGG` (CG is deleted).
    """)
st.markdown("<br>", unsafe_allow_html=True) # Add some space

# DNA sequence input
dna_input_area = st.text_area(
    "🧬 Input DNA Sequence (A, T, G, C only):",
    value=st.session_state.dna_input,
    height=100,
    max_chars=1000,
    help="✍️ Type or paste your DNA sequence here. It will be automatically converted to uppercase."
)
st.session_state.dna_input = dna_input_area.upper()

# Random DNA generator button
if st.button("✨ Generate Random DNA Sequence"):
    random_seq = generate_random_dna(50)
    st.session_state.dna_input = random_seq
    st.session_state.position = 1
    st.rerun()

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
    if mutation_type == "insertion":
        current_max_pos = max_pos + 1
    else:
        current_max_pos = max_pos

    position = st.number_input(
        "📍 Mutation Position (1-based index):",
        min_value=1,
        max_value=current_max_pos,
        value=min(st.session_state.position, current_max_pos),
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
    elif mutation_type == "insertion" and (pos_idx < 0 or pos_idx > len(seq)):
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
        biological_impact_details = ""
        mutation_length_for_highlight = 1

        # --- Perform Mutation ---
        if mutation_type == "point":
            mutated_seq, original_base = apply_point_mutation(seq, pos_idx, st.session_state.new_base)
            if mutated_seq == seq:
                consequence_summary = f"No actual change occurred! 🔄 The base '{original_base}' was already at position {st.session_state.position}. Sequence remains identical. ✅"
                biological_impact_details = "Since no change occurred, there is no biological impact on the sequence."
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

        # --- Biological Impact Classification ---
        if mutation_type == "point" and mutated_seq != seq: # Only classify if an actual change happened
            orig_codon, _ = get_codon_at_position(original_mrna, pos_idx)
            mut_codon, _ = get_codon_at_position(mutated_mrna, pos_idx)

            if orig_codon and mut_codon: # Ensure codons can be extracted
                orig_aa = GENETIC_CODE.get(orig_codon, 'X')
                mut_aa = GENETIC_CODE.get(mut_codon, 'X')

                if orig_aa == mut_aa:
                    biological_impact_details = f"""
                    **Type: Silent Mutation 🤫**
                    * The base change from `{seq[pos_idx]}` to `{st.session_state.new_base}` at DNA position {st.session_state.position} altered the mRNA codon from `{orig_codon}` to `{mut_codon}`.
                    * However, both codons translate to the **same amino acid ({orig_aa})**.
                    * This typically has **no direct biological impact** on the resulting protein due to the redundancy of the genetic code.
                    """
                elif mut_aa == '*':
                    biological_impact_details = f"""
                    **Type: Nonsense Mutation 💀**
                    * The base change from `{seq[pos_idx]}` to `{st.session_state.new_base}` at DNA position {st.session_state.position} altered the mRNA codon from `{orig_codon}` to `{mut_codon}`.
                    * The new codon `{mut_codon}` is a **premature stop codon**.
                    * This leads to a **truncated (shortened) protein**, which is usually **non-functional** and can have severe biological consequences.
                    """
                elif orig_aa != 'X' and mut_aa != 'X': # Exclude 'X' for unknown codons
                    biological_impact_details = f"""
                    **Type: Missense Mutation 🧬**
                    * The base change from `{seq[pos_idx]}` to `{st.session_state.new_base}` at DNA position {st.session_state.position} altered the mRNA codon from `{orig_codon}` to `{mut_codon}`.
                    * This results in a change from amino acid `{orig_aa}` to `{mut_aa}`.
                    * The biological impact can vary:
                        * **Conservative Missense:** If the new amino acid is chemically similar to the original, the impact might be minimal.
                        * **Non-conservative Missense:** If the new amino acid is chemically different, it can significantly alter the protein's shape, function, or stability, potentially causing disease.
                    """
                else: # Fallback for edge cases with invalid codons
                     biological_impact_details = f"""
                    **Type: Point Mutation (Unclassified) ❓**
                    * A base change occurred at position {st.session_state.position}.
                    * Could not precisely classify the amino acid change (e.g., due to an invalid codon).
                    * Typically leads to a change in a single amino acid or premature stop codon.
                    """
            else: # If codon extraction failed
                 biological_impact_details = f"""
                **Type: Point Mutation (Context Issue) ❓**
                * A base change occurred at position {st.session_state.position}.
                * Could not classify the amino acid change precisely, possibly due to the mutation occurring at the very end of the sequence or if the sequence is too short to form a codon.
                """

        elif mutation_type == "insertion":
            if len(st.session_state.insert_seq) % 3 == 0:
                biological_impact_details = f"""
                **Type: In-frame Insertion ➕ (Length: {len(st.session_state.insert_seq)})**
                * The number of inserted bases ({len(st.session_state.insert_seq)}) is a multiple of 3.
                * The genetic "reading frame" of the sequence is **maintained** after the insertion point.
                * This results in the **addition of new amino acids** to the protein sequence.
                * The impact can range from **minor** (if the added amino acids don't disrupt protein folding or active sites) to **significant** (if a crucial region is affected, leading to altered or non-functional proteins).
                """
            else:
                biological_impact_details = f"""
                **Type: Frameshift Insertion 🚨 (Length: {len(st.session_state.insert_seq)})**
                * The number of inserted bases ({len(st.session_state.insert_seq)}) is NOT a multiple of 3.
                * The genetic "reading frame" of the sequence is **shifted** from the point of insertion onwards.
                * This drastically changes all downstream codons, leading to a completely different amino acid sequence.
                * Often results in a **premature stop codon** and a **non-functional protein**, leading to severe consequences.
                """
        elif mutation_type == "deletion":
            if st.session_state.del_length % 3 == 0:
                biological_impact_details = f"""
                **Type: In-frame Deletion ➖ (Length: {st.session_state.del_length})**
                * The number of deleted bases ({st.session_state.del_length}) is a multiple of 3.
                * The genetic "reading frame" of the sequence is **maintained** after the deletion point.
                * This results in the **removal of amino acids** from the protein sequence.
                * The impact can range from **minor** (if the deleted amino acids don't disrupt protein folding or active sites) to **significant** (if a crucial region is affected, leading to altered or non-functional proteins).
                """
            else:
                biological_impact_details = f"""
                **Type: Frameshift Deletion 🚨 (Length: {st.session_state.del_length)})**
                * The number of deleted bases ({st.session_state.del_length}) is NOT a multiple of 3.
                * The genetic "reading frame" of the sequence is **shifted** from the point of deletion onwards.
                * This drastically changes all downstream codons, leading to a completely different amino acid sequence.
                * Often results in a **premature stop codon** and a **non-functional protein**, leading to severe consequences.
                """
        # --- Display Results ---
        orig_highlight_dna, mut_highlight_dna = highlight_mutation_dna(
            seq,
            mutated_seq,
            pos_idx,
            mutation_type,
            length=mutation_length_for_highlight
        )

        # New: Highlight mRNA and Protein
        orig_highlight_mrna, mut_highlight_mrna = highlight_mrna_protein_diff(
            original_mrna,
            mutated_mrna,
            pos_idx, # Pass the DNA mutation index for reference
            mutation_type,
            mutation_length_for_highlight
        )
        orig_highlight_protein, mut_highlight_protein = highlight_mrna_protein_diff(
            original_protein,
            mutated_protein,
            pos_idx, # Pass the DNA mutation index for reference
            mutation_type,
            mutation_length_for_highlight
        )


        st.markdown("<hr style='border: 1px dashed #44475a;'>", unsafe_allow_html=True)

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

        st.markdown("---")

        st.markdown("### 📊 mRNA Sequences (Transcribed from DNA):")
        st.markdown(f"**Original mRNA:** <div class='output-box'>{orig_highlight_mrna}</div>", unsafe_allow_html=True)
        st.markdown(f"**Mutated mRNA:** <div class='output-box'>{mut_highlight_mrna}</div>", unsafe_allow_html=True)


        st.markdown("### 🧪 Protein Sequences (Translated from mRNA):")
        st.markdown(f"**Original Protein:** <div class='output-box'>{orig_highlight_protein}</div>", unsafe_allow_html=True)
        st.markdown(f"**Mutated Protein:** <div class='output-box'>{mut_highlight_protein}</div>", unsafe_allow_html=True)

        st.markdown("<hr style='border: 1px dashed #44475a;'>", unsafe_allow_html=True)

        st.markdown("### 🔬 Summary of Mutation:")
        st.markdown(
            f"""
            <div style='background-color:#1A1A1A;color:#FF7777;padding:15px;border-radius:10px;font-weight:bold; border: 1px solid #FF7777;'>
                {consequence_summary}
            </div>
            """,
            unsafe_allow_html=True
        )

        st.markdown("### 💡 Potential Biological Impact:")
        with st.expander("Click to learn more about the biological impact of this mutation type"):
            st.markdown(
                f"<div class='biological-impact-box'>{biological_impact_details}</div>",
                unsafe_allow_html=True
            )
