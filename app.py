import streamlit as st
import os
import time
from PhrankCalculator import PhrankCalculator  # Your existing class


# Static file paths
OBO_PATH = "data/hp.obo"  # Update with your actual path
HPOA_PATH = "data/phenotype.hpoa"  # Update with your actual path

# Placeholder ClinPhen API integration (replace with actual implementation)
def clinphen_to_hpo(clinical_notes):
    """Simulate ClinPhen API call to convert clinical notes to HPO terms"""
    # In reality, this would call the ClinPhen API/package
    return {
        'HP:0000726', 'HP:0000750', 'HP:0001250', 
        'HP:0001263', 'HP:0002172', 'HP:0012447'
    }

# Streamlit App Configuration
st.set_page_config(
    page_title="HPO Disease Matcher",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="auto"
)

# Session State Initialization
if 'phrank' not in st.session_state:
    try:
        st.session_state.phrank = PhrankCalculator(
            obo_file=OBO_PATH,
            hpoa_file=HPOA_PATH
        )
    except Exception as e:
        st.error(f"Error loading HPO data: {str(e)}")
        st.session_state.phrank = None

if 'patient_terms' not in st.session_state:
    st.session_state.patient_terms = set()
if 'clinical_notes' not in st.session_state:
    st.session_state.clinical_notes = ""
if 'hpo_generated' not in st.session_state:
    st.session_state.hpo_generated = False

# # Sidebar - Data Configuration
# with st.sidebar:
#     st.header("Data Configuration")
#     obo_file = st.file_uploader("Upload HPO OBO File", type=['obo'])
#     hpoa_file = st.file_uploader("Upload Phenotype HPOA File", type=['hpoa'])
    
#     if obo_file and hpoa_file:
#         with st.spinner("Loading HPO data..."):
#             # Save uploaded files
#             obo_path = os.path.join("data", "hp.obo")
#             hpoa_path = os.path.join("data", "phenotype.hpoa")
            
#             with open(obo_path, "wb") as f:
#                 f.write(obo_file.getbuffer())
#             with open(hpoa_path, "wb") as f:
#                 f.write(hpoa_file.getbuffer())
            
#             # Initialize Phrank calculator
#             try:
#                 st.session_state.phrank = PhrankCalculator(
#                     obo_file=obo_path,
#                     hpoa_file=hpoa_path
#                 )
#                 st.success("HPO data loaded successfully!")
#             except Exception as e:
#                 st.error(f"Error loading HPO data: {str(e)}")

# Main Application Interface
st.title("Phenotype-Driven Disease Matching")
st.markdown("""
This application helps match patient phenotypes to known genetic disorders using the Human Phenotype Ontology (HPO).
""")

# Patient Input Section
input_col, vis_col = st.columns([2, 1])

with input_col:
    st.header("Patient Phenotype Input")
    
    # Clinical Notes Input (for HPO generation)
    if not st.session_state.patient_terms and not st.session_state.hpo_generated:
        with st.expander("Clinical Notes to HPO Conversion", expanded=True):
            st.markdown("""
            **Don't have HPO terms?**  
            Paste clinical notes below to generate HPO terms automatically.
            """)
            
            clinical_notes = st.text_area(
                "Clinical Notes (English only):",
                height=200,
                key="clinical_notes"
            )
            
        if st.button("Generate HPO Terms"):
            if clinical_notes.strip():
                with st.spinner("Analyzing clinical notes..."):
                    time.sleep(2)
                    generated_terms = clinphen_to_hpo(clinical_notes)
                    st.session_state.patient_terms = generated_terms
                    st.session_state.hpo_generated = True
                    st.rerun()  # Changed from st.experimental_user()

    # Direct HPO Input
    with st.expander("Manual HPO Term Entry", expanded=bool(st.session_state.patient_terms)):
        hpo_input = st.text_input(
            "Add HPO terms (comma-separated, e.g., HP:0001234,HP:0005678):"
        )
        
        if st.button("Add Terms"):
            new_terms = {t.strip() for t in hpo_input.split(',') if t.strip()}
            valid_terms = set()
            
            if st.session_state.phrank:
                for term in new_terms:
                    if term in st.session_state.phrank.ontology:
                        valid_terms.add(term)
                    else:
                        st.warning(f"Invalid HPO term skipped: {term}")
                
                st.session_state.patient_terms.update(valid_terms)
                st.rerun()  # Changed from st.experimental_user()

    # Current Terms Display
    if st.session_state.patient_terms:
        st.subheader("Current Phenotype Profile")
        term_cols = st.columns(3)
        
        for idx, term in enumerate(sorted(st.session_state.patient_terms)):
            with term_cols[idx % 3]:
                term_name = st.session_state.phrank.ontology[term].name if term in st.session_state.phrank.ontology else "Unknown term"
                st.markdown(f"""
                **{term}**  
                *{term_name}*  
                """)
        
        if st.button("Clear All Terms"):
            st.session_state.patient_terms = set()
            st.experimental_user()

# Analysis Section
if st.session_state.phrank and st.session_state.patient_terms:
    st.header("Analysis Results")
    
    # Run Analysis
    with st.spinner("Analyzing phenotype matches..."):
        diseases = st.session_state.phrank.rank_diseases(
            st.session_state.patient_terms, 
            top_n=10
        )
        suggestions = st.session_state.phrank.suggest_terms(
            diseases, 
            st.session_state.patient_terms
        )
    
    # Display Results
    score_threshold = 15.0  # Adjust based on validation
    top_score = diseases[0][0] if diseases else 0
    
    if top_score >= score_threshold:
        st.success(f"Confident Match Found (Score: {top_score:.1f})")
    else:
        st.warning(f"Borderline Match (Score: {top_score:.1f}) - Consider Adding More Terms")
    
    # Top Matches Table
    st.subheader("Top Disease Matches")
    for idx, (score, disease_id, name) in enumerate(diseases[:5]):
        with st.container():
            cols = st.columns([1, 4, 2])
            cols[0].markdown(f"**#{idx+1}**")
            cols[1].markdown(f"**{name}**  \n`{disease_id}`")
            cols[2].markdown(f"**Score:** {score:.2f}")
            st.divider()
    
    # Term Suggestions
    if top_score < score_threshold and suggestions:
        st.subheader("Recommended Additional Terms")
        st.markdown("Consider adding these high-impact terms from top matches:")
        
        rec_cols = st.columns(3)
        for idx, (term, ic) in enumerate(suggestions):
            with rec_cols[idx % 3]:
                term_name = st.session_state.phrank.ontology[term].name
                with st.container():
                    st.markdown(f"""
                    **{term}**  
                    *{term_name}*  
                    Information Content: `{ic:.2f}`
                    """)
                    if st.button(f"Add {term}", key=f"add_{term}"):
                        st.session_state.patient_terms.add(term)
                        st.experimental_user()
    
    # Visualization
    with vis_col:
        st.header("Ontology Visualization")
        viz_type = st.radio("Select Visualization:", 
                          ["Disease Network", "HPO Structure"])
        
        with st.spinner("Generating visualization..."):
            if viz_type == "Disease Network":
                img_path = st.session_state.phrank.visualize_disease_similarity_network(
                    st.session_state.patient_terms
                )
            else:
                img_path = st.session_state.phrank.visualize_dag()
            
            st.image(img_path, use_column_width=True)

# Initial State Messaging
elif not st.session_state.phrank:
    st.info("Please load HPO data files in the sidebar to begin analysis")

elif not st.session_state.patient_terms:
    st.info("Please input patient phenotypes using clinical notes or manual entry above")