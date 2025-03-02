import streamlit as st
import os
import time
import requests
import json
from PhrankCalculator import PhrankCalculator  # Your existing class

# Static file paths
OBO_PATH = "data/hp.obo"  # Update with your actual path
HPOA_PATH = "data/phenotype.hpoa"  # Update with your actual path

# Llama integration for symptom recommendations
def get_llama_recommendations(patient_terms, ontology):
    """
    Get symptom recommendations from locally running Llama model
    """
    try:
        # Convert HPO terms to their human-readable names for better context
        phenotype_names = []
        for term in patient_terms:
            if term in ontology:
                phenotype_names.append(f"{term} ({ontology[term].name})")
            else:
                phenotype_names.append(term)
                
        # Prepare prompt for Llama
        prompt = f"""
        Based on the following patient phenotypes:
        {', '.join(phenotype_names)}
        
        Suggest 8-10 additional symptoms or phenotypic features that might be present in this patient.
        Make your suggestions specific and clinically relevant.
        For each suggestion, provide a brief explanation of why it might be relevant.
        Format each suggestion as: "Symptom: [brief explanation]"
        """
        
        # Call local Llama API
        response = requests.post(
            "http://localhost:11434/api/generate",
            json={
                "model": "llama3",
                "prompt": prompt,
                "stream": False
            }
        )
        
        if response.status_code == 200:
            result = response.json()
            suggestions = result.get('response', '')
            
            # Parse the suggestions from Llama
            raw_suggestions = [s.strip() for s in suggestions.split('\n') if s.strip() and ':' in s]
            
            # Map suggestions to HPO terms using similarity matching
            # In a real implementation, you'd use a more sophisticated mapping algorithm
            mapped_suggestions = map_to_hpo_terms(raw_suggestions, ontology)
            
            return mapped_suggestions
        else:
            st.error(f"Error connecting to Llama API: {response.status_code}")
            return []
    except Exception as e:
        st.error(f"Error getting recommendations: {str(e)}")
        return []

def map_to_hpo_terms(raw_suggestions, ontology):
    """
    Map natural language symptom descriptions to HPO terms
    This is a simplified approach - in a real app, you'd use a more sophisticated mapping
    """
    # This is a simplified mapping - in reality you would:
    # 1. Use a more sophisticated text similarity or NLP approach
    # 2. Or have a dedicated API/service for this mapping
    
    # Simulate mapping with a simplified approach
    mapped_results = []
    
    # Get a list of all HPO terms and their names for matching
    all_terms = {term_id: term_obj.name.lower() for term_id, term_obj in ontology.items()}
    
    for suggestion in raw_suggestions:
        # Extract the symptom part (before the explanation)
        if ':' in suggestion:
            symptom_text = suggestion.split(':', 1)[0].strip().lower()
            explanation = suggestion.split(':', 1)[1].strip()
            
            # Find potential matches in HPO
            matches = []
            for term_id, term_name in all_terms.items():
                # Simple matching score based on word overlap
                suggestion_words = set(symptom_text.split())
                term_words = set(term_name.split())
                overlap = len(suggestion_words.intersection(term_words))
                
                if overlap > 0:
                    matches.append((term_id, overlap, term_name))
            
            # Sort by overlap score and take the best match
            if matches:
                matches.sort(key=lambda x: x[1], reverse=True)
                best_match = matches[0]
                mapped_results.append({
                    "term_id": best_match[0],
                    "term_name": best_match[2],
                    "original_suggestion": symptom_text,
                    "explanation": explanation,
                    "confidence": min(best_match[1] / len(symptom_text.split()), 0.95)  # Simple confidence score
                })
            
    # Return only the highest confidence and unique matches
    seen_terms = set()
    final_results = []
    
    for match in sorted(mapped_results, key=lambda x: x["confidence"], reverse=True):
        if match["term_id"] not in seen_terms and match["confidence"] > 0.3:
            seen_terms.add(match["term_id"])
            final_results.append(match)
            if len(final_results) >= 10:
                break
                
    return final_results

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
if 'llama_recommendations' not in st.session_state:
    st.session_state.llama_recommendations = []
if 'similarity_score' not in st.session_state:
    st.session_state.similarity_score = 0
if 'show_results' not in st.session_state:
    st.session_state.show_results = False

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
                    # Placeholder for actual ClinPhen implementation
                    time.sleep(2)
                    generated_terms = {
                        'HP:0000726', 'HP:0000750', 'HP:0001250', 
                        'HP:0001263', 'HP:0002172', 'HP:0012447'
                    }
                    st.session_state.patient_terms = generated_terms
                    st.session_state.hpo_generated = True
                    st.rerun()

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
                
                # Recalculate similarity score
                if st.session_state.patient_terms:
                    with st.spinner("Updating disease matches..."):
                        diseases = st.session_state.phrank.rank_diseases(
                            st.session_state.patient_terms, 
                            top_n=10
                        )
                        if diseases:
                            st.session_state.similarity_score = diseases[0][0]
                            if st.session_state.similarity_score >= 60:
                                st.session_state.show_results = True
                
                # Get new recommendations based on updated terms
                if st.session_state.phrank and st.session_state.patient_terms:
                    with st.spinner("Generating new recommendations..."):
                        st.session_state.llama_recommendations = get_llama_recommendations(
                            st.session_state.patient_terms,
                            st.session_state.phrank.ontology
                        )
                
                st.rerun()

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
            st.session_state.similarity_score = 0
            st.session_state.show_results = False
            st.rerun()
    
    # Llama Recommendations
    if st.session_state.phrank and st.session_state.patient_terms and not st.session_state.show_results:
        st.subheader("Recommended Symptoms")
        
        if not st.session_state.llama_recommendations:
            if st.button("Get Recommendations"):
                with st.spinner("Generating recommendations..."):
                    st.session_state.llama_recommendations = get_llama_recommendations(
                        st.session_state.patient_terms,
                        st.session_state.phrank.ontology
                    )
                    st.rerun()
        else:
            st.markdown("""
            Based on current phenotypes, consider these additional symptoms. 
            Select relevant ones to add to the patient profile:
            """)
            
            # Create a container for recommendations with scrolling
            with st.container():
                recommendations_to_add = []
                
                for idx, rec in enumerate(st.session_state.llama_recommendations):
                    term_id = rec["term_id"]
                    
                    # Skip if already in patient terms
                    if term_id in st.session_state.patient_terms:
                        continue
                        
                    with st.container():
                        col1, col2 = st.columns([5, 1])
                        
                        with col1:
                            st.markdown(f"""
                            **{rec['term_name']}** ({term_id})  
                            *Original suggestion: {rec['original_suggestion']}*  
                            {rec['explanation']}
                            """)
                        
                        with col2:
                            if st.checkbox("Select", key=f"rec_{idx}"):
                                recommendations_to_add.append(term_id)
                        
                        st.divider()
                
                # Add selected recommendations
                if recommendations_to_add:
                    if st.button("Add Selected Terms"):
                        st.session_state.patient_terms.update(recommendations_to_add)
                        
                        # Recalculate similarity score
                        with st.spinner("Updating disease matches..."):
                            diseases = st.session_state.phrank.rank_diseases(
                                st.session_state.patient_terms, 
                                top_n=10
                            )
                            if diseases:
                                st.session_state.similarity_score = diseases[0][0]
                                if st.session_state.similarity_score >= 60:
                                    st.session_state.show_results = True
                        
                        # Clear recommendations to get fresh ones
                        st.session_state.llama_recommendations = []
                        st.rerun()
            
            if st.button("Refresh Recommendations"):
                st.session_state.llama_recommendations = []
                st.rerun()
        
        # Show current similarity score if available
        if st.session_state.similarity_score > 0:
            progress_color = "red" if st.session_state.similarity_score < 40 else "orange" if st.session_state.similarity_score < 60 else "green"
            st.markdown(f"""
            <div style="margin-top: 20px;">
                <p>Current similarity score: <span style="color:{progress_color}; font-weight:bold;">{st.session_state.similarity_score:.1f}/100</span></p>
                <div style="width:100%; background-color:#ddd; border-radius:5px;">
                    <div style="width:{min(st.session_state.similarity_score, 100)}%; background-color:{progress_color}; height:20px; border-radius:5px;"></div>
                </div>
                <p style="font-size:0.8em; margin-top:5px;">Need to reach 60 to view matches</p>
            </div>
            """, unsafe_allow_html=True)

# Analysis Section
if st.session_state.phrank and st.session_state.patient_terms and st.session_state.show_results:
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
    score_threshold = 60.0
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
    if suggestions:
        st.subheader("Additional Critical Terms")
        st.markdown("These terms are strongly associated with top matches:")
        
        rec_cols = st.columns(3)
        for idx, (term, ic) in enumerate(suggestions[:9]):  # Show top 9 suggestions
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
                        st.rerun()

# Visualization
with vis_col:
    st.header("Ontology Visualization")
    viz_type = st.radio("Select Visualization:", 
                      ["Disease Network", "HPO Structure"])
    
    with st.spinner("Generating visualization..."):
        if st.session_state.phrank and st.session_state.patient_terms:
            if viz_type == "Disease Network":
                img_path = st.session_state.phrank.visualize_disease_similarity_network(
                    st.session_state.patient_terms
                )
            else:
                img_path = st.session_state.phrank.visualize_dag()
            
            st.image(img_path, use_column_width=True)
        else:
            st.info("Add patient phenotypes to view visualizations")

# Initial State Messaging
if not st.session_state.phrank:
    st.info("Please ensure HPO data files are correctly loaded to begin analysis")
elif not st.session_state.patient_terms:
    st.info("Please input patient phenotypes using clinical notes or manual entry above")