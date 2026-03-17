"""
RAPTOR Dashboard - Hybrid Recommender (Module 4)

Pipeline recommendation using BOTH rule-based and ML approaches.
Provides comprehensive recommendations with comparison and agreement analysis.

Author: Ayeh Bolouki
Version: 2.2.0
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

# Import RAPTOR recommender modules
try:
    from raptor.recommender import (
        PipelineRecommender,
        Recommendation,
        PipelineInfo,
        recommend_pipeline,
        get_pipeline_info,
        list_pipelines
    )
    RULE_BASED_AVAILABLE = True
except ImportError as e:
    RULE_BASED_AVAILABLE = False
    rule_import_error = str(e)

try:
    from raptor.ml_recommender import (
        MLPipelineRecommender,
        MLRecommendation,
        quick_ml_recommend
    )
    ML_RECOMMENDER_AVAILABLE = True
except ImportError as e:
    ML_RECOMMENDER_AVAILABLE = False
    ml_import_error = str(e)

try:
    from raptor.profiler import DataProfile
    PROFILER_AVAILABLE = True
except ImportError:
    PROFILER_AVAILABLE = False

# Page config
st.set_page_config(
    page_title="RAPTOR - Hybrid Recommender",
    page_icon="🤖",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# Main content
st.title("Pipeline Recommender — Module 4")
st.caption("Pipeline recommendation using rule-based and ML approaches")

# Help section
with st.expander("ℹ️ How to use this module"):
    st.markdown("""
    **Purpose:** Get the best DE pipeline recommendation for your data using TWO complementary approaches.
    
    ## **Dual Recommendation System**
    
    ### **1. Rule-Based Recommender** Always Available
    **Approach:** Literature-based heuristic rules
    
    **Decision Factors:**
    - Sample size (most critical)
    - BCV (Biological Coefficient of Variation)
    - Outlier presence
    - Low count proportion
    - Design complexity
    
    **Advantages:**
    - Always works (no training needed)
    - Fast and interpretable
    - Based on benchmarking research
    - Explains reasoning clearly
    
    ### **2. ML-Based Recommender** If Model Available
    **Approach:** Random Forest trained on synthetic benchmarks
    
    **Training Process:**
    1. Generate diverse synthetic RNA-seq datasets
    2. Run DESeq2, edgeR, limma-voom, Wilcoxon on each
    3. Evaluate against ground truth
    4. Train classifier on data features
    
    **Advantages:**
    - More accurate (if well-trained)
    - Learns from actual performance
    - Provides confidence scores
    - Shows feature importances
    
    ### **3. Comparison Mode** 🔍 Best of Both Worlds
    **When both available:**
    - See both recommendations side-by-side
    - Check if they agree (validation!)
    - Understand different reasoning
    - Make informed decision
    
    **Agreement Analysis:**
    - **Both Agree:** High confidence - use recommended pipeline
    - **Disagree:** Investigate further, consider alternatives
    - **Close Scores:** Multiple good options available
    
    ---
    
    ## 📚 **Available Pipelines**
    
    ### 🥇 **DESeq2** (~60% of publications)
    - **Best for:** General use, batch effects, small samples
    - **Model:** Negative binomial + shrinkage
    - **Strengths:** Most conservative, well-documented
    
    ### 🥈 **edgeR** (~25% of publications)
    - **Best for:** Low counts, overdispersion, small samples
    - **Model:** Negative binomial + empirical Bayes
    - **Strengths:** Excellent for low counts, fast
    
    ### 🥉 **limma-voom** (~10% of publications)
    - **Best for:** Large samples, complex designs, speed
    - **Model:** Linear model + precision weights
    - **Strengths:** Very fast, efficient for large datasets
    
    ### 🔹 **Wilcoxon** (~5% of publications)
    - **Best for:** Non-parametric, no distributional assumptions
    - **Model:** Rank-based test
    - **Strengths:** Robust to outliers, simple
    
    ---
    
    ## **Recommendation Factors**
    
    **Sample Size** - Most Critical
    - Small (n<3): Conservative methods
    - Medium (n=3-7): Standard methods
    - Large (n≥8): All methods work
    - Very large (n>20): Speed matters (limma-voom)
    
    **BCV (Biological Coefficient of Variation)**
    - Low (<0.2): All methods work
    - Moderate (0.2-0.4): DESeq2 or edgeR
    - High (>0.4): edgeR handles better
    
    **Outliers**
    - Present + small samples: edgeR robust
    - Present + large samples: limma-voom robust
    - Absent: Standard methods
    
    **Low Count Genes**
    - Many (>30%): edgeR preferred
    - Few: All appropriate
    
    ---
    
    ## 🚀 **Workflow**
    ```
    Module 3 (Profile) → Extract 32 features
          ↓
    Module 4 (Recommend) → Rule-Based + ML
          ↓
    Compare & Choose → Best Pipeline
          ↓
    Module 7 (DE Import) → Run Analysis
    ```
    """)

# Check module availability
if not RULE_BASED_AVAILABLE:
    st.error(f"""
    **Rule-Based Recommender not available**
    
    Error: {rule_import_error}
    
    Please ensure RAPTOR is properly installed:
    ```bash
    cd RAPTOR/
    pip install -e .
    ```
    """)
    st.stop()

st.markdown("---")

# Section 1: Check for Data Profile
st.markdown("## 1. Data Profile")

profile = None
profile_source = None

# Check for profile from Module 3
if st.session_state.get('m3_complete', False):
    st.success("Using data profile from Module 3 (Profiler)")
    profile = st.session_state.get('m3_profile')
    profile_source = "module3"
    
    # Show key metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Samples", profile.n_samples)
        st.metric("Groups", profile.n_groups)
    
    with col2:
        st.metric("BCV", f"{profile.bcv:.3f}")
        st.caption(f"({profile.bcv_category})")
    
    with col3:
        st.metric("Min Group Size", profile.min_group_size)
        st.metric("Balance", f"{profile.sample_balance:.2f}")
    
    with col4:
        st.metric("Sparsity", f"{profile.sparsity*100:.1f}%")
        st.metric("Overdispersion", profile.overdispersion_category.upper())

else:
    st.info("No profile from Module 3. Please run Data Profiler first, or upload profile below.")
    
    # Option to upload profile JSON
    with st.expander("Upload profile JSON (from Module 3 export)"):
        uploaded_profile = st.file_uploader(
            "Upload data profile (JSON)",
            type=['json'],
            help="Use the JSON file exported from Module 3"
        )
        
        if uploaded_profile:
            try:
                import json
                profile_dict = json.load(uploaded_profile)
                
                # Convert to DataProfile if possible
                if PROFILER_AVAILABLE:
                    from dataclasses import fields
                    field_names = {f.name for f in fields(DataProfile)}
                    profile_data = {k: v for k, v in profile_dict.items() if k in field_names}
                    profile = DataProfile(**profile_data)
                else:
                    profile = profile_dict
                
                profile_source = "upload"
                st.success("Profile loaded successfully")
                
            except Exception as e:
                st.error(f"Error loading profile: {str(e)}")
                profile = None

st.markdown("---")

# Section 2: Check ML Model Availability
st.markdown("## 2. Recommender Setup")

ml_model_available = False
ml_model_path = None

# Check for ML model in common locations
possible_paths = [
    Path("models/ml_recommender.pkl"),
    Path("models/ml_recommender_random_forest.pkl"),
    Path("raptor/models/ml_recommender.pkl"),
    Path.cwd() / "models" / "ml_recommender.pkl"
]

for path in possible_paths:
    if path.exists():
        ml_model_available = True
        ml_model_path = str(path)
        break

# Display status
col1, col2 = st.columns(2)

with col1:
    st.markdown("### Rule-Based Recommender")
    if RULE_BASED_AVAILABLE:
        st.success("**Available** - Always ready")
        st.caption("Uses literature-based heuristics")
    else:
        st.error("Not available")

with col2:
    st.markdown("### ML-Based Recommender")
    if ML_RECOMMENDER_AVAILABLE and ml_model_available:
        st.success(f"**Available** - Model found")
        st.caption(f"Model: {ml_model_path}")
    elif ML_RECOMMENDER_AVAILABLE and not ml_model_available:
        st.warning("**Module Available** - No trained model")
        st.caption("Train a model with CLI: `raptor train-ml-recommender`")
    else:
        st.info("Not available")
        st.caption("Install: pip install scikit-learn joblib")

# Choose mode
if profile is not None:
    st.markdown("### Recommendation Mode")
    
    # Determine available modes
    modes = ["🔹 Rule-Based Only"]
    if ML_RECOMMENDER_AVAILABLE and ml_model_available:
        modes.extend(["ML-Based Only", "Compare Both (Recommended)"])
    
    mode = st.radio(
        "Select recommendation mode:",
        modes,
        index=len(modes) - 1 if len(modes) > 1 else 0,  # Default to compare if available
        help="Compare mode shows both recommendations side-by-side"
    )
    
    # Run button
    if st.button("Get Recommendation", type="primary", use_container_width=True):
        
        with st.spinner("Generating recommendations..."):
            
            try:
                # Rule-based recommendation
                if "Rule-Based" in mode or "Compare" in mode:
                    
                    status_text = st.empty()
                    status_text.text("Running rule-based analysis...")
                    
                    # Initialize rule-based recommender
                    rule_recommender = PipelineRecommender(profile)
                    rule_recommendation = rule_recommender.get_recommendation()
                    
                    # Store in session state
                    st.session_state['m4_rule_recommendation'] = rule_recommendation
                    st.session_state['m4_rule_complete'] = True
                    
                    status_text.text("Rule-based recommendation complete")
                
                # ML-based recommendation
                if ("ML-Based" in mode or "Compare" in mode) and ML_RECOMMENDER_AVAILABLE and ml_model_available:
                    
                    status_text = st.empty()
                    status_text.text("Running ML analysis...")
                    
                    # Initialize ML recommender
                    ml_recommender = MLPipelineRecommender(model_path=ml_model_path)
                    ml_recommendation = ml_recommender.recommend(profile)
                    
                    # Store in session state
                    st.session_state['m4_ml_recommendation'] = ml_recommendation
                    st.session_state['m4_ml_complete'] = True
                    
                    status_text.text("ML recommendation complete")
                
                # Set completion flags
                st.session_state['m4_complete'] = True
                st.session_state['m4_mode'] = mode
                st.session_state['m4_error'] = False
                
                st.success("""
                **Recommendation complete!**
                
                Scroll down to see results and analysis.
                """)
                
            except Exception as e:
                st.error(f"""
                **Error during recommendation**
                
                {str(e)}
                
                **Troubleshooting:**
                - Ensure profile has all required fields
                - For ML: Check model file is valid
                - Verify RAPTOR installation
                """)
                st.session_state['m4_error'] = True

# Section 3: Results
if st.session_state.get('m4_complete', False):
    st.markdown("---")
    st.markdown("## 3. Recommendations")
    
    mode = st.session_state.get('m4_mode', '')
    
    # Rule-Based Results
    if st.session_state.get('m4_rule_complete', False):
        rec = st.session_state['m4_rule_recommendation']
        
        st.markdown("### Rule-Based Recommendation")
        
        # Primary recommendation
        col1, col2, col3 = st.columns([2, 1, 1])
        
        with col1:
            st.markdown(f"#### 🥇 **{rec.primary_pipeline}**")
            st.markdown(f"**Reason:** {rec.primary_reason}")
        
        with col2:
            st.metric("Confidence", f"{rec.primary_score:.0f}%")
        
        with col3:
            # Get pipeline info
            info = get_pipeline_info(rec.primary_pipeline)
            if info:
                st.caption(f"{info.description}")
        
        # Alternative recommendation
        with st.expander("🥈 Alternative Recommendation"):
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown(f"**Pipeline:** {rec.alternative_pipeline}")
                st.markdown(f"**Reason:** {rec.alternative_reason}")
            
            with col2:
                st.metric("Confidence", f"{rec.alternative_score:.0f}%")
        
        # Decision factors
        with st.expander("Decision Factors"):
            factors_df = pd.DataFrame([rec.decision_factors]).T
            factors_df.columns = ['Value']
            factors_df.index.name = 'Factor'
            st.dataframe(factors_df, use_container_width=True)
        
        # Warnings
        if rec.warnings:
            with st.expander("Warnings & Recommendations"):
                for warning in rec.warnings:
                    st.warning(warning)
        
        # All scores
        with st.expander("All Pipeline Scores"):
            scores_df = pd.DataFrame([rec.all_scores]).T
            scores_df.columns = ['Score']
            scores_df.index.name = 'Pipeline'
            scores_df = scores_df.sort_values('Score', ascending=False)
            
            fig = go.Figure(go.Bar(
                x=scores_df.index,
                y=scores_df['Score'],
                marker_color=['gold' if i == 0 else 'silver' if i == 1 else 'lightblue' 
                             for i in range(len(scores_df))],
                text=scores_df['Score'].apply(lambda x: f"{x:.0f}%"),
                textposition='outside'
            ))
            
            fig.update_layout(
                title="Pipeline Score Comparison",
                xaxis_title="Pipeline",
                yaxis_title="Score (%)",
                height=400,
                showlegend=False
            )
            
            st.plotly_chart(fig, use_container_width=True)
        
        # R code
        with st.expander("💻 R Code (Primary Pipeline)"):
            st.code(rec.r_code_primary, language='r')
            
            st.download_button(
                "Download R Script",
                rec.r_code_primary,
                f"{rec.primary_pipeline.lower()}_analysis.R",
                "text/plain"
            )
    
    # ML Results
    if st.session_state.get('m4_ml_complete', False):
        ml_rec = st.session_state['m4_ml_recommendation']
        
        st.markdown("### ML-Based Recommendation")
        
        # Primary recommendation
        col1, col2, col3 = st.columns([2, 1, 1])
        
        with col1:
            st.markdown(f"#### 🥇 **{ml_rec.pipeline}**")
            st.caption(f"Model: {ml_rec.model_version}")
        
        with col2:
            st.metric("Confidence", f"{ml_rec.confidence*100:.1f}%")
        
        with col3:
            # Confidence category
            if ml_rec.confidence >= 0.8:
                st.success("High Confidence")
            elif ml_rec.confidence >= 0.6:
                st.info("Medium Confidence")
            else:
                st.warning("Low Confidence")
        
        # Alternatives
        with st.expander("🥈 Alternative Predictions"):
            alt_data = []
            for alt_pipeline, alt_prob in ml_rec.alternatives[:5]:
                alt_data.append({
                    'Pipeline': alt_pipeline,
                    'Probability': f"{alt_prob*100:.1f}%",
                    'Confidence': alt_prob
                })
            
            alt_df = pd.DataFrame(alt_data)
            st.dataframe(alt_df[['Pipeline', 'Probability']], use_container_width=True)
            
            # Chart
            fig = go.Figure(go.Bar(
                x=alt_df['Pipeline'],
                y=alt_df['Confidence'],
                marker_color=['gold' if i == 0 else 'lightblue' for i in range(len(alt_df))],
                text=alt_df['Probability'],
                textposition='outside'
            ))
            
            fig.update_layout(
                title="ML Prediction Probabilities",
                xaxis_title="Pipeline",
                yaxis_title="Probability",
                height=400,
                showlegend=False
            )
            
            st.plotly_chart(fig, use_container_width=True)
        
        # Feature importances
        with st.expander("Key Contributing Features"):
            st.markdown("**Features that most influenced the ML prediction:**")
            
            feat_data = []
            for feat, importance in list(ml_rec.feature_importances.items())[:10]:
                feat_data.append({
                    'Feature': feat,
                    'Importance': importance
                })
            
            feat_df = pd.DataFrame(feat_data)
            
            fig = go.Figure(go.Bar(
                x=feat_df['Importance'],
                y=feat_df['Feature'],
                orientation='h',
                marker_color='steelblue'
            ))
            
            fig.update_layout(
                title="Top 10 Feature Importances",
                xaxis_title="Importance",
                yaxis_title="Feature",
                height=400,
                showlegend=False
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    # Comparison Analysis (if both available)
    if st.session_state.get('m4_rule_complete', False) and st.session_state.get('m4_ml_complete', False):
        st.markdown("---")
        st.markdown("### Comparison Analysis")
        
        rec = st.session_state['m4_rule_recommendation']
        ml_rec = st.session_state['m4_ml_recommendation']
        
        # Check agreement
        agreement = rec.primary_pipeline == ml_rec.pipeline
        
        if agreement:
            st.success(f"""
            **BOTH RECOMMENDERS AGREE!**
            
            **Recommended Pipeline:** {rec.primary_pipeline}
            
            **Confidence:**
            - Rule-Based: {rec.primary_score:.0f}%
            - ML-Based: {ml_rec.confidence*100:.1f}%
            
            **Conclusion:** High confidence recommendation. Proceed with {rec.primary_pipeline}!
            """)
        else:
            st.warning(f"""
            **RECOMMENDERS DISAGREE**
            
            **Rule-Based recommends:** {rec.primary_pipeline} ({rec.primary_score:.0f}%)
            **ML-Based recommends:** {ml_rec.pipeline} ({ml_rec.confidence*100:.1f}%)
            
            **What this means:**
            - Multiple pipelines may perform similarly on your data
            - Consider both options and their reasoning
            - Check if recommended pipelines are close in scores
            
            **Recommendation:** Review both analyses below and choose based on your priorities.
            """)
        
        # Side-by-side comparison
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### Rule-Based")
            st.markdown(f"**Primary:** {rec.primary_pipeline}")
            st.markdown(f"**Score:** {rec.primary_score:.0f}%")
            st.markdown(f"**Reason:** {rec.primary_reason}")
            
            if rec.all_scores:
                st.markdown("**All Scores:**")
                for pipeline, score in sorted(rec.all_scores.items(), key=lambda x: x[1], reverse=True)[:3]:
                    st.caption(f"• {pipeline}: {score:.0f}%")
        
        with col2:
            st.markdown("#### ML-Based")
            st.markdown(f"**Primary:** {ml_rec.pipeline}")
            st.markdown(f"**Confidence:** {ml_rec.confidence*100:.1f}%")
            st.markdown(f"**Model:** {ml_rec.model_version}")
            
            if ml_rec.alternatives:
                st.markdown("**All Predictions:**")
                for pipeline, prob in ml_rec.alternatives[:3]:
                    st.caption(f"• {pipeline}: {prob*100:.1f}%")
        
        # Agreement matrix
        with st.expander("Detailed Score Comparison"):
            # Combine scores
            pipelines = list(rec.all_scores.keys())
            
            comparison_data = []
            for pipeline in pipelines:
                rule_score = rec.all_scores.get(pipeline, 0)
                ml_prob = 0
                if pipeline == ml_rec.pipeline:
                    ml_prob = ml_rec.confidence
                else:
                    for alt_pipe, alt_prob in ml_rec.alternatives:
                        if alt_pipe == pipeline:
                            ml_prob = alt_prob
                            break
                
                comparison_data.append({
                    'Pipeline': pipeline,
                    'Rule-Based Score': f"{rule_score:.0f}%",
                    'ML Probability': f"{ml_prob*100:.1f}%",
                    'Rule_Score_Val': rule_score,
                    'ML_Prob_Val': ml_prob * 100
                })
            
            comp_df = pd.DataFrame(comparison_data)
            comp_df = comp_df.sort_values('Rule_Score_Val', ascending=False)
            
            st.dataframe(comp_df[['Pipeline', 'Rule-Based Score', 'ML Probability']], use_container_width=True)
            
            # Comparison chart
            fig = go.Figure()
            
            fig.add_trace(go.Bar(
                name='Rule-Based',
                x=comp_df['Pipeline'],
                y=comp_df['Rule_Score_Val'],
                marker_color='steelblue'
            ))
            
            fig.add_trace(go.Bar(
                name='ML-Based',
                x=comp_df['Pipeline'],
                y=comp_df['ML_Prob_Val'],
                marker_color='orange'
            ))
            
            fig.update_layout(
                title="Rule-Based vs ML-Based Scores",
                xaxis_title="Pipeline",
                yaxis_title="Score (%)",
                height=500,
                barmode='group'
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    # Export recommendations
    st.markdown("---")
    st.markdown("### Export Recommendations")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.session_state.get('m4_rule_complete', False):
            rec = st.session_state['m4_rule_recommendation']
            summary_text = rec.summary()
            
            st.download_button(
                "Rule-Based Summary",
                summary_text,
                "rule_based_recommendation.txt",
                "text/plain"
            )
    
    with col2:
        if st.session_state.get('m4_ml_complete', False):
            ml_rec = st.session_state['m4_ml_recommendation']
            ml_summary = ml_rec.summary()
            
            st.download_button(
                "ML Summary",
                ml_summary,
                "ml_recommendation.txt",
                "text/plain"
            )
    
    with col3:
        if st.session_state.get('m4_rule_complete', False) and st.session_state.get('m4_ml_complete', False):
            # Create comparison report
            rec = st.session_state['m4_rule_recommendation']
            ml_rec = st.session_state['m4_ml_recommendation']
            
            comparison_report = f"""
RAPTOR Hybrid Recommendation Report
====================================

RULE-BASED RECOMMENDATION:
{rec.summary()}

ML-BASED RECOMMENDATION:
{ml_rec.summary()}

AGREEMENT ANALYSIS:
-------------------
Rule-Based Primary: {rec.primary_pipeline} ({rec.primary_score:.0f}%)
ML-Based Primary: {ml_rec.pipeline} ({ml_rec.confidence*100:.1f}%)

Agreement: {'YES ✓' if rec.primary_pipeline == ml_rec.pipeline else 'NO ✗'}

CONCLUSION:
-----------
{'Both recommenders agree on ' + rec.primary_pipeline + '. High confidence recommendation.' 
 if rec.primary_pipeline == ml_rec.pipeline 
 else 'Recommenders disagree. Consider both ' + rec.primary_pipeline + ' and ' + ml_rec.pipeline + '.'}
            """
            
            st.download_button(
                "Comparison Report",
                comparison_report,
                "hybrid_recommendation_report.txt",
                "text/plain"
            )
    
    # Next steps
    st.markdown("---")
    st.info(f"""
    **Recommendation complete!**
    
    **Recommended Pipeline:** {st.session_state['m4_rule_recommendation'].primary_pipeline if st.session_state.get('m4_rule_complete', False) else 'See above'}
    
    **Next steps:**
    1. Review the recommendation and reasoning
    2. Check agreement analysis (if both methods used)
    3. Download R code for your chosen pipeline
    4. Proceed to **Module 7: DE Import** to run the analysis
    5. Or use the R code directly in your own workflow
    
    **Note:** Module 7 will use your chosen pipeline to perform differential expression analysis.
    """)

else:
    if profile is None:
        st.info("Run Module 3 (Data Profiler) first to generate a data profile, then return here for recommendations.")
    else:
        st.info("Click 'Get Recommendation' above to see pipeline recommendations.")
