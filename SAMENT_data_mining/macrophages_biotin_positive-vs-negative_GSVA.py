#!/usr/bin/env python
# coding: utf-8

# In[3]:
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
from plotly.io import to_image

# Load the data (Cached for performance)
@st.cache_data
def load_data():
    # Get the absolute path of the current file (macrophages_biotin_positive-vs-negative_GSVA.py)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    # Construct the absolute path to the CSV file
    file_path = os.path.join(dir_path, 'macrophages_biotin_positive-vs-negative_GSVA.csv')
    
    if not os.path.exists(file_path):
        st.error(f"File {file_path} not found. Please check the file path.")
        return None
    
    df = pd.read_csv(file_path)
    df = df.set_index(df.columns[0])
    df.index = df.index.str.strip()  # Remove leading/trailing spaces from pathway names
    df['-log10(adj.P.Val)'] = -np.log10(df['P.Value'])
    return df

df = load_data()

# Define a function to categorize pathways based on keywords and logic (AND/OR)
def get_category(row, keywords=[], logic='AND'):
    pathway_name = row.name.replace('_', ' ').upper()
    
    # Ensure keywords are processed correctly (uppercased and stripped)
    keywords = [kw.upper().strip() for kw in keywords if kw.strip() != '']
    
    if logic == 'AND':
        # All keywords must be present in the pathway name
        if all(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'
    elif logic == 'OR':
        # Any keyword must be present in the pathway name
        if any(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'

    # Default logic for significance and GSVA score
    if row['GSVA_score'] > 0.5 and row['-log10(adj.P.Val)'] > -np.log10(0.05):
        return 'upregulated'
    elif row['GSVA_score'] < -0.5 and row['-log10(adj.P.Val)'] > -np.log10(0.05):
        return 'downregulated'
    else:
        return 'non-significant'

# Function to update the plot based on keyword input
def update_plot(keywords=[], logic='AND', width=800, height=600, interactive=True):
    # Apply category based on keyword search
    df['category'] = df.apply(get_category, axis=1, keywords=keywords, logic=logic)
    
    # Define color palette
    palette = {
        'keyword_match': '#32CD32',  # Green for pathways matching the search
        'upregulated': '#FF6347',  # Red for upregulated pathways
        'downregulated': '#1E90FF',  # Blue for downregulated pathways
        'non-significant': '#A9A9A9'  # Grey for non-significant pathways
    }

    # Create figure object
    fig = go.Figure()

    # Plot non-significant pathways
    non_significant_df = df[df['category'] == 'non-significant']
    fig.add_trace(go.Scatter(
        x=non_significant_df['GSVA_score'], 
        y=non_significant_df['-log10(adj.P.Val)'], 
        mode='markers',
        marker=dict(
            size=8,  # Dot size for other pathways
            color=palette['non-significant'],  # Grey color for non-significant pathways
            opacity=0.8,  # Set transparency
            line=dict(
                width=0.5,  # Set border thickness
                color='black'  # Set border color
            )
        ),
        text=[name for name in non_significant_df.index],  # Always interactive
        hoverinfo='text',
        name='Non-Significant'
    ))

    # Plot upregulated pathways
    upregulated_df = df[df['category'] == 'upregulated']
    fig.add_trace(go.Scatter(
        x=upregulated_df['GSVA_score'], 
        y=upregulated_df['-log10(adj.P.Val)'], 
        mode='markers',
        marker=dict(
            size=8,  # Dot size for upregulated pathways
            color=palette['upregulated'],  # Red color for upregulated pathways
            opacity=0.8,  # Set transparency
            line=dict(
                width=0.5,  # Set border thickness
                color='black'  # Set border color
            )
        ),
        text=[name for name in upregulated_df.index],  # Always interactive
        hoverinfo='text',
        name='Upregulated'
    ))

    # Plot downregulated pathways
    downregulated_df = df[df['category'] == 'downregulated']
    fig.add_trace(go.Scatter(
        x=downregulated_df['GSVA_score'], 
        y=downregulated_df['-log10(adj.P.Val)'], 
        mode='markers',
        marker=dict(
            size=8,  # Dot size for downregulated pathways
            color=palette['downregulated'],  # Blue color for downregulated pathways
            opacity=0.8,  # Set transparency
            line=dict(
                width=0.5,  # Set border thickness
                color='black'  # Set border color
            )
        ),
        text=[name for name in downregulated_df.index],  # Always interactive
        hoverinfo='text',
        name='Downregulated'
    ))

    # Plot keyword matching pathways
    keyword_df = df[df['category'] == 'keyword_match'].sort_values('P.Value')  # Sort by P.Value
    
    if interactive:
        # Interactive: Show name on hover
        fig.add_trace(go.Scatter(
            x=keyword_df['GSVA_score'], 
            y=keyword_df['-log10(adj.P.Val)'], 
            mode='markers',
            marker=dict(
                size=15,  # Dot size for keyword-matching pathways
                color=palette['keyword_match'],  # Green color for keyword-matching pathways
                opacity=0.8,  # Set transparency
                line=dict(
                    width=0.5,  # Set border thickness
                    color='black'  # Set border color
                )
            ),
            text=[name for name in keyword_df.index],  # Interactive: show name on hover
            hoverinfo='text',
            name='Keyword Matched Pathways'
        ))
    else:
        # Non-interactive: Always show name
        for i, (index, row) in enumerate(keyword_df.iterrows()):
            fig.add_trace(go.Scatter(
                x=[row['GSVA_score']], 
                y=[row['-log10(adj.P.Val)']], 
                mode='text+markers',
                marker=dict(
                    size=15,  # Dot size for keyword-matching pathways
                    color=palette['keyword_match'],  # Green color for keyword-matching pathways
                    opacity=0.8,  # Set transparency
                    line=dict(
                        width=0.5,  # Set border thickness
                        color='black'  # Set border color
                    )
                ),
                text=f"{i+1}",  # Number the pathway
                hoverinfo='text',
                name=f"Keyword Matched Pathways {i+1}"
            ))
    
    # Set layout with transparent background and white plot background
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',  # Transparent background
        plot_bgcolor='rgba(255,255,255,1)',  # White plot background
        title='Interactive Volcano Plot with Keyword Search',
        xaxis_title='GSVA Score',
        yaxis_title='-log10(adj.P.Val)',
        title_font_size=18,
        width=width,  # Custom width for figure size
        height=height,  # Custom height for figure size
        legend_title_text='Pathway Categories'  # Custom legend title
    )

    return fig, keyword_df

if df is not None:
    # Sidebar input for keywords and logic
    st.sidebar.header('Search Parameters')
    num_keywords = st.sidebar.number_input('Number of Keywords', min_value=1, max_value=10, value=2)

    # Gather keywords dynamically in the sidebar
    keywords = [st.sidebar.text_input(f'Keyword {i+1}') for i in range(num_keywords)]

    # Allow user to select logic (AND or OR)
    logic = st.sidebar.selectbox('Logic', ['AND', 'OR'])

    # Filter out empty keywords
    keywords = [kw for kw in keywords if kw.strip() != '']

    # Allow user to adjust figure size
    fig_width = st.sidebar.slider('Figure Width', min_value=400, max_value=1200, value=800, step=50)
    fig_height = st.sidebar.slider('Figure Height', min_value=400, max_value=1000, value=600, step=50)

    # Allow user to toggle interactive mode for keyword-matched pathways
    interactive_keywords = st.sidebar.radio('Keyword-Matched Pathways Interactive?', ('Yes', 'No'))

    # Show plot in Streamlit app
    fig, keyword_df = update_plot(keywords, logic, width=fig_width, height=fig_height, interactive=(interactive_keywords == 'Yes'))

    st.plotly_chart(fig)

    if interactive_keywords == 'No' and not keyword_df.empty:
        # Show a table of the keyword-matched pathways sorted by P.Value
        st.write("### Keyword-Matched Pathways")
        st.dataframe(keyword_df[['P.Value']].reset_index().rename(columns={'index': 'Pathway'}))

    # Display search info
    st.write(f"Keywords used: {keywords}")
    st.write(f"Logic used: {logic}")
    st.write(f"Figure size: {fig_width} x {fig_height}")
    st.write(f"Keyword-Matched Pathways Interactive: {interactive_keywords}")

    # Allow user to download the plot as PNG or PDF
    st.sidebar.header('Download Plot')
    download_format = st.sidebar.radio('Download Format', ('PNG', 'PDF'))

    if st.sidebar.button('Download'):
        if download_format == 'PNG':
            file_bytes = to_image(fig, format='png', engine="kaleido", scale=3)  # 300 DPI
            st.sidebar.download_button(label='Download as PNG', data=file_bytes, file_name='plot.png', mime='image/png')
        elif download_format == 'PDF':
            file_bytes = to_image(fig, format='pdf', engine="kaleido", scale=3)  # 300 DPI
            st.sidebar.download_button(label='Download as PDF', data=file_bytes, file_name='plot.pdf', mime='application/pdf')
