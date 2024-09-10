#!/usr/bin/env python
# coding: utf-8

# In[3]:
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os

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

if df is not None:
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
        'upregulated': '#FF6347',    # Red for upregulated pathways
        'downregulated': '#1E90FF',  # Blue for downregulated pathways
        'non-significant': '#A9A9A9' # Grey for non-significant pathways
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
            size=8,
            color=palette['non-significant'],
            opacity=0.8,
            line=dict(
                width=0.5,
                color='black'
            )
        ),
        text=[f'<span style="color:{palette["non-significant"]};">{name}</span>' for name in non_significant_df.index],
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
            size=8,
            color=palette['upregulated'],
            opacity=0.8,
            line=dict(
                width=0.5,
                color='black'
            )
        ),
        text=[f'<span style="color:{palette["upregulated"]};">{name}</span>' for name in upregulated_df.index],
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
            size=8,
            color=palette['downregulated'],
            opacity=0.8,
            line=dict(
                width=0.5,
                color='black'
            )
        ),
        text=[f'<span style="color:{palette["downregulated"]};">{name}</span>' for name in downregulated_df.index],
        hoverinfo='text',
        name='Downregulated'
    ))

    # Plot keyword matching pathways (interactive or not)
    keyword_df = df[df['category'] == 'keyword_match']
    
    if interactive:
        # If interactive, show pathway names on hover
        fig.add_trace(go.Scatter(
            x=keyword_df['GSVA_score'], 
            y=keyword_df['-log10(adj.P.Val)'], 
            mode='markers',
            marker=dict(
                size=15,
                color=palette['keyword_match'],
                opacity=0.8,
                line=dict(
                    width=0.5,
                    color='black'
                )
            ),
            text=[f'<span style="color:{palette["keyword_match"]};">{name}</span>' for name in keyword_df.index],
            hoverinfo='text',
            name=', '.join(keywords)  # Use keywords as the legend entry
        ))
    else:
        # If not interactive, display pathway names directly on the plot
        fig.add_trace(go.Scatter(
            x=keyword_df['GSVA_score'], 
            y=keyword_df['-log10(adj.P.Val)'], 
            mode='markers+text',
            marker=dict(
                size=15,
                color=palette['keyword_match'],
                opacity=0.8,
                line=dict(
                    width=0.5,
                    color='black'
                )
            ),
            text=[str(i+1) for i in range(len(keyword_df))],  # Number the pathways
            textposition='top center',
            textfont=dict(color='black', size=12),  # Bold black text for numbers
            hoverinfo='text',
            name=', '.join(keywords)  # Use keywords as the legend entry
        ))
    
    # Set layout with transparent background and white plot background
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',  # Transparent background
        plot_bgcolor='rgba(255,255,255,1)',  # White plot background
        title='Interactive Volcano Plot with Keyword Search',
        xaxis_title='GSVA Score',
        yaxis_title='-log10(adj.P.Val)',
        title_font_size=18,
        width=width,
        height=height,
        legend_title_text='Pathway Categories'
    )

    return fig
