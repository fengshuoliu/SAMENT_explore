#!/usr/bin/env python
# coding: utf-8

# In[5]:


import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go

# Load the data (Cached for performance)
@st.cache
def load_data():
    df = pd.read_csv('macrophages_biotin_positive-vs-negative_GSVA.csv')  # Adjust the file path as needed
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
def update_plot(keywords=[], logic='AND', width=800, height=600):
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
        text=[f'<span style="color:{palette["non-significant"]};">{name}</span>' for name in non_significant_df.index],  # Hover text with pathway names
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
        text=[f'<span style="color:{palette["upregulated"]};">{name}</span>' for name in upregulated_df.index],  # Hover text with pathway names
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
        text=[f'<span style="color:{palette["downregulated"]};">{name}</span>' for name in downregulated_df.index],  # Hover text with pathway names
        hoverinfo='text',
        name='Downregulated'
    ))

    # Plot keyword matching pathways on top
    keyword_df = df[df['category'] == 'keyword_match']
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
        text=[f'<span style="color:{palette["keyword_match"]};">{name}</span>' for name in keyword_df.index],  # Hover text with pathway names
        hoverinfo='text',
        name='Keyword Matched Pathways'
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

    return fig

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

# Show plot in Streamlit app
st.plotly_chart(update_plot(keywords, logic, width=fig_width, height=fig_height))

# Display search info
st.write(f"Keywords used: {keywords}")
st.write(f"Logic used: {logic}")
st.write(f"Figure size: {fig_width} x {fig_height}")

