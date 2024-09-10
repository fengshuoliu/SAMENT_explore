#!/usr/bin/env python
# coding: utf-8

# In[3]:
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
from plotly.io import to_image

# Check if kaleido is installed for image export
try:
    import kaleido
    kaleido_available = True
except ImportError:
    kaleido_available = False

# Load the data (Cached for performance)
@st.cache_data
def load_data():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, 'macrophages_biotin_positive-vs-negative_GSVA.csv')
    
    if not os.path.exists(file_path):
        st.error(f"File {file_path} not found. Please check the file path.")
        return None
    
    df = pd.read_csv(file_path)
    df = df.set_index(df.columns[0])
    df.index = df.index.str.strip()
    df['-log10(adj.P.Val)'] = -np.log10(df['P.Value'])
    return df

df = load_data()

# Function to categorize pathways
def get_category(row, keywords=[], logic='AND'):
    pathway_name = row.name.replace('_', ' ').upper()
    keywords = [kw.upper().strip() for kw in keywords if kw.strip() != '']
    
    if logic == 'AND':
        if all(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'
    elif logic == 'OR':
        if any(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'
    
    if row['GSVA_score'] > 0.5 and row['-log10(adj.P.Val)'] > -np.log10(0.05):
        return 'upregulated'
    elif row['GSVA_score'] < -0.5 and row['-log10(adj.P.Val)'] > -np.log10(0.05):
        return 'downregulated'
    else:
        return 'non-significant'

# Function to update the plot
def update_plot(keywords=[], logic='AND', width='100%', height=800, interactive=True):
    df['category'] = df.apply(get_category, axis=1, keywords=keywords, logic=logic)
    palette = {'keyword_match': '#32CD32', 'upregulated': '#FF6347', 'downregulated': '#1E90FF', 'non-significant': '#A9A9A9'}
    fig = go.Figure()
    
    # Plot non-significant pathways
    non_significant_df = df[df['category'] == 'non-significant']
    fig.add_trace(go.Scatter(x=non_significant_df['GSVA_score'], y=non_significant_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['non-significant'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[name for name in non_significant_df.index], hoverinfo='text', name='Non-Significant'))
    
    # Plot upregulated pathways
    upregulated_df = df[df['category'] == 'upregulated']
    fig.add_trace(go.Scatter(x=upregulated_df['GSVA_score'], y=upregulated_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['upregulated'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[f'<span style="color:{palette["upregulated"]};">{name}</span>' for name in upregulated_df.index], hoverinfo='text', name='Upregulated'))

    # Plot downregulated pathways
    downregulated_df = df[df['category'] == 'downregulated']
    fig.add_trace(go.Scatter(x=downregulated_df['GSVA_score'], y=downregulated_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['downregulated'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[f'<span style="color:{palette["downregulated"]};">{name}</span>' for name in downregulated_df.index], hoverinfo='text', name='Downregulated'))
    
    # Sort keyword-matched pathways by P.Value
    keyword_df = df[df['category'] == 'keyword_match'].sort_values('P.Value')
    
    if interactive:
        # Plot keyword-matched pathways interactively
        fig.add_trace(go.Scatter(x=keyword_df['GSVA_score'], y=keyword_df['-log10(adj.P.Val)'], mode='markers',
                                 marker=dict(size=15, color=palette['keyword_match'], opacity=0.8, line=dict(width=0.5, color='black')),
                                 text=[f'<span style="color:{palette["keyword_match"]};">{name}</span>' for name in keyword_df.index], hoverinfo='text', name=', '.join(keywords)))
    else:
        # Plot numbered keyword-matched pathways
        for i, (index, row) in enumerate(keyword_df.iterrows()):
            showlegend = i == 0  # Only show legend for the first trace
            fig.add_trace(go.Scatter(x=[row['GSVA_score']], y=[row['-log10(adj.P.Val)']], mode='text+markers',
                                     marker=dict(size=15, color=palette['keyword_match'], opacity=0.8, line=dict(width=0.5, color='black')),
                                     text=f"<b style='color:black;'>{i+1}</b>",  # Bold and black color for numbers
                                     hoverinfo='text', name=f"{', '.join(keywords)}" if showlegend else None, showlegend=showlegend))

    # Set layout
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(255,255,255,1)', title='Interactive Volcano Plot',
                      xaxis_title='GSVA Score', yaxis_title='-log10(adj.P.Val)', title_font_size=18, width=width, height=height,
                      legend_title_text='Pathway Categories', autosize=True)
    
    return fig, keyword_df

if df is not None:
    # Sidebar input for keywords and logic
    st.sidebar.header('Search Parameters')
    num_keywords = st.sidebar.number_input('Number of Keywords', min_value=1, max_value=10, value=2)
    keywords = [st.sidebar.text_input(f'Keyword {i+1}') for i in range(num_keywords)]
    logic = st.sidebar.selectbox('Logic', ['AND', 'OR'])
    keywords = [kw for kw in keywords if kw.strip() != '']
    fig_width = st.sidebar.slider('Figure Width', min_value=400, max_value=1200, value=1000, step=50)
    fig_height = st.sidebar.slider('Figure Height', min_value=400, max_value=1000, value=800, step=50)
    interactive_keywords = st.sidebar.radio('Keyword-Matched Pathways Interactive?', ('Yes', 'No'))
    
    # Update plot based on user input
    fig, keyword_df = update_plot(keywords, logic, width=fig_width, height=fig_height, interactive=(interactive_keywords == 'Yes'))
    
    # Display plot
    st.plotly_chart(fig, use_container_width=True)

    # Display the table for keyword-matched pathways if interactive is set to No
    if interactive_keywords == 'No' and not keyword_df.empty:
        st.write("### Keyword-Matched Pathways")
        st.dataframe(keyword_df[['P.Value']].reset_index().rename(columns={'index': 'Pathway'}))

    st.write(f"Keywords used: {keywords}")
    st.write(f"Logic used: {logic}")
    st.write(f"Figure size: {fig_width} x {fig_height}")
    st.write(f"Keyword-Matched Pathways Interactive: {interactive_keywords}")

    # Download plot as PNG or PDF
    st.sidebar.header('Download Plot')
    download_format = st.sidebar.radio('Download Format', ('PNG', 'PDF'))
    
    if st.sidebar.button('Download'):
        if kaleido_available:
            if download_format == 'PNG':
                file_bytes = to_image(fig, format='png', engine="kaleido", scale=3)  # 300 DPI
                st.sidebar.download_button(label='Download as PNG', data=file_bytes, file_name='plot.png', mime='image/png')
            elif download_format == 'PDF':
                file_bytes = to_image(fig, format='pdf', engine="kaleido", scale=3)  # 300 DPI
                st.sidebar.download_button(label='Download as PDF', data=file_bytes, file_name='plot.pdf', mime='application/pdf')
        else:
            st.sidebar.error("Image export requires the 'kaleido' package. Please install it by adding 'kaleido' to your requirements.txt file.")
