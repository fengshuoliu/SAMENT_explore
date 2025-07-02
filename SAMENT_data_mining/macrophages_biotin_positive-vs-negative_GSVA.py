import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
from plotly.io import to_image
import textwrap
import streamlit.components.v1 as components

# Insert the Google Analytics tracking code
components.html("""
<!-- Google tag (gtag.js) -->
<script async src="https://www.googletagmanager.com/gtag/js?id=G-ST1X16N2KC"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'G-ST1X16N2KC');
</script>
""", height=0, width=0)

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

# MSigDB Prefixes (splitting prefixes properly)
prefixes = ['HALLMARK', 'CHR', 'KEGG', 'REACTOME', 'BIOCARTA', 'CGP', 'MIR', 'TFT', 
            'CGN', 'CM', 'GOBP', 'GOCC', 'GOMF', 'C6', 'ONCOGENIC', 'C7', 'C8']

# Function to categorize pathways dynamically based on user-set cutoffs
def get_category(row, neg_cutoff, pos_cutoff, significance_cutoff, keywords=[], exclude_keywords=[], logic='AND'):
    pathway_name = row.name.replace('_', ' ').upper()
    keywords = [kw.upper().strip() for kw in keywords if kw.strip() != '']
    exclude_keywords = [kw.upper().strip() for kw in exclude_keywords if kw.strip() != '']
    
    # Exclude pathways that match any excluded keywords
    if any(keyword in pathway_name for keyword in exclude_keywords):
        return 'non-significant'
    
    if logic == 'AND':
        if all(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'
    elif logic == 'OR':
        if any(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'
    
    if row['GSVA_score'] > pos_cutoff and row['-log10(adj.P.Val)'] > significance_cutoff:
        return 'upregulated'
    elif row['GSVA_score'] < neg_cutoff and row['-log10(adj.P.Val)'] > significance_cutoff:
        return 'downregulated'
    else:
        return 'non-significant'

# Function to wrap long text in legend
def wrap_text(text, width=30):
    return '<br>'.join(textwrap.wrap(text, width=width))

# Function to update the plot
def update_plot(keywords=[], exclude_keywords=[], logic='AND',
                neg_cutoff=-0.2, pos_cutoff=0.2, significance_cutoff=1.0,
                width='100%', height=800, interactive=True):
    df['category'] = df.apply(
    lambda row: get_category(
          row,
          neg_cutoff=neg_cutoff,
          pos_cutoff=pos_cutoff,
          significance_cutoff=significance_cutoff,
          keywords=keywords,
          exclude_keywords=exclude_keywords,
          logic=logic
        ),
        axis=1
    )

    palette = {'keyword_match': '#32CD32', 'upregulated': '#FF6347', 'downregulated': '#1E90FF', 'non-significant': '#A9A9A9'}
    fig = go.Figure()
    
    # Plot non-significant pathways
    non_significant_df = df[df['category'] == 'non-significant']
    fig.add_trace(go.Scatter(x=non_significant_df['GSVA_score'], y=non_significant_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['non-significant'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[name for name in non_significant_df.index], hoverinfo='text', name='Non-Significant'))
    
    # Plot up-regulated pathways
    upregulated_df = df[df['category'] == 'upregulated']
    fig.add_trace(go.Scatter(x=upregulated_df['GSVA_score'], y=upregulated_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['upregulated'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[f'<span style="color:{palette["upregulated"]};">{name}</span>' for name in upregulated_df.index], hoverinfo='text', name='Up-regulated'))

    # Plot down-regulated pathways
    downregulated_df = df[df['category'] == 'downregulated']
    fig.add_trace(go.Scatter(x=downregulated_df['GSVA_score'], y=downregulated_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['downregulated'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[f'<span style="color:{palette["downregulated"]};">{name}</span>' for name in downregulated_df.index], hoverinfo='text', name='Down-regulated'))
    
    # Sort keyword-matched pathways by P.Value
    keyword_df = df[df['category'] == 'keyword_match'].sort_values('P.Value')
    
    if interactive:
        # Plot keyword-matched pathways interactively
        fig.add_trace(go.Scatter(x=keyword_df['GSVA_score'], y=keyword_df['-log10(adj.P.Val)'], mode='markers',
                                 marker=dict(size=15, color=palette['keyword_match'], opacity=0.8, line=dict(width=0.5, color='black')),
                                 text=[f'<span style="color:{palette["keyword_match"]};">{name}</span>' for name in keyword_df.index],
                                 hoverinfo='text', name=wrap_text(', '.join(keywords), width=20)))
    else:
        # Plot numbered keyword-matched pathways (Start numbering from 1)
        for i, (index, row) in enumerate(keyword_df.iterrows(), start=1):
            showlegend = i == 1  # Only show legend for the first trace
            fig.add_trace(go.Scatter(x=[row['GSVA_score']], y=[row['-log10(adj.P.Val)']], mode='text+markers',
                                     marker=dict(size=15, color=palette['keyword_match'], opacity=0.8, line=dict(width=0.5, color='black')),
                                     text=f"<b style='color:black;'>{i}</b>",  # Bold and black color for numbers, starting from 1
                                     hoverinfo='text', name=wrap_text(f"{', '.join(keywords)}") if showlegend else None, showlegend=showlegend))

    # Add vertical and horizontal dashed lines based on user inputs
    fig.update_layout(
        shapes=[
            dict(
                type="line",
                x0=neg_cutoff,
                y0=0,
                x1=neg_cutoff,
                y1=df['-log10(adj.P.Val)'].max(),
                line=dict(
                    color="grey",
                    width=2,
                    dash="dash",
                ),
                layer='above'
            ),
            dict(
                type="line",
                x0=pos_cutoff,
                y0=0,
                x1=pos_cutoff,
                y1=df['-log10(adj.P.Val)'].max(),
                line=dict(
                    color="grey",
                    width=2,
                    dash="dash",
                ),
                layer='above'
            ),
            dict(
                type="line",
                x0=-df['GSVA_score'].max(),
                y0=significance_cutoff,
                x1=df['GSVA_score'].max(),
                y1=significance_cutoff,
                line=dict(
                    color="grey",
                    width=2,
                    dash="dash",
                ),
                layer='above'
            )
        ]
    )

    # Set layout with wrapped legend text
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(255,255,255,1)',
        title='Macrophage: Biotin Positive vs Negative',
        xaxis_title='GSVA Score',
        yaxis_title='-log10(adj.P.Val)',
        title_font_size=18,
        width=width,
        height=height,
        legend_title_text='Pathway Categories',
        autosize=True,
        legend=dict(
            x=1, y=1, 
            bgcolor='rgba(255,255,255,0)',  # Transparent background for legend
            bordercolor='rgba(255,255,255,0)',  # Transparent border
        )
    )
    
    return fig, keyword_df

if df is not None:
    # Sidebar input for cutoffs and other parameters
    st.sidebar.header('Cutoffs')
    neg_cutoff = st.sidebar.number_input('Negative GSVA Score Cutoff', value=-0.2, step=0.1)
    pos_cutoff = st.sidebar.number_input('Positive GSVA Score Cutoff', value=0.2, step=0.1)
    significance_cutoff = st.sidebar.number_input('Significance Cutoff', value=1.0, step=0.1)
    
    # Sidebar input for keywords and logic
    st.sidebar.header('Search Parameters')
    num_keywords = st.sidebar.number_input('Number of Keywords', min_value=1, max_value=10, value=2)
    keywords = [st.sidebar.text_input(f'Keyword {i+1}') for i in range(num_keywords)]
    logic = st.sidebar.selectbox('Logic', ['AND', 'OR'])
    exclude_keywords = st.sidebar.text_area('Exclude Keywords (comma-separated)').split(',')
    keywords = [kw for kw in keywords if kw.strip() != '']
    
    fig_width = st.sidebar.slider('Figure Width', min_value=400, max_value=1600, value=1000, step=50)
    fig_height = st.sidebar.slider('Figure Height', min_value=400, max_value=1200, value=800, step=50)
    
    # Prefix selection for filtering pathways
    st.sidebar.header('Filter Pathways by Prefix')
    selected_prefixes = st.sidebar.multiselect('Select Prefixes', prefixes, default=prefixes)
    
    # Filter the data based on selected prefixes
    if selected_prefixes:
        pattern = '|'.join(selected_prefixes)  # Create a pattern to match selected prefixes
        df = df[df.index.str.contains(pattern)]

    interactive_keywords = st.sidebar.radio('Keyword-Matched Pathways Interactive?', ('Yes', 'No'), index=0)  # Default to 'No'
    
    # Update plot based on user input
    fig, keyword_df = update_plot(keywords, exclude_keywords, logic, neg_cutoff, pos_cutoff, significance_cutoff, width=fig_width, height=fig_height, interactive=(interactive_keywords == 'Yes'))
    
    # Display plot
    st.plotly_chart(fig, use_container_width=True)

    # Display the table for keyword-matched pathways if interactive is set to No
    if interactive_keywords == 'No' and not keyword_df.empty:
        st.write("### Keyword-Matched Pathways")
        keyword_df_display = keyword_df[['P.Value']].reset_index().rename(columns={'index': 'Pathway'})
        keyword_df_display.index += 1  # Ensure the table starts numbering from 1
        st.dataframe(keyword_df_display)
                          
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
