import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import re
from scipy.stats import zscore
import matplotlib.lines as mlines
from adjustText import adjust_text  # Optional for reducing label overlap
from matplotlib.lines import Line2D

# =============================================================================
# ANALYSIS OF THE DATA
# =============================================================================
pathogenic_df = pd.read_csv("final_pathogenic_df.csv")
non_pathogenic_df = pd.read_csv("non_pathogenic_df.csv")


ethnic_map = {
    'British': 'White',
    'Any other white background': 'White',
    'Irish': 'White',
    'White': 'White',
    'Indian': 'Asian',
    'Pakistani': 'Asian',
    'Bangladeshi': 'Asian',
    'Chinese': 'Asian',
    'Any other Asian background': 'Asian',
    'Asian or Asian British': 'Asian',
    'African': 'Black or African American',
    'Caribbean': 'Black or African American',
    'Any other Black background': 'Black or African American',
    'Black or Black British': 'Black or African American',
    'White and Asian': 'Other',
    'White and Black Caribbean': 'Other',
    'White and Black African': 'Other',
    'Any other mixed background': 'Other',
    'Mixed': 'Other',
    'Other ethnic group': 'Other',
    'Do not know': 'Other',
    'Prefer not to answer': 'Other'
}

# Apply the map to BOTH DataFrames **right after reading/creating them**
pathogenic_df['Ethnic Category'] = pathogenic_df['Ethnic background'].map(ethnic_map)
non_pathogenic_df['Ethnic Category'] = non_pathogenic_df['Ethnic background'].map(ethnic_map)

# Add a label to each group for clarity
pathogenic_df['Group'] = 'Pathogenic'
non_pathogenic_df['Group'] = 'Non-Pathogenic'

# Combine for overall stats
combined_df = pd.concat([pathogenic_df, non_pathogenic_df], ignore_index=True)

# Total participants
n_pathogenic = len(pathogenic_df)
n_non_pathogenic = len(non_pathogenic_df)
n_total = len(combined_df)

# Age stats
age_mean_path = pathogenic_df['Age at recruitment'].mean()
age_median_path = pathogenic_df['Age at recruitment'].median()

age_mean_non_path = non_pathogenic_df['Age at recruitment'].mean()
age_median_non_path = non_pathogenic_df['Age at recruitment'].median()

# Sex counts
sex_counts_path = pathogenic_df['Sex'].value_counts()
sex_counts_non_path = non_pathogenic_df['Sex'].value_counts()
sex_counts_total = combined_df['Sex'].value_counts()

# Use the **recoded column** for ethnic counts!
ethnic_counts_path = pathogenic_df['Ethnic Category'].value_counts()
ethnic_counts_non_path = non_pathogenic_df['Ethnic Category'].value_counts()
ethnic_counts_total = combined_df['Ethnic Category'].value_counts()

# Prepare the summary rows
summary = []

# Add total participants
summary.append(['Number of participants', f'{n_pathogenic}', f'{n_non_pathogenic}'])

# Add age mean and median
summary.append(['Mean age at recruitment', f'{age_mean_path:.1f}', f'{age_mean_non_path:.1f}'])
summary.append(['Median age at recruitment', f'{age_median_path:.1f}', f'{age_median_non_path:.1f}'])

# Add sex breakdown with local and global percentages
for sex in sex_counts_total.index:
    path_n = sex_counts_path.get(sex, 0)
    non_path_n = sex_counts_non_path.get(sex, 0)
    local_path_pct = path_n / n_pathogenic * 100 if n_pathogenic > 0 else 0
    local_non_path_pct = non_path_n / n_non_pathogenic * 100 if n_non_pathogenic > 0 else 0
    global_pct = sex_counts_total[sex] / n_total * 100 if n_total > 0 else 0
    
    path_str = f'{path_n} ({local_path_pct:.1f}%, {global_pct:.1f}%)'
    non_path_str = f'{non_path_n} ({local_non_path_pct:.1f}%, {global_pct:.1f}%)'
    
    summary.append([f'Sex: {sex}', path_str, non_path_str])

for eth in ethnic_counts_total.index:
    path_n = ethnic_counts_path.get(eth, 0)
    non_path_n = ethnic_counts_non_path.get(eth, 0)
    local_path_pct = path_n / n_pathogenic * 100 if n_pathogenic > 0 else 0
    local_non_path_pct = non_path_n / n_non_pathogenic * 100 if n_non_pathogenic > 0 else 0
    global_pct = ethnic_counts_total[eth] / n_total * 100 if n_total > 0 else 0
    
    path_str = f'{path_n} ({local_path_pct:.1f}%, {global_pct:.1f}%)'
    non_path_str = f'{non_path_n} ({local_non_path_pct:.1f}%, {global_pct:.1f}%)'
    
    summary.append([f'Ethnic background: {eth}', path_str, non_path_str])

# Create the final summary table
summary_df = pd.DataFrame(summary, columns=['Characteristic', 'Pathogenic', 'Non-Pathogenic'])


# Define the biomarker columns
biomarker_columns = [
    'C-reactive protein', 
    'Rheumatoid factor',
    'Neutrophill count', 
    'Monocyte count', 
    'White blood cell (leukocyte) count',
    'Lymphocyte count',
    'Red blood cell (erythrocyte) count'
]
# Create an empty dictionary to store p-values
comparison_results = {}
# Loop through each biomarker to perform the appropriate test
for biomarker in biomarker_columns:
    # Extract data for pathogenic and non-pathogenic
    pathogenic_data = pathogenic_df[biomarker].dropna()
    non_pathogenic_data = non_pathogenic_df[biomarker].dropna()

    # Check if there are enough data points in each group to perform the test
    if len(pathogenic_data) < 3 or len(non_pathogenic_data) < 3:
        print(f"Skipping {biomarker} due to insufficient data in one or both groups.")
        continue  # Skip biomarker if not enough data

    # Calculate the means for reference
    pathogenic_averages = pathogenic_data.mean()
    non_pathogenic_averages = non_pathogenic_data.mean()

    # Perform Mann-Whitney U test
    stat, p_value = stats.mannwhitneyu(pathogenic_data, non_pathogenic_data, alternative='two-sided')
    # Format the p-value to three significant figures
    formatted_p_value = f"{p_value:.3g}"
    # Store the results
    comparison_results[biomarker] = {'p-value': formatted_p_value, 
                                      'Pathogenic': pathogenic_averages, 
                                      'Non-Pathogenic': non_pathogenic_averages}

# Convert the results into a DataFrame for easy viewing
comparison_results_df = pd.DataFrame(comparison_results).T
print(comparison_results_df)



# -------------------------------------------------------------------------------
# Split data into gene counts and consequences as before
consequence_counts = pathogenic_df.groupby(['SYMBOL', 'Consequence']).size().unstack(fill_value=0)
gene_totals = consequence_counts.sum(axis=1)
sorted_gene_totals = gene_totals.sort_values(ascending=False)
sorted_consequence_counts = consequence_counts.loc[sorted_gene_totals.index]
total_variants = gene_totals.sum()
top_consequences = sorted_consequence_counts.columns[sorted_consequence_counts.sum(axis=0) > 0]

fig, ax = plt.subplots(figsize=(15, 8), dpi=300)

# Plot
sorted_consequence_counts[top_consequences].plot(
    kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.6  # narrower bars
)

# Set y-axis limit
ax.set_ylim(0, 60)
ax.set_xlim(-1, len(sorted_consequence_counts.index) - 0.5)
# Align x-ticks with bars and space out labels
ax.set_xticks(range(len(sorted_consequence_counts.index)))
ax.set_xticklabels(sorted_consequence_counts.index, rotation=45, ha='right', fontsize=10)

# Annotations
for i, gene in enumerate(sorted_consequence_counts.index):
    total = gene_totals[gene]
    if total <= 60:
        percentage = (total / total_variants) * 100
        ax.text(i, total + 1.5, f'{percentage:.1f}%\n({total})',
                ha='center', va='bottom', fontsize=8, color='black')

# Labels
ax.set_xlabel('Gene')
ax.set_ylabel('Count of Variants')

# Legend
handles, labels = ax.get_legend_handles_labels()
new_labels = [label.replace('&', '\n&') for label in labels]
ax.legend(handles, new_labels, title='Consequence',
          bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

plt.tight_layout()
plt.savefig('Consequences_genes.svg', bbox_inches='tight')
plt.show()



df = pathogenic_df.copy()

# Drop rows where VAF or gnomADe_AF is missing
df = df.dropna(subset=['VAF', 'gnomADe_AF'])

# Filter genes of interest
filtered_df = df[df['VAF'] <= 0.25]
unique_genes = filtered_df['SYMBOL'].unique()

# Marker and color maps
markers = ['o', 's', '^', 'v', 'D', 'P', 'X', '*', 'H', 'd', '<', '>', 'p']
marker_map = {gene: markers[i % len(markers)] for i, gene in enumerate(unique_genes)}
palette = sns.color_palette("husl", len(unique_genes))
color_map = dict(zip(unique_genes, palette))

# --- Plot styling ---
sns.set_theme(style="whitegrid")

fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

# Scatter plot
ax.scatter(
    df['gnomADe_AF'],
    df['VAF'],
    color='steelblue',
    edgecolors='black',
    s=80,
    linewidths=0.5,
    alpha=0.7
)

# Log scale x-axis
ax.set_xscale('log')
ax.set_xlabel('gnomAD Allele Frequency', fontsize=14, weight='bold')
ax.set_ylabel('Variant Allele Fraction (VAF)', fontsize=14, weight='bold')

# Y-axis limits and formatting
ax.set_ylim(0, 0.7)
ax.set_yticks(np.arange(0, 0.8, 0.1))
ax.tick_params(axis='both', which='major', labelsize=12)

# Get safe x-axis limits after setting log scale
xmin, xmax = ax.get_xlim()
threshold = 0.25

# Horizontal threshold line
ax.axhline(y=threshold, color='red', linestyle='--', linewidth=1.5)

# Shaded regions
ax.fill_between(
    x=[xmin, xmax],
    y1=0,
    y2=threshold,
    color='red',
    alpha=0.1,
    label='Likely Mosaic'
)

ax.fill_between(
    x=[xmin, xmax],
    y1=threshold,
    y2=ax.get_ylim()[1],
    color='blue',
    alpha=0.05,
    label='Likely Germline'
)

# Annotation text
ax.text(
    xmax * 0.95, threshold - 0.02,
    'Likely Mosaic',
    color='red',
    fontsize=12,
    ha='right',
    va='top',
    weight='bold'
)

ax.text(
    xmax * 0.95, threshold + 0.02,
    'Likely Germline',
    color='blue',
    fontsize=12,
    ha='right',
    va='bottom',
    weight='bold'
)

# Clean up spines
sns.despine()

# Optional legend for regions
ax.legend(frameon=False, fontsize=12, loc='upper right')

# Tight layout for clean saving
plt.tight_layout()
# Save as PDF with high DPI and tight layout
plt.savefig('vaf_scatterplot_genes.svg', bbox_inches='tight')

# Show or save
plt.show()


# -------------------------------------------------------------------------------


# Combine AF columns and drop rows without VAF or gnomAD AF
df = pathogenic_df.copy()
df_patho = df.dropna(subset=['VAF', 'SYMBOL', 'HGVSc']).copy()
df_nonpatho = non_pathogenic_df.dropna(subset=['VAF', 'SYMBOL', 'HGVSc']).copy()

# 2. Add source labels
df_patho['source'] = 'pathogenic'
df_nonpatho['source'] = 'non_pathogenic'

# 3. Add a unique identifier (e.g., SYMBOL + HGVSc)
df_patho['variant_id'] = df_patho['SYMBOL'] + '|' + df_patho['HGVSc']
df_nonpatho['variant_id'] = df_nonpatho['SYMBOL'] + '|' + df_nonpatho['HGVSc']

# 4. Combine for population z-scoring
combined_df = pd.concat([df_patho, df_nonpatho], ignore_index=True)

# 5. Z-score per gene
combined_df['VAF_zscore'] = combined_df.groupby('SYMBOL')['VAF'].transform(zscore)

# 6. Mark outliers
threshold = 2
combined_df['is_outlier'] = combined_df['VAF_zscore'].abs() > threshold

# 7. Extract only outliers from pathogenic
outliers_df = combined_df[(combined_df['is_outlier']) & (combined_df['source'] == 'pathogenic')].copy()

# 8. Extract normal (non-outlier) pathogenic variants by variant_id
outlier_ids = set(outliers_df['variant_id'])
normal_df = df_patho[~df_patho['variant_id'].isin(outlier_ids)].copy()

# 8. Sort SYMBOLs alphabetically and enforce categorical order
sorted_symbols = sorted(df_patho['SYMBOL'].unique())
normal_df['SYMBOL'] = pd.Categorical(normal_df['SYMBOL'], categories=sorted_symbols, ordered=True)
outliers_df['SYMBOL'] = pd.Categorical(outliers_df['SYMBOL'], categories=sorted_symbols, ordered=True)

# 9. Plotting
fig, ax = plt.subplots(figsize=(12, 14), dpi=300)

# Define palette
palette = dict(zip(sorted_symbols, sns.color_palette('tab20', len(sorted_symbols))))

# Plot normal points
sns.scatterplot(data=normal_df, x='VAF', y='SYMBOL', hue='SYMBOL',
                palette=palette, s=80, ax=ax, marker='o', legend=False)

# Plot outliers
sns.scatterplot(data=outliers_df, x='VAF', y='SYMBOL', hue='SYMBOL',
                palette=palette, s=120, ax=ax, marker='X', legend=False)

# 10. Custom legend for outliers
legend_data = outliers_df.groupby('SYMBOL')['HGVSc'].apply(
    lambda x: ',\n'.join(sorted(set([v.split(':')[-1] for v in x])))
)
legend_labels = [
    f"{symbol}:\n{variants}"
    for symbol, variants in legend_data.items()
    if variants and variants.strip()
]
handles = [mlines.Line2D([0], [0], color='w', label=label) for label in legend_labels]

ax.legend(handles=handles, title='Outlier Variants (|z| > 2)', loc='upper left',
          bbox_to_anchor=(1.05, 1), ncol=2, fontsize=10)

# 11. Axes and layout
ax.set_xlabel('VAF', fontsize=12)
ax.set_ylabel('SYMBOL', fontsize=12)
ax.set_yticks(sorted_symbols)  # ensures only valid symbols show

plt.subplots_adjust(left=0.1, right=0.75, top=0.95, bottom=0.1)
plt.show()


# -------------------------------------------------------------------------------

def plot_variant_outliers(pathogenic_df, non_pathogenic_df, variable_col, variable_label, threshold_value=None):
    df_patho = pathogenic_df.dropna(subset=['VAF', 'SYMBOL', variable_col, 'HGVSc']).copy()
    df_nonpatho = non_pathogenic_df.dropna(subset=['VAF', 'SYMBOL', variable_col, 'HGVSc']).copy()

    # Add source labels
    df_patho['source'] = 'pathogenic'
    df_nonpatho['source'] = 'non_pathogenic'

    # Unique variant ID
    df_patho['variant_id'] = df_patho['SYMBOL'] + '|' + df_patho['HGVSc']
    df_nonpatho['variant_id'] = df_nonpatho['SYMBOL'] + '|' + df_nonpatho['HGVSc']

    # Combine and clean
    combined_df = pd.concat([df_patho, df_nonpatho], ignore_index=True)
    combined_df[variable_col] = pd.to_numeric(combined_df[variable_col], errors='coerce')
    combined_df.dropna(subset=[variable_col, 'VAF'], inplace=True)

    # Z-score and outlier marking
    combined_df['VAF_zscore'] = zscore(combined_df[variable_col])
    combined_df['is_outlier'] = combined_df['VAF_zscore'].abs() > 3

    # Extract short HGVS annotations
    combined_df['HGVSc_short'] = combined_df['HGVSc'].apply(lambda x: x.split(':')[-1] if pd.notna(x) else '')
    combined_df['HGVSp_short'] = combined_df['HGVSp'].apply(lambda x: x.split(':')[-1] if pd.notna(x) else None)

    # Use HGVSp_short if available, else HGVSc_short
    combined_df['HGVS_label'] = combined_df.apply(
        lambda row: row['HGVSp_short'] if pd.notnull(row['HGVSp_short']) else row['HGVSc_short'],
        axis=1
    )

    # Extract pathogenic outliers
    outliers_df = combined_df[(combined_df['is_outlier']) & (combined_df['source'] == 'pathogenic')].copy()

    # VAF categorization
    def categorize_vaf(vaf):
        if vaf < 0.3:
            return 'Mosaic (VAF < 0.3)'
        elif vaf < 0.7:
            return 'Germline Heterozygous (0.3 ≤ VAF < 0.7)'
        else:
            return 'Germline Homozygous (VAF ≥ 0.7)'

    df_patho['VAF_category'] = df_patho['VAF'].apply(categorize_vaf)
    outliers_df['VAF_category'] = outliers_df['VAF'].apply(categorize_vaf)

    # Color palette
    vaf_palette = {
        'Mosaic (VAF < 0.3)': '#1f77b4',
        'Germline Heterozygous (0.3 ≤ VAF < 0.7)': '#2ca02c',
        'Germline Homozygous (VAF ≥ 0.7)': '#d62728'
    }

    # Categorical y-axis ordering
    sorted_symbols = sorted(df_patho['SYMBOL'].unique())
    df_patho['SYMBOL'] = pd.Categorical(df_patho['SYMBOL'], categories=sorted_symbols, ordered=True)
    outliers_df['SYMBOL'] = pd.Categorical(outliers_df['SYMBOL'], categories=sorted_symbols, ordered=True)

    # Plot setup
    fig, ax = plt.subplots(figsize=(14, 12), dpi=300)

    sns.scatterplot(
        data=df_patho,
        x=variable_col,
        y='SYMBOL',
        hue='VAF_category',
        palette=vaf_palette,
        edgecolor='black',
        linewidth=0.5,
        alpha=0.8,
        legend=False,
        ax=ax
    )

    # Plot outliers with annotation
    texts = []
    for _, row in outliers_df.iterrows():
        ax.scatter(
            row[variable_col],
            row['SYMBOL'],
            color=vaf_palette[row['VAF_category']],
            marker='X',
            s=150,
            edgecolor='black',
            linewidth=0.5
        )
        texts.append(
            ax.annotate(
                row['HGVS_label'],  # Use combined HGVSp if available else HGVSc
                xy=(row[variable_col], row['SYMBOL']),
                xytext=(row[variable_col] + 0.2, row['SYMBOL']),
                textcoords='data',
                fontsize=8,
                color='black',
                ha='left',
                va='center',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8, edgecolor='gray'),
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.6)
            )
        )

    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

    # Optional threshold line
    if threshold_value is not None:
        ax.axvline(threshold_value, color='red', linestyle='--', label=f'{variable_label} Threshold ({threshold_value})')

    ax.set_xlabel(variable_label)
    ax.set_ylabel('Gene (SYMBOL)')
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Outlier legend
    symbols_with_outliers = outliers_df['SYMBOL'].unique()
    handles = []
    for symbol in sorted(symbols_with_outliers):
        symbol_outliers = outliers_df[outliers_df['SYMBOL'] == symbol]
        handles.append(
            mlines.Line2D([], [], color='black', linestyle='None', markersize=10, label=f"{symbol}:")
        )
        for _, row in symbol_outliers.iterrows():
            handles.append(
                mlines.Line2D([], [], color=vaf_palette[row['VAF_category']], marker='X',
                              linestyle='None', markersize=5,
                              label=f"{row['HGVS_label']}, {variable_label}:({row[variable_col]:.2f})")
            )

    outlier_legend = ax.legend(
        handles=handles,
        title='Outlier Variants (|z| > 3)',
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        ncol=1,
        fontsize=10
    )

    vaf_legend_handles = [
        mlines.Line2D([], [], color=color, marker='o', linestyle='None', markersize=8, label=label)
        for label, color in vaf_palette.items()
    ]
    ax.add_artist(outlier_legend)
    ax.legend(
        handles=vaf_legend_handles,
        title='Likely VAF Category',
        bbox_to_anchor=(1, 0.35),
        loc='upper left',
        fontsize=10
    )

    plt.subplots_adjust(left=0.1, right=0.55, top=0.95, bottom=0.1)
    plt.tight_layout()
    plt.show()


plot_variant_outliers(
    pathogenic_df=pathogenic_df,
    non_pathogenic_df=non_pathogenic_df,
    variable_col='C-reactive protein',
    variable_label='CRP',
    threshold_value=2.59
)

plot_variant_outliers(
    pathogenic_df=pathogenic_df,
    non_pathogenic_df=non_pathogenic_df,
    variable_col='Rheumatoid factor',
    variable_label='RF',
    threshold_value=24.54
)

plot_variant_outliers(
    pathogenic_df=pathogenic_df,
    non_pathogenic_df=non_pathogenic_df,
    variable_col='White blood cell (leukocyte) count',
    variable_label='WBC',
    threshold_value=6.89
)

plot_variant_outliers(
    pathogenic_df=pathogenic_df,
    non_pathogenic_df=non_pathogenic_df,
    variable_col='Monocyte count',
    variable_label='MC',
    threshold_value=0.48
)

plot_variant_outliers(
    pathogenic_df=pathogenic_df,
    non_pathogenic_df=non_pathogenic_df,
    variable_col='Neutrophill count',
    variable_label='NC',
    threshold_value=4.23
)

plot_variant_outliers(
    pathogenic_df=pathogenic_df,
    non_pathogenic_df=non_pathogenic_df,
    variable_col='Lymphocyte count',
    variable_label='LC',
    threshold_value=1.97
)

plot_variant_outliers(
    pathogenic_df=pathogenic_df,
    non_pathogenic_df=non_pathogenic_df,
    variable_col='Red blood cell (erythrocyte) count',
    variable_label='RBC',
    threshold_value=4.52
)



#-----------------------------------------------------------------------------------------------------

def plot_crp_vs_vaf(pathogenic_df, non_pathogenic_df, gene_symbol):
    # Combine dataframes and label source
    pathogenic_df = pathogenic_df.copy()
    non_pathogenic_df = non_pathogenic_df.copy()
    pathogenic_df['source'] = 'pathogenic'
    non_pathogenic_df['source'] = 'non_pathogenic'
    ukb_df = pd.concat([pathogenic_df, non_pathogenic_df], ignore_index=True)

    # Filter for the given gene and drop missing data
    gene_df = ukb_df[ukb_df['SYMBOL'] == gene_symbol].copy()
    gene_df = gene_df.dropna(subset=['VAF', 'C-reactive protein', 'Sex', 'HGVSc'])

    # Calculate CRP z-score within each sex
    gene_df['CRP_zscore'] = gene_df.groupby('Sex')['C-reactive protein'].transform(zscore)

    # Identify outliers
    gene_df['is_outlier'] = gene_df['CRP_zscore'].abs() > 3

    # Extract short HGVS annotations
    gene_df['HGVSc_short'] = gene_df['HGVSc'].apply(lambda x: x.split(':')[-1])
    gene_df['HGVSp_short'] = gene_df['HGVSp'].apply(lambda x: x.split(':')[-1] if pd.notnull(x) else None)

    # Use HGVSp_short if available, else fallback to HGVSc_short
    gene_df['HGVS_label'] = gene_df.apply(
        lambda row: row['HGVSp_short'] if pd.notnull(row['HGVSp_short']) else row['HGVSc_short'],
        axis=1
    )

    sex_colors = {'Male': '#1f77b4', 'Female': '#ff7f0e'}

    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=300, sharey=True)

    for i, sex in enumerate(['Male', 'Female']):
        ax = axes[i]

        # Add shaded regions with labels only on first subplot
        if i == 0:
            ax.axvspan(0, 0.3, facecolor='lightblue', alpha=0.3, label='Likely Mosaic (VAF < 0.3)')
            ax.axvspan(0.3, 0.7, facecolor='lightgreen', alpha=0.3, label='Likely Germline Het (0.3 ≤ VAF < 0.7)')
            ax.axvspan(0.7, 1, facecolor='lightcoral', alpha=0.3, label='Likely Germline Homo (VAF ≥ 0.7)')
        else:
            ax.axvspan(0, 0.3, facecolor='lightblue', alpha=0.3)
            ax.axvspan(0.3, 0.7, facecolor='lightgreen', alpha=0.3)
            ax.axvspan(0.7, 1, facecolor='lightcoral', alpha=0.3)

        subset = gene_df[gene_df['Sex'] == sex]

        sns.scatterplot(
            data=subset,
            x='VAF',
            y='C-reactive protein',
            color=sex_colors[sex],
            edgecolor='black',
            s=80,
            ax=ax,
            alpha=0.7,
            legend=False
        )

        outliers = subset[subset['is_outlier']]
        ax.scatter(
            outliers['VAF'],
            outliers['C-reactive protein'],
            color='red',
            edgecolor='black',
            marker='X',
            s=120,
            label='Outlier (|z| > 3)' if i == 0 else None
        )

        # Annotate outliers with adjusted text to avoid overlap
        texts = []
        for _, row in outliers.iterrows():
            text = ax.text(
                row['VAF'],
                row['C-reactive protein'] + 0.05,
                f"{row['HGVS_label']}, z={row['CRP_zscore']:.1f}",
                fontsize=7,
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8)
            )
            texts.append(text)

        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle='-', color='black', lw=0.8)
        )
        ax.set_title(f'{gene_symbol}: {sex}', fontsize=13, fontweight='bold')
        ax.set_xlabel('Allele Frequency (AF)', fontsize=11)
        if i == 0:
            ax.set_ylabel('C-reactive Protein (CRP)', fontsize=11)
        else:
            ax.set_ylabel('')
        ax.grid(True)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(-0.02, 1), fontsize=9)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

plot_crp_vs_vaf(pathogenic_df, non_pathogenic_df, gene_symbol='UBA1')




def plot_crp_vs_vaf_by_clinsig(pathogenic_df, non_pathogenic_df, gene_symbol):
    # Combine and tag data sources
    pathogenic_df = pathogenic_df.copy()
    non_pathogenic_df = non_pathogenic_df.copy()
    pathogenic_df['source'] = 'pathogenic'
    non_pathogenic_df['source'] = 'non_pathogenic'
    df = pd.concat([pathogenic_df, non_pathogenic_df], ignore_index=True)

    # Filter for the gene of interest
    df = df[df['SYMBOL'] == gene_symbol].copy()
    df = df.dropna(subset=['VAF', 'C-reactive protein', 'Sex', 'HGVSc', 'CLIN_SIG'])

    # Extract first CLIN_SIG term
    df['CLIN_SIG_FIRST'] = df['CLIN_SIG'].apply(lambda x: x.split('&')[0].lower().strip())

    # Map CLIN_SIG_FIRST to clinical categories
    def map_clinsig(sig):
        if 'benign' in sig:
            return 'B/LB'
        elif 'uncertain' in sig:
            return 'VUS'
        elif 'pathogenic' in sig:
            return 'P/LP'
        else:
            return None

    df['CLIN_CATEGORY'] = df['CLIN_SIG_FIRST'].apply(map_clinsig)
    df = df.dropna(subset=['CLIN_CATEGORY'])

    # Calculate z-scores within each sex group
    df['CRP_zscore'] = df.groupby('Sex')['C-reactive protein'].transform(zscore)
    df['is_outlier'] = df['CRP_zscore'].abs() > 3

    # Extract short HGVS annotations
    df['HGVSc_short'] = df['HGVSc'].apply(lambda x: x.split(':')[-1])
    df['HGVSp_short'] = df['HGVSp'].apply(lambda x: x.split(':')[-1] if pd.notnull(x) else None)

    # Use HGVSp_short if available, otherwise fall back to HGVSc_short
    df['HGVS_label'] = df.apply(
        lambda row: row['HGVSp_short'] if pd.notnull(row['HGVSp_short']) else row['HGVSc_short'],
        axis=1
    )

    # Set up plot
    categories = ['B/LB', 'VUS', 'P/LP']
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=300, sharey=True)

    for i, category in enumerate(categories):
        ax = axes[i]
        subset = df[df['CLIN_CATEGORY'] == category]

        # Add shaded regions
        ax.axvspan(0, 0.3, facecolor='lightblue', alpha=0.3, label='VAF < 0.3')
        ax.axvspan(0.3, 0.7, facecolor='lightgreen', alpha=0.3, label='0.3 ≤ VAF < 0.7')
        ax.axvspan(0.7, 1, facecolor='lightcoral', alpha=0.3, label='VAF ≥ 0.7')

        # Plot all points
        sns.scatterplot(
            data=subset,
            x='VAF',
            y='C-reactive protein',
            hue='Sex',
            palette={'Male': '#1f77b4', 'Female': '#ff7f0e'},
            edgecolor='black',
            s=80,
            ax=ax,
            alpha=0.7,
            legend=False
        )

        # Highlight outliers
        outliers = subset[subset['is_outlier']]
        ax.scatter(
            outliers['VAF'],
            outliers['C-reactive protein'],
            color='red',
            edgecolor='black',
            marker='X',
            s=120,
            label='Outlier (|z| > 3)'
        )

        # Annotate outliers with adjusted text to avoid overlap
        texts = []
        for _, row in outliers.iterrows():
            text = ax.text(
                row['VAF'],
                row['C-reactive protein'] + 0.05,
                f"{row['HGVS_label']}, z={row['CRP_zscore']:.1f}",
                fontsize=7,
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8)
            )
            texts.append(text)

        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle='-', color='black', lw=0.8)
        )

        ax.set_title(f"{gene_symbol}: {category}", fontsize=13, fontweight='bold')
        ax.set_xlabel('Allele Frequency (AF)', fontsize=11)
        if i == 0:
            ax.set_ylabel('C-reactive Protein (CRP)', fontsize=11)
        else:
            ax.set_ylabel('')
        ax.grid(True)

    # Combine legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.01, 1), fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


plot_crp_vs_vaf_by_clinsig(pathogenic_df, non_pathogenic_df, gene_symbol='NLRP3')
