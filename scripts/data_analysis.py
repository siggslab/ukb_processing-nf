import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import re


GENE_LIST = {
    "STING1", "LSM11", "RNU7-1", "CDC42", "STAT2", 
    "ATAD3A", "C2orf69", "RIPK1", "NCKAP1L", 
    "HCK1", "PSMB9", "IKBKG", "TBK1", "UBA1", "ADA2",
    "TREX1", "SAMHD1", "IFIH1", "DNASE2", "LSM11", "ACP5",
    "POLA1", "USP18", "OAS1", "MEFV", "MVK", "NLRP3",
    "NLRP12", "NLRC4", "PLCG2", "NLRP1", "TNFRSF1A",
    "PSTPIP1", "NOD2", "ADAM17", "LPIN2", "IL1RN",
    "IL36RN","SLC29A3", "CARD14", "SH3BP2", "PSMB8",
    "PSMG2", "COPA", "OTULIN", "TNFAIP3", "AP1S3",
    "ALPI", "TRIM22", "HAVCR2", "SYK", "HCK"
}


# Load your data into a pandas DataFrame
df1 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/10_UKB_mosaic_variants.tsv', sep='\t')
df2 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/11-24_UKB_mosiac_variants.tsv', sep='\t')
df3 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/15-28_UKB_mosiac_variants.tsv', sep='\t')
# Concatenate the dataframes
merged_df = pd.concat([df1, df2, df3], ignore_index=True)

# Ensure that 'POS' is a float
merged_df['POS'] = pd.to_numeric(merged_df['POS'], errors='coerce')

iuis_iei = pd.read_excel('C:/Users/Ricky/Garvan/ukb_processing-nf/assets/IUIS IEI list for web site March 2022.xlsx')

vexas_uba1_bed = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/assets/vexas_variants.bed', sep='\t', header=None, names=['CHROM', 'START', 'END'])


# Extract UBA1 regions from the BED file (chrX with specified coordinates)
uba1_regions = vexas_uba1_bed[vexas_uba1_bed ['CHROM'] == 'chrX'][['START', 'END']].values.tolist()

def replace_info_with_annotations(df):
    """
    Replaces the INFO column with extracted annotation columns in the DataFrame.
    The new columns start at the position of the original INFO column, shifting others.

    Parameters:
        df (pd.DataFrame): Original DataFrame containing an INFO column.

    Returns:
        pd.DataFrame: Updated DataFrame with new annotation columns replacing INFO.
    """
    # Extract annotations from INFO column
    annotations = df['INFO'].apply(extract_annotations)
    annotations_df = pd.DataFrame(annotations.tolist())

    # Find the position of the INFO column
    info_index = df.columns.get_loc('INFO')

    # Drop the INFO column
    df = df.drop(columns=['INFO'])

    # Insert the new annotation columns at the INFO column's position
    for i, col in enumerate(annotations_df.columns):
        df.insert(info_index + i, col, annotations_df[col])

    return df

def filter_annotations(df, iuis_iei, uba1_regions):
    """
    Filters the DataFrame based on criteria:
    - gnomADe_AF < 0.001 or gnomADg_AF < 0.001
    - IMPACT in ['MODERATE', 'HIGH']
    - CLIN_SIG in ['pathogenic', 'likely_pathogenic']
    - Includes rows where SYMBOL is 'UBA1', regardless of other criteria
    - Only include genes that are associated with Autosomal Dominant (AD) inheritance

    Parameters:
        df (pd.DataFrame): Input DataFrame with extracted annotations.
        iuis_iei (pd.DataFrame): DataFrame containing inheritance data.

    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    # Define filter conditions for gnomAD, impact, and clinical significance
    df = updated_df
    gnomad_filter = (
        (df['gnomADe_AF'].notna() & (df['gnomADe_AF'] < 0.001)) |
        (df['gnomADe_AF'].isna() & df['gnomADg_AF'].notna() & (df['gnomADg_AF'] < 0.001))
    )
    impact_filter = df['IMPACT'].isin(['MODERATE', 'HIGH'])
    clin_sig_filter = df['CLIN_SIG'].isin(['pathogenic', 'likely_pathogenic'])
    
    # Include rows where SYMBOL is 'UBA1', regardless of other criteria
    uba1_filter = (
        (df['SYMBOL'] == 'UBA1') & 
        df['CHROM'].str.contains('chrX') & 
        df['POS'].apply(lambda x: any(start <= x <= end for start, end in uba1_regions))
    )

    #uba1_filter = df['SYMBOL'] == 'UBA1'
    # Gene list filter: genes in GENE_LIST must pass the filter, except UBA1
    gene_list_filter = df['SYMBOL'].isin(GENE_LIST)
    
    # Extract genes associated with Autosomal Dominant (AD) inheritance from iuis_iei DataFrame
    ad_genes = iuis_iei[iuis_iei['Inheritance'].str.contains('AD', case=False, na=False)]['Genetic defect']
    
    # Filter the DataFrame to only include genes with AD inheritance
    ad_inheritance_filter = df['SYMBOL'].isin(ad_genes)

    # Combine all conditions: gnomad, impact, clin_sig, AD inheritance, and gene list (UBA1 exception)
    conditions = gnomad_filter & impact_filter & clin_sig_filter & gene_list_filter & ad_inheritance_filter
    #conditions = gnomad_filter & impact_filter & clin_sig_filter & gene_list_filter
    combined_filter = conditions | uba1_filter  
    # In case no rows pass the uba1_filter, still retain the rows that pass the conditions filter
    if uba1_filter.sum() == 0:  # If no UBA1 rows pass the filter
        filtered_df = df[conditions]  # Keep rows that pass the other conditions
    else:
        # Apply the filter to the DataFrame
        filtered_df = df[combined_filter]
    
    return filtered_df

def extract_annotations(info):
    # Regex pattern to extract each required field from the INFO column
    pattern = {
        'Consequence': r'Consequence=([^|;]+)',
        'IMPACT': r'IMPACT=([^|;]+)',
        'SYMBOL': r'SYMBOL=([^|;]+)',
        'gnomADg_AF': r'gnomADg_AF=([^|;]+)',
        'gnomADe_AF': r'gnomADe_AF=([^|;]+)',
        'CLIN_SIG': r'CLIN_SIG=([^|;]+)'
    }
    
    # Extract each field using regex
    annotations = {}
    for key, pat in pattern.items():
        match = re.search(pat, info)
        if match:
            value = match.group(1)
            # Convert gnomADg_AF and gnomADe_AF to floats if needed
            if key in ['gnomADg_AF', 'gnomADe_AF']:
                value = value.replace('_', '-')  # Replace underscore with dash
                try:
                    value = float(value)  # Convert to float
                except ValueError:
                    value = None  # Handle invalid float values
            annotations[key] = value
        else:
            annotations[key] = None  # If no match, store None
    
    return annotations

updated_df = replace_info_with_annotations(merged_df)
pathogenic_df = filter_annotations(updated_df, iuis_iei, uba1_regions)

pathogenic_df.to_csv('pathogenic_df.csv', index=False)
 
# Remove rows in `filtered_df` from `updated_df`
non_pathogenic_df = pd.concat([updated_df, pathogenic_df, pathogenic_df]).drop_duplicates(keep=False)

# Group by 'Participant ID' and randomly select one row from each group
non_pathogenic_df_single_variant = non_pathogenic_df.groupby('Participant ID').apply(lambda x: x.sample(n=1)).reset_index(drop=True)

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



# Create an empty dictionary to store the results
shapiro_results = {}

# Loop through each biomarker and apply the Shapiro-Wilk test for both datasets
for biomarker in biomarker_columns:
    # Drop missing values for both datasets
    pathogenic_data = pathogenic_df[biomarker].dropna()
    non_pathogenic_data = non_pathogenic_df[biomarker].dropna()
    
    # Perform Shapiro-Wilk test on both datasets
    pathogenic_stat, pathogenic_p_value = stats.shapiro(pathogenic_data)
    non_pathogenic_stat, non_pathogenic_p_value = stats.shapiro(non_pathogenic_data)
    
    # Store results for both datasets
    shapiro_results[biomarker] = {
        'Pathogenic_Statistic': pathogenic_stat, 
        'Pathogenic_p-value': pathogenic_p_value,
        'Non-Pathogenic_Statistic': non_pathogenic_stat, 
        'Non-Pathogenic_p-value': non_pathogenic_p_value
    }

# Convert the results into a DataFrame for easier viewing
shapiro_results_df = pd.DataFrame(shapiro_results).T

# Add conclusions to the results for both datasets
shapiro_results_df['Pathogenic_Conclusion'] = shapiro_results_df['Pathogenic_p-value'].apply(
    lambda x: 'Normally distributed' if x > 0.05 else 'Not normally distributed'
)

shapiro_results_df['Non-Pathogenic_Conclusion'] = shapiro_results_df['Non-Pathogenic_p-value'].apply(
    lambda x: 'Normally distributed' if x > 0.05 else 'Not normally distributed'
)

# Separate biomarkers into normally distributed and non-normally distributed for both datasets
normal_biomarkers_pathogenic = shapiro_results_df[shapiro_results_df['Pathogenic_Conclusion'] == 'Normally distributed'].index
non_normal_biomarkers_pathogenic = shapiro_results_df[shapiro_results_df['Pathogenic_Conclusion'] == 'Not normally distributed'].index

normal_biomarkers_non_pathogenic = shapiro_results_df[shapiro_results_df['Non-Pathogenic_Conclusion'] == 'Normally distributed'].index
non_normal_biomarkers_non_pathogenic = shapiro_results_df[shapiro_results_df['Non-Pathogenic_Conclusion'] == 'Not normally distributed'].index

# Create an empty dictionary to store p-values
comparison_results = {}
# Loop through each biomarker to perform the appropriate test
for biomarker in biomarker_columns:
    # Extract data for pathogenic and non-pathogenic
    pathogenic_data = pathogenic_df[biomarker].dropna()
    non_pathogenic_data = non_pathogenic_df_single_variant[biomarker].dropna()
    
    pathogenic_averages = pathogenic_df[biomarker].mean()
    non_pathogenic_averages = non_pathogenic_df_single_variant[biomarker].mean()

    # Perform Mann-Whitney U test if the biomarker is not normally distributed
    stat, p_value = stats.mannwhitneyu(pathogenic_data, non_pathogenic_data)
    
    # Store the results
    comparison_results[biomarker] = {'p-value': p_value, 'Pathogenic': pathogenic_averages, 'Non_Pathogenic': non_pathogenic_averages}

# Convert the results into a DataFrame for easy viewing
comparison_results_df = pd.DataFrame(comparison_results).T
print(comparison_results_df)


# Create an empty dictionary to store p-values
comparison_results = {}
# Loop through each biomarker to perform the appropriate test
for biomarker in biomarker_columns:
    # Extract data for pathogenic and non-pathogenic
    pathogenic_data = pathogenic_df[biomarker].dropna()
    non_pathogenic_data = non_pathogenic_df[biomarker].dropna()
    
    pathogenic_averages = pathogenic_df[biomarker].mean()
    non_pathogenic_averages = non_pathogenic_df[biomarker].median()
    
    # Perform Mann-Whitney U test if the biomarker is not normally distributed
    stat, p_value = stats.mannwhitneyu(pathogenic_data, non_pathogenic_data)
    
    # Store the results
    comparison_results[biomarker] = {'p-value': p_value, 'Pathogenic': pathogenic_averages, 'Non_Pathogenic': non_pathogenic_averages}

# Convert the results into a DataFrame for easy viewing
comparison_results_df = pd.DataFrame(comparison_results).T
print(comparison_results_df)





# Ensure 'Age at recruitment' is numeric and handle missing values
pathogenic_df['Age at recruitment'] = pd.to_numeric(pathogenic_df['Age at recruitment'], errors='coerce')
non_pathogenic_df['Age at recruitment'] = pd.to_numeric(non_pathogenic_df['Age at recruitment'], errors='coerce')

# Remove rows with missing age values
pathogenic_df = pathogenic_df.dropna(subset=['Age at recruitment'])
non_pathogenic_df = non_pathogenic_df.dropna(subset=['Age at recruitment'])

# Create a combined DataFrame for easier plotting
pathogenic_df['Group'] = 'Pathogenic'
non_pathogenic_df['Group'] = 'Non-Pathogenic'
combined_df = pd.concat([pathogenic_df[['Age at recruitment', 'Group']], non_pathogenic_df[['Age at recruitment', 'Group']]])

# Plot age distribution using seaborn (boxplot)
plt.figure(figsize=(10, 6))
sns.boxplot(data=combined_df, x='Group', y='Age at recruitment', palette='Set1')
plt.title('Age Distribution Comparison: Pathogenic vs Non-Pathogenic', fontsize=15)
plt.xlabel('Group', fontsize=12)
plt.ylabel('Age at Recruitment', fontsize=12)
plt.show()

# Optionally, display descriptive statistics (mean, median, etc.)
pathogenic_stats = pathogenic_df['Age at recruitment'].describe()
non_pathogenic_stats = non_pathogenic_df['Age at recruitment'].describe()

print("Pathogenic Group Age Distribution:\n", pathogenic_stats)
print("\nNon-Pathogenic Group Age Distribution:\n", non_pathogenic_stats)




# Create a new DataFrame to count the occurrences of each consequence per gene
consequence_counts = pathogenic_df.groupby(['SYMBOL', 'Consequence']).size().unstack(fill_value=0)

# Calculate total counts for each gene
gene_totals = consequence_counts.sum(axis=1)

# Get the top 10 genes by total count (descending order)
top_10_genes = gene_totals.sort_values(ascending=False).head(10).index

# Filter consequence counts to include only the top 10 genes
consequence_counts = consequence_counts.loc[top_10_genes]

# Calculate the total number of variants in the dataset
total_variants = gene_totals.sum()

# Get the unique consequences that appear in the top 10 genes
top_consequences = consequence_counts.columns[consequence_counts.sum(axis=0) > 0]

# Plotting the stacked bar chart with higher resolution and wider bars for more space
fig, ax = plt.subplots(figsize=(12, 8), dpi=300)  # Increase DPI for higher quality
consequence_counts[top_consequences].plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.6)  # Adjust bar width here

# Annotating each bar with the percentage distribution of the gene in the dataset
for i, gene in enumerate(consequence_counts.index):
    total = gene_totals[gene]
    percentage = (total / total_variants) * 100
    # Get the height of the stacked bars for the gene
    height = gene_totals[gene]
    ax.text(i, height + 0.05 * height, 
            f'{percentage:.1f}%\n({total})', 
            ha='center', va='bottom', fontsize=8, color='black')

# Customize plot
plt.title('Distribution of Variant Consequence Across LP/P Variants')
plt.xlabel('Gene')
plt.ylabel('Count of Variants')
plt.xticks(rotation=90)  # Rotate x-axis labels for readability
plt.legend(title='Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')

# Adjust y-axis to add more space at the top
max_y = gene_totals[top_10_genes].max()
ax.set_ylim(0, max_y * 1.15)  # Add 15% space to the top of the y-axis
# Increase space between the y-axis and first stacked bar
plt.subplots_adjust(left=0.1)  # Adjust left margin here (increase the value to add more space)
# Increase space at the top for annotations
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust right margin for legend and space at the top

# Show plot
plt.show()




# Group by Gene and Ethnic background, then count occurrences
ancestry_counts = pathogenic_df.groupby(['SYMBOL', 'Ethnic background']).size().unstack(fill_value=0)

# Sort by the total count of variants across all ethnic backgrounds for each gene (descending order)
top_10_genes = ancestry_counts.sum(axis=1).sort_values(ascending=False).head(10).index
ancestry_counts = ancestry_counts.loc[top_10_genes]

# Normalize by dividing by the sum of each row to get percentages
ancestry_percentages = ancestry_counts.div(ancestry_counts.sum(axis=1), axis=0) * 100

# Get the unique ethnic backgrounds that appear in the top 10 genes
top_ethnic_backgrounds = ancestry_percentages.columns[ancestry_percentages.sum(axis=0) > 0]

# Plotting the horizontal bar chart
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)  # Increase DPI for better image quality
ancestry_percentages[top_ethnic_backgrounds].plot(kind='barh', stacked=True, ax=ax, colormap='tab20')

# Customize plot
plt.title('Ancestry Breakdown of Gene Variants (Top 10 Genes)')
plt.xlabel('Percentage of Variants')
plt.ylabel('Gene')
plt.legend(title='Ethnic Background', bbox_to_anchor=(1.05, 1), loc='upper left')

# Increase space at the top for readability
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Show plot
plt.show()




# Function to compare allele frequencies using Mann-Whitney U test
def compare_allele_frequencies_mannwhitney(pathogenic_df, non_pathogenic_df, column):
    path_values = pathogenic_df[column].dropna()
    non_path_values = non_pathogenic_df[column].dropna()

    if len(path_values) < 2:
        print(f"Not enough data in {column} for meaningful comparison.")
        return

    stat, p_value = stats.mannwhitneyu(path_values, non_path_values, alternative="two-sided")

    print(f"Comparison for {column}:")
    print(f"  - Median (Pathogenic): {path_values.mean():.6f}")
    print(f"  - Median (Non-Pathogenic): {non_path_values.mean():.6f}")
    print(f"  - Mann-Whitney U statistic: {stat:.3f}, p-value: {p_value:.6f}")
    print("-" * 50)

# Run the test for both allele frequency columns
compare_allele_frequencies_mannwhitney(pathogenic_df, non_pathogenic_df_single_variant, 'gnomADg_AF')
compare_allele_frequencies_mannwhitney(pathogenic_df, non_pathogenic_df_single_variant, 'gnomADe_AF')

