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


# Get the IUIS list and get only AD or XL variants based off Inheritance 
iuis_iei = pd.read_excel("path/to/IUIS_OMIM_220628.xlsx")
AD_iuis = iuis_iei[iuis_iei['Inheritance.iuis'].str.contains('AD|XL', na=False)]
iuis_iei = AD_iuis
iuis_iei['OMIM'] = iuis_iei["OMIM_Phenotype"].fillna(iuis_iei["OMIM_gene"])
# Convert 'OMIM' column to numeric, forcing errors to NaN
iuis_iei['OMIM'] = pd.to_numeric(iuis_iei['OMIM'], errors='coerce')
# Drop rows where 'OMIM' is NaN or empty
iuis_iei = iuis_iei.dropna(subset=['OMIM'])
# Reset index after dropping rows
iuis_iei = iuis_iei.reset_index(drop=True)


# Load the data
df1 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/10-20_UKB_mosaic_variants.tsv', sep='\t')
df2 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/21-30_UKB_mosaic_variants.tsv', sep='\t')
df3 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/31-40_UKB_mosaic_variants.tsv', sep='\t')
df4 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/41-50_UKB_mosaic_variants.tsv', sep='\t')
df5 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/51-60_UKB_mosaic_variants.tsv', sep='\t')
df6 = pd.read_csv('C:/Users/Ricky/Garvan/ukb_processing-nf/extra_UKB_mosaic_variants.tsv', sep='\t')


# Ensure that 'POS' is a float
df1['POS'] = pd.to_numeric(df1['POS'], errors='coerce')
df2['POS'] = pd.to_numeric(df2['POS'], errors='coerce')
df3['POS'] = pd.to_numeric(df3['POS'], errors='coerce')
df4['POS'] = pd.to_numeric(df4['POS'], errors='coerce')
df5['POS'] = pd.to_numeric(df5['POS'], errors='coerce')
df6['POS'] = pd.to_numeric(df6['POS'], errors='coerce')



# Get the filtered Clinvar file
clinvar = pd.read_csv("path/to/clinvar_2022.tsv", sep="\t")


# Function to extract unique OMIM IDs
def extract_unique_omims(info_str):
    # Use regex to find all OMIM IDs (numeric values after OMIM:)
    omims = re.findall(r'OMIM:(\d+)', info_str)
    # Return unique OMIM IDs
    return list(set(omims))


# Extract the OMIM codes from the clinvar dataframe
clinvar['OMIMS'] = clinvar['INFO'].apply(extract_unique_omims)
clinvar = clinvar.drop(['ID', 'QUAL', 'FILTER'], axis=1)
clinvar['CHROM'] = 'chr' + clinvar['CHROM'].astype(str)
check_dup_clinvar = clinvar
check_dup_clinvar['OMIMS'] = check_dup_clinvar['OMIMS'].apply(lambda x: tuple(x) if isinstance(x, list) else x)
check_dup_clinvar = check_dup_clinvar.drop_duplicates()

# Get the UBA1 whitelist variants
vexas_uba1_bed = pd.read_csv('path/to/vexas_variants.bed', sep='\t', header=None, names=['CHROM', 'START', 'END'])
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

def filter_annotations(updated_df, iuis_iei, uba1_regions):
    """
    Filters the DataFrame based on criteria:
    - gnomADe_AF < 0.001 or gnomADg_AF < 0.001
    - IMPACT in ['MODERATE', 'HIGH']
    - CLIN_SIG in ['pathogenic', 'likely_pathogenic']
    - Includes rows where SYMBOL is 'UBA1', regardless of other criteria
    - Inheritance should already be filtered as AD|XL prior to calling this function

    Parameters:
        updated_df (pd.DataFrame): Input DataFrame with extracted annotations.
        iuis_iei (pd.DataFrame): DataFrame containing inheritance data already filtered to AD|XL.

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
    clin_sig_filter = df['CLIN_SIG'].str.contains(r'\b(pathogenic|likely_pathogenic)\b', case=False, na=False)
    
    # Include rows where SYMBOL is 'UBA1', regardless of other criteria
    uba1_filter = (
        (df['SYMBOL'] == 'UBA1') & 
        df['CHROM'].str.contains('chrX') & 
        df['POS'].apply(lambda x: any(start <= x <= end for start, end in uba1_regions))
    )
   
    # Now directly use the 'iuis_iei' DataFrame to filter based on AD|XL genes
    ad_genes = iuis_iei['Gene']  # Already filtered to AD|XL

    # Filter the DataFrame to only include genes with AD inheritance (since it's already filtered)
    ad_inheritance_filter = df['SYMBOL'].isin(ad_genes)

    # Combine all conditions: gnomad, impact, clin_sig, AD inheritance, and gene list (UBA1 exception)
    conditions = gnomad_filter & impact_filter & clin_sig_filter & ad_inheritance_filter
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
        'CLIN_SIG': r'CLIN_SIG=([^|;]+)',
        'HGVSc': r'HGVSc=([^|;]+)',
        'HGVSp': r'HGVSp=([^|;]+)',
        'VARIANT_CLASS': r'VARIANT_CLASS=([^|;]+)'
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

updated_df1 = replace_info_with_annotations(df1)
updated_df2 = replace_info_with_annotations(df2)
updated_df3 = replace_info_with_annotations(df3)
updated_df4 = replace_info_with_annotations(df4)
updated_df5 = replace_info_with_annotations(df5)
updated_df6 = replace_info_with_annotations(df6)


pathogenic_df1 = filter_annotations(updated_df1, iuis_iei, uba1_regions)
pathogenic_df2 = filter_annotations(updated_df2, iuis_iei, uba1_regions)
pathogenic_df3 = filter_annotations(updated_df3, iuis_iei, uba1_regions)
pathogenic_df4 = filter_annotations(updated_df4, iuis_iei, uba1_regions)
pathogenic_df5 = filter_annotations(updated_df5, iuis_iei, uba1_regions)
pathogenic_df6 = filter_annotations(updated_df6, iuis_iei, uba1_regions)


full_pathogenic_df = pd.concat([pathogenic_df1, 
                           pathogenic_df2, 
                           pathogenic_df3, 
                           pathogenic_df4, 
                           pathogenic_df5, 
                           pathogenic_df6], ignore_index=True)


full_pathogenic_df.to_csv('full_pathogenic_df.csv', index=False)



iuis_info =  iuis_iei[["Gene", "Disease", "OMIM", "Associated_features"]]


# Merge using an outer join to capture all possible matches
merged_df = full_pathogenic_df.merge(clinvar, 
                        left_on=['CHROM', 'POS', 'REF', 'ALT'], 
                        right_on=['CHROM', 'POS', 'REF', 'ALT'], 
                        how='inner')  # Regular merge


pathogenic_df_explode = merged_df.explode('OMIMS')
pathogenic_df_explode['OMIMS'] = pd.to_numeric(pathogenic_df_explode['OMIMS'])


merged_df = pathogenic_df_explode.merge(iuis_info , left_on=['SYMBOL','OMIMS'], right_on=['Gene','OMIM'], how='inner')


# Drop 'Genetic defect' column since SYMBOL already exists
merged_df = merged_df.drop(columns=["Gene"])
merged_df = merged_df.drop(columns=["OMIMS"])
# Move 'Disease ' and 'Associated features' to the 12th and 13th columns
cols = list(merged_df.columns)
cols.insert(12, cols.pop(cols.index("Disease")))
cols.insert(13, cols.pop(cols.index("OMIM")))
cols.insert(14, cols.pop(cols.index("Associated_features")))
cols.insert(15, cols.pop(cols.index("INFO")))
merged_df = merged_df[cols]
pathogenic_df = merged_df


# First, group by the relevant columns
grouped = pathogenic_df.groupby(['Participant ID', 'CHROM', 'POS', 'REF', 'ALT'])
# For each group, concatenate the specified columns
def concat_unique_values(group):
    # Keep the ALT, CHROM, POS, REF columns unchanged
    group = group.copy()
    # Concatenate values for the other specified columns
    for column in ['Disease', 'OMIM', 'Associated_features']:
        # If the column is OMIM, convert to string before joining
        if column == 'OMIM':
            group[column] = '| '.join(group[column].dropna().astype(int).astype(str).unique())  # Convert OMIM to int and then to str
        else:
            group[column] = '| '.join(group[column].dropna().unique())  # Drop NaN and join unique values for others
    return group.iloc[0]  # Return the first row in the group after concatenating
# Apply the function to concatenate values in each group
concatenated_df = grouped.apply(concat_unique_values).reset_index(drop=True)
pathogenic_df = concatenated_df
dups = concatenated_df[concatenated_df.duplicated(subset=['Participant ID', 'CHROM', 'POS', 'REF', 'ALT'], keep=False)]


# Function to check if an indel is larger than 2bp
def is_large_indel(row):
    return max(len(row["REF"]), len(row["ALT"])) > 2

# Filter out large indels
filtered_pathogenic_df = pathogenic_df[~pathogenic_df.apply(is_large_indel, axis=1)]

# Manually Whitelist this AIRE variant
AD_AIRE_mask = (
    (filtered_pathogenic_df["SYMBOL"] != "AIRE") |  # Keep non-AIRE rows
    (
        (filtered_pathogenic_df["CHROM"] == "chr12") & 
        (filtered_pathogenic_df["POS"] == 44289686) & 
        (filtered_pathogenic_df["REF"] == "G") & 
        (filtered_pathogenic_df["ALT"] == "T")
    )
)

# Apply the filter
filtered_df = filtered_pathogenic_df[AD_AIRE_mask]

# Manually exlucde all ACD variants
filtered_pathogenic_df = filtered_df[filtered_df["SYMBOL"] != "ACD"]
pathogenic_df = filtered_pathogenic_df

# Extract the Clinrevstat column to filter for high confidence ClinVar Annotations
pathogenic_df['CLNREVSTAT'] = pathogenic_df["INFO"].str.extract(r'CLNREVSTAT=([^;]+)')
filtered_pathogenic_df = pathogenic_df[~pathogenic_df["CLNREVSTAT"].str.startswith("no", na=False)]
pathogenic_df = filtered_pathogenic_df


# Remove the column and reinsert at position 15 (16th column)
col = pathogenic_df.pop("CLNREVSTAT")
pathogenic_df.insert(16, "CLNREVSTAT", col)
# Split FORMAT_INFO by ":" and extract the AF value (3rd field, index 2)
pathogenic_df["VAF"] = pathogenic_df["FORMAT_INFO"].str.split(":").str[2]
# Insert as column 17 (index 16 because Python is 0-based)
vaf_col = pathogenic_df.pop("VAF")
pathogenic_df.insert(17, "VAF", vaf_col)
# Ensure the VAF column is numeric
pathogenic_df['VAF'] = pd.to_numeric(pathogenic_df['VAF'], errors='coerce')  # Coerce any errors to NaN
pathogenic_df.iloc[:, 10] = pathogenic_df['INFO'].str.extract(r'CLNSIG=([^;]+)')


# Define whitelist of genes to always keep
whitelist_genes = ['UBA1']

pathogenic_mask = pathogenic_df['CLIN_SIG'].str.contains(r'\b(pathogenic|likely_pathogenic)\b', case=False, na=False)
clnrevstat_mask = ~pathogenic_df["CLNREVSTAT"].isin([
    "criteria_provided,_single_submitter",
    "criteria_provided,_conflicting_classifications"
])

whitelist_mask = pathogenic_df['SYMBOL'].isin(whitelist_genes)

# Apply filters
allowed_pathogenic = pathogenic_mask & clnrevstat_mask
combined_mask = allowed_pathogenic | whitelist_mask

# Final filtered DataFrame
pathogenic_df_filtered = pathogenic_df[combined_mask]

col = pathogenic_df_filtered.pop("HGVSp")
pathogenic_df_filtered.insert(12, "HGVSp", col)

pathogenic_df = pathogenic_df_filtered

pathogenic_df.to_csv('pathogenic_df.csv', index=False)


# Define list of HGVSc variants to exclude for G6PD
G6PD_excluded_variants = [
    "ENST00000393562.10:c.1388G>A",
    "ENST00000393562.10:c.1003G>A",
    "ENST00000393562.10:c.968T>C"
]

# Apply the filter to remove those rows
pathogenic_df_filtered = pathogenic_df_filtered[~((pathogenic_df_filtered['SYMBOL'] == "G6PD") & (pathogenic_df_filtered['HGVSc'].isin(G6PD_excluded_variants)))]

CFH_excluded_variants = [
    "ENST00000367429.9:c.3572C>T"
]

# Apply the filter to remove those rows
pathogenic_df_filtered = pathogenic_df_filtered[~((pathogenic_df_filtered['SYMBOL'] == "CFH") & (pathogenic_df_filtered['HGVSc'].isin(CFH_excluded_variants)))]


# Add CARD14 as the whitelist
card14_df = full_pathogenic_df[
    (full_pathogenic_df['SYMBOL'] == 'CARD14') & 
    (full_pathogenic_df['POS'] == 80182790)
].copy()

# Merge with iuis on SYMBOL == Gene to bring in the additional columns
card14_df = card14_df.merge(iuis_info, left_on='SYMBOL', right_on='Gene', how='left')
# Identify any columns in pathogenic_df not present in card14_df
missing_cols = [col for col in pathogenic_df_filtered.columns if col not in card14_df.columns]

# Add those missing columns as empty (NaN)
for col in missing_cols:
    card14_df[col] = pd.NA

# Reorder card14_df columns to match pathogenic_df
card14_df = card14_df[pathogenic_df_filtered.columns]

# Append to pathogenic_df
pathogenic_df_filtered = pd.concat([pathogenic_df_filtered, card14_df], ignore_index=True)
pathogenic_df = filtered_pathogenic_df
pathogenic_df.to_csv('filtered_pathogenic_df.csv', index=False)


# Use the first 8 columns for filtering
common_columns = updated_df1.columns[:8]  # assuming both DataFrames have same first 8 columns


#Filter the updated_df DataFrames to only keep rows that are not in pathogenic_df_filtered
# Perform anti-join to keep only rows in updated_df1 NOT in pathogenic_df
non_pathogenic_df1 = updated_df1.merge(
    pathogenic_df[common_columns],
    on=list(common_columns),
    how='left',
    indicator=True
).query('_merge == "left_only"').drop(columns=['_merge'])
non_pathogenic_df2 = updated_df2.merge(
    pathogenic_df[common_columns],
    on=list(common_columns),
    how='left',
    indicator=True
).query('_merge == "left_only"').drop(columns=['_merge'])
non_pathogenic_df3 = updated_df3.merge(
    pathogenic_df[common_columns],
    on=list(common_columns),
    how='left',
    indicator=True
).query('_merge == "left_only"').drop(columns=['_merge'])
non_pathogenic_df4 = updated_df4.merge(
    pathogenic_df[common_columns],
    on=list(common_columns),
    how='left',
    indicator=True
).query('_merge == "left_only"').drop(columns=['_merge'])
non_pathogenic_df5 = updated_df5.merge(
    pathogenic_df[common_columns],
    on=list(common_columns),
    how='left',
    indicator=True
).query('_merge == "left_only"').drop(columns=['_merge'])
non_pathogenic_df6 = updated_df6.merge(
    pathogenic_df[common_columns],
    on=list(common_columns),
    how='left',
    indicator=True
).query('_merge == "left_only"').drop(columns=['_merge'])

# Group by 'Participant ID' and randomly select one row from each group
non_pathogenic_df_single_variant1 = non_pathogenic_df1.groupby('Participant ID').apply(lambda x: x.sample(n=1)).reset_index(drop=True)
non_pathogenic_df_single_variant2 = non_pathogenic_df2.groupby('Participant ID').apply(lambda x: x.sample(n=1)).reset_index(drop=True)
non_pathogenic_df_single_variant3 = non_pathogenic_df3.groupby('Participant ID').apply(lambda x: x.sample(n=1)).reset_index(drop=True)
non_pathogenic_df_single_variant4 = non_pathogenic_df4.groupby('Participant ID').apply(lambda x: x.sample(n=1)).reset_index(drop=True)
non_pathogenic_df_single_variant5 = non_pathogenic_df5.groupby('Participant ID').apply(lambda x: x.sample(n=1)).reset_index(drop=True)
non_pathogenic_df_single_variant6 = non_pathogenic_df6.groupby('Participant ID').apply(lambda x: x.sample(n=1)).reset_index(drop=True)


non_pathogenic_df = pd.concat([non_pathogenic_df_single_variant1, 
                           non_pathogenic_df_single_variant2, 
                           non_pathogenic_df_single_variant3, 
                           non_pathogenic_df_single_variant4, 
                           non_pathogenic_df_single_variant5, 
                           non_pathogenic_df_single_variant6], ignore_index=True)

non_pathogenic_df.to_csv('non_pathogenic_df.csv', index=False)


# Further filtering manual exlusion of G6PD and Females from X-linked variants

pathogenic_df = pd.read_csv("filtered_pathogenic_df.csv")
non_pathogenic_df = pd.read_csv("non_pathogenic_df.csv")
# Split FORMAT_INFO by ":" and extract the AF value (3rd field, index 2)
non_pathogenic_df["VAF"] = non_pathogenic_df["FORMAT_INFO"].str.split(":").str[2]
non_pathogenic_df['VAF'] = pd.to_numeric(non_pathogenic_df['VAF'], errors='coerce')  # Coerce any errors to NaN
# Step 1: Extract G6PD rows from pathogenic_df
g6pd_rows = pathogenic_df[pathogenic_df['SYMBOL'] == 'G6PD']
# Step 2: Append these rows to non_pathogenic_df
non_pathogenic_df = pd.concat([non_pathogenic_df, g6pd_rows], ignore_index=True)
# Optionally, ensure the column order is consistent (if required)
non_pathogenic_df = non_pathogenic_df[pathogenic_df.columns]
# Step 3: If you want to drop G6PD from pathogenic_df as well
pathogenic_df = pathogenic_df[pathogenic_df['SYMBOL'] != 'G6PD']


# Step 1: Extract rows on chrX where Sex is female from pathogenic_df
chrX_female_rows = pathogenic_df[(pathogenic_df['CHROM'] == 'chrX') & (pathogenic_df['Sex'].str.lower() == 'female')]
# Step 2: Append these rows to non_pathogenic_df
non_pathogenic_df = pd.concat([non_pathogenic_df, chrX_female_rows], ignore_index=True)
# Ensure the column order is consistent (optional)
non_pathogenic_df = non_pathogenic_df[pathogenic_df.columns]
# Step 3: Remove these rows from pathogenic_df
pathogenic_df = pathogenic_df[~((pathogenic_df['CHROM'] == 'chrX') & (pathogenic_df['Sex'].str.lower() == 'female'))]


pathogenic_df.to_csv('final_pathogenic_df.csv', index=False)
non_pathogenic_df.to_csv('non_pathogenic_df.csv', index=False)


# =============================================================================
# ANALYSIS OF THE DATA
# =============================================================================
pathogenic_df = pd.read_csv("final_pathogenic_df.csv")
non_pathogenic_df = pd.read_csv("non_pathogenic_df.csv")


