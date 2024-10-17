import pandas as pd
import matplotlib.pyplot as plt

# Load the combined data file
combined_file = '/N/project/compgen/PGSCalc/scoring_results/combined_pgs_scores.tsv'
df = pd.read_csv(combined_file, sep='\t')

# Calculate alt_effect as total_variants - multiallelic - haplotype - ref_effect
df['alt_effect'] = df['total_variants'] - df['multiallelic'] - df['haplotype'] - df['ref_effect']

# Calculate the total as ref_effect + alt_effect
df['total_effect'] = df['ref_effect'] + df['alt_effect']

# Calculate the proportion of ref_effect and alt_effect relative to their total
df['ref_effect_percentage'] = df['ref_effect'] / df['total_effect']
df['alt_effect_percentage'] = df['alt_effect'] / df['total_effect']

# Sort data based on ref_effect proportion (you can adjust sorting criteria if needed)
df_sorted = df.sort_values(by='ref_effect_percentage', ascending=False)

# Plotting stacked bar for ref and alt effect proportions, with the total going to 1
plt.figure(figsize=(10, 6))

# Plot stacked bar: ref_effect_percentage on bottom, alt_effect_percentage on top
plt.bar(range(len(df_sorted)), df_sorted['ref_effect_percentage'], width=1.0, label='Ref Effect')
plt.bar(range(len(df_sorted)), df_sorted['alt_effect_percentage'], bottom=df_sorted['ref_effect_percentage'], width=1.0, label='Alt Effect')

plt.title('Proportion of Ref and Alt Effects (Stacked to 1)')
plt.xlabel('PGS Models (sorted)')
plt.ylabel('Proportion (0-1)')
plt.legend(loc='upper right')
plt.tight_layout()

# Save the plot
plt.savefig('/N/project/compgen/PGSCalc/scoring_results/ref_alt_proportion_stacked_total_1.png')

print("Stacked proportion plot with total set to 1 saved successfully.")
