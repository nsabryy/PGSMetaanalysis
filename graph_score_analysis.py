import pandas as pd
import matplotlib.pyplot as plt

# Load the combined data
combined_file = '/N/project/compgen/PGSCalc/scoring_results/combined_pgs_scores.tsv'
df = pd.read_csv(combined_file, sep='\t')

df['alt_effect'] = df['total_variants'] - df['multiallelic'] - df['haplotype'] - df['ref_effect']
df['total_effect'] = df['ref_effect'] + df['alt_effect']

df['ref_effect_percentage'] = df['ref_effect'] / df['total_effect']
df['alt_effect_percentage'] = df['alt_effect'] / df['total_effect']

df_sorted = df.sort_values(by='ref_effect_percentage', ascending=False)

plt.figure(figsize=(10, 6))
plt.bar(range(len(df_sorted)), df_sorted['ref_effect_percentage'], width=1.0, label='Ref Effect')
plt.bar(range(len(df_sorted)), df_sorted['alt_effect_percentage'], bottom=df_sorted['ref_effect_percentage'], width=1.0, label='Alt Effect')

plt.title('Proportion of Ref and Alt Effects (Stacked to 1)')
plt.xlabel('PGS Models (sorted)')
plt.ylabel('Proportion (0-1)')
plt.legend(loc='upper right')
plt.tight_layout()

plt.savefig('/N/project/compgen/PGSCalc/scoring_results/ref_alt_proportion_stacked_total_1.png')

print("Stacked proportion plot with total set to 1 saved successfully.")
