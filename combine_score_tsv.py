import os
import pandas as pd

directory = '/N/project/compgen/PGSCalc/scoring_results/score_analysis'

combined_df = pd.DataFrame()

for root, dirs, files in os.walk(directory):
    if 'archive' in root:
        continue
    
    for file in files:
        if file.endswith('.tsv'):
            file_path = os.path.join(root, file)
            df = pd.read_csv(file_path, sep='\t')
            combined_df = pd.concat([combined_df, df])

combined_output_file = '/N/project/compgen/PGSCalc/scoring_results/combined_pgs_scores.tsv'
combined_df.to_csv(combined_output_file, sep='\t', index=False)

print(f"Combined files saved to {combined_output_file}")