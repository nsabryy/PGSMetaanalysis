import subprocess
import os
import pandas as pd
import logging
import pickle
import hashlib
import requests 
import pysam
import argparse
import csv
# # Get the logger instance
# logger = logging.getLogger(__name__)

class ScorePackage:
    def __init__(self, pgs_ids=None, trait_id=None, trait_include_children=False, genome_build='GRCh38', scoring_file_directory=None, override=True):
        logging.info("Initializing score package.")
        logging.info(f"trait: {trait_id}")
        self.output_dir = "/N/project/compgen/PGSCalc/src/PGScore/scoring_files"
        self.scores = []
        self.combined_scores = {}

        if scoring_file_directory:
            self._load_local_scores(scoring_file_directory)
        else:
            if pgs_ids:
                subset_hash = hashlib.md5(','.join(sorted(pgs_ids)).encode()).hexdigest()[:8]
            elif trait_id:
                subset_hash = hashlib.md5(trait_id.encode()).hexdigest()[:8]

            if subset_hash:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                save_path = os.path.join(script_dir, "../PGScore/score_package_objects", f"{subset_hash}_{genome_build}_ref.pkl")
                save_path = os.path.abspath(save_path)
            else:
                logging.warning("Unable to generate unique ID for subset. ScorePackage not cacheable.")

            if not override and os.path.exists(save_path):
                self.load(save_path)
                logging.info("Using previously saved ScorePackage object.")
                
            else:
                self.pgs_ids = pgs_ids # Don't know this if using PGS ID
                self.trait_id = trait_id
                self.trait_include_children = trait_include_children
                self.genome_build = genome_build
                self._download_scores(override)
                self._combine_scores()
                if save_path:
                    self.save(save_path)
            logging.info(f"Scoring object built for ids: {self.pgs_ids}")

    def save(self, save_path):
        with open(save_path, 'wb') as f:
            pickle.dump(self, f)
        logging.info(f"VCFQuery object saved to {save_path}.")

    def load(self, save_path):
        with open(save_path, 'rb') as f:
            obj = pickle.load(f)
        self.__dict__.update(obj.__dict__)
        logging.info(f"VCFQuery object loaded from {save_path}.")

    def _download_scores(self, override):
        logging.info("Starting the download of PGS scores.")

        if self.trait_id:
            logging.info(f"Pulling associated traits with {self.trait_id}")
            base_url = "https://www.pgscatalog.org/rest/trait/"
            url = f"{base_url}{self.trait_id}?include_children={self.trait_include_children}"
            response = requests.get(url)
            
            if response.status_code == 200:
                data = response.json()
                pgs_ids = data.get("associated_pgs_ids", [])
                logging.info(f"Received PGS IDs: {pgs_ids}")
                # TODO if we want to change so you can do a trait and input pgs ids, change here!
                self.pgs_ids = list(pgs_ids)
            
            command = f"python -m pgscatalog_utils.download.download_scorefile -t {self.trait_id} -b {self.genome_build} -o {self.output_dir}"
        else:
            ids_str = ' '.join(self.pgs_ids)
            command = f"python -m pgscatalog_utils.download.download_scorefile -i {ids_str} -o {self.output_dir} -b {self.genome_build}"
        # Run the command
        try:
            subprocess.run(command, shell=True, check=True)
            logging.info("Download completed successfully.")
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while downloading scores: {e}")
            return

        # Create PGScore objects
        for pgs_id in self.pgs_ids:
            # TODO double check this - harmonized is getting fetched
            scoring_file_path = os.path.join(self.output_dir, f"{pgs_id}_hmPOS_{self.genome_build}.txt.gz")
            if os.path.exists(scoring_file_path):
                logging.info(f"Creating PGScore object for {pgs_id}.")
                try:
                    self.scores.append(PGScore(pgs_id, scoring_file_path, self.genome_build, override))
                except Exception as e:
                    logging.error(f"Error creating PGScore object for {pgs_id}: {e}")
            else:
                logging.warning(f"Warning: {scoring_file_path} does not exist.")

    def _load_local_scores(self, directory):
        logging.info("Loading local scorefiles.")
        for filename in os.listdir(directory):
            if filename.endswith('.txt') or filename.endswith('.txt.gz'):
                pgs_id = filename.split('_')[0]
                scoring_file_path = os.path.join(directory, filename)
                if os.path.exists(scoring_file_path):
                    logging.info(f"Creating PGScore object for {pgs_id}.")
                    self.scores.append(PGScore(pgs_id, scoring_file_path))
                else:
                    logging.warning(f"Warning: {scoring_file_path} does not exist.")

    def _combine_scores(self):
        for score in self.scores:
            snp = score.snps
            # print(f"adding {snp}")
            for chrom, pos_dict in snp.items():
                if chrom not in self.combined_scores.keys():
                    self.combined_scores[chrom] = {}
                for pos, scoring_info in pos_dict.items():
                    if pos not in self.combined_scores[chrom].keys():
                        self.combined_scores[chrom][pos] = []
                    self.combined_scores[chrom][pos].extend(snp[chrom][pos])

class PGScore:
    def __init__(self, pgs_id, scoring_file_path, build="GRCh38", check_ref=None, override=True):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        save_path = os.path.join(script_dir, "../PGScore/score_objects", f"{pgs_id}_{build}_ref.pkl")
        save_path = os.path.abspath(save_path)

        if not override and os.path.exists(save_path):
            self.load(save_path)
            logging.info(f"Using previously saved PGScore object for {pgs_id} in {build}.")
        else:
            self.pgs_id = pgs_id
            self.genome_build = build
            self.scoring_file_path = scoring_file_path
            self.total_snps = self._extract_total_snps()
            self.multiallelic = 0
            self.haplotype = 0
            self.ambiguous_strand = 0
            self.ref_effect = 0
            self.data = self._parse_scoring()
            self.snps = self._extract_snps(check_ref)
            self.save(save_path)

        def _extract_total_snps(self):
            """Extracts the total number of SNPs from the scoring file header."""
            logging.info(f"Extracting total SNP count from the header of {self.scoring_file_path}.")
            try:
                if self.scoring_file_path.endswith('.gz'):
                    with gzip.open(self.scoring_file_path, 'rt') as file:
                        for line in file:
                            if line.startswith("#variants_number="):
                                return int(line.strip().split('=')[1])
                else:
                    with open(self.scoring_file_path, 'r') as file:
                        for line in file:
                            if line.startswith("#variants_number="):
                                return int(line.strip().split('=')[1])
            except Exception as e:
                logging.error(f"Error reading total SNPs from header in {self.scoring_file_path}: {e}")
                raise
            return None

    def _parse_scoring(self):
        logging.info(f"Parsing scoring file for {self.pgs_id}.")
        if not os.path.exists(self.scoring_file_path):
            logging.info(f"Scoring file {self.scoring_file_path} not found. Downloading...")
            self.download_score()
        try:
            if self.scoring_file_path.endswith('.gz'):
                data = pd.read_csv(self.scoring_file_path, compression='gzip', sep='\t', comment='#')
            else:
                data = pd.read_csv(self.scoring_file_path, sep='\t', comment='#')
        except Exception as e:
            logging.error(f"Error reading scoring file {self.scoring_file_path}: {e}")
            raise
        return data

    def download_score(self):
        ids_str = self.pgs_id
        output_dir = os.path.dirname(self.scoring_file_path)
        command = f"python -m pgscatalog_utils.download.download_scorefile -i {ids_str} -o {output_dir} -b {self.genome_build}"
        try:
            subprocess.run(command, shell=True, check=True)
            logging.info("Download completed successfully.")
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while downloading scores: {e}")
            return

    def get_reference_allele(self, fasta_file, chrom, pos):
        """
        Fetches the reference allele from the given FASTA file at the specified chromosome and position.

        Parameters:
        fasta_file (str): Path to the reference FASTA file (can be gzipped).
        chrom (str): Chromosome (e.g., 'chr22').
        pos (int): Position (1-based).

        Returns:
        str: The reference allele at the specified location.
        """
        # Load the reference genome using pysam
        fasta = pysam.FastaFile(fasta_file)
        
        # Fetch the reference allele (0-based indexing in pysam, so subtract 1 from position)
        try:
            ref_allele = fasta.fetch(chrom, pos - 1, pos)
        except Exception as e:
            # print(f"Fetching ref allele failed for {chrom}:{pos}. Trying with chr{chrom}.")
            try:
                ref_allele = fasta.fetch(f"chr{chrom}", pos - 1, pos)
            except ValueError as e2:
                raise ValueError(f"Failed to fetch reference allele for chrom {chrom} at pos {pos}. Error: {e2}")
  
        # Close the FASTA file after fetching
        fasta.close()
        
        return ref_allele

    def _extract_snps(self, check_ref):
        logging.info(f"Extracting SNPs from scoring file for {self.pgs_id}.")
        snps = {}
        haplotype_skipped = False
        # TODO need more checking for complex data types
        for _, row in self.data.iterrows():
            if 'is_haplotype' in row and row['is_haplotype']:
                self.haplotype += 1
                haplotype_skipped = True
                continue

            try:
                # Chromosome parsing
                try:
                    try:
                        snp_chrom = str(int(row['hm_chr']))
                    except ValueError:
                        if row['hm_chr'] in ['X', 'Y']:
                            snp_chrom = row['hm_chr']
                        else:
                            raise ValueError(f"Unexpected chromosome value: {row['hm_chr']}")
                except Exception as e:
                    logging.error(f"Error parsing chromosome in row {row}: {e}")

                # Position parsing
                try:
                    snp_pos = int(row['hm_pos'])
                except Exception as e:
                    logging.error(f"Error parsing position in row {row}: {e}")

                # SNP dictionary initialization
                try:
                    if snp_chrom not in snps.keys():
                        snps[snp_chrom] = {}
                    if snp_pos not in snps[snp_chrom].keys():
                        snps[snp_chrom][snp_pos] = []
                except Exception as e:
                    logging.error(f"Error initializing SNP dictionary for chrom {snp_chrom}, pos {snp_pos}: {e}")

                # Construct the current SNP data
                try:
                    cur = {
                        'score': self.pgs_id,
                        'effect_allele': row['effect_allele'],
                        'other_allele': tuple(row['other_allele'].split(',')) if 'other_allele' in row else None,
                        'rsID': row['hm_rsID'],
                        'effect_weight': row['effect_weight']
                    }
                    snps[snp_chrom][snp_pos].extend([cur])
                except Exception as e:
                    logging.error(f"For PGS {self.pgs_id}: Error adding SNP data for chrom {snp_chrom}, pos {snp_pos}: {e}")

                # Multiallelic check
                try:
                    if len(set(cur['effect_allele'])) > 1:  
                        self.multiallelic += 1
                    if cur['other_allele'] and len(set(cur['other_allele'])) > 1:  
                        self.multiallelic += 1
                except Exception as e:
                    logging.error(f"For PGS {self.pgs_id}: Error checking multiallelic alleles for chrom {snp_chrom}, pos {snp_pos}: {e}")

                # Reference effect check (with FASTA lookup)
                if check_ref:
                    try:
                        ref_allele = self.get_reference_allele(check_ref, snp_chrom, snp_pos)

                        if isinstance(ref_allele, list):
                            ref_allele = ref_allele[0]

                        if isinstance(cur['effect_allele'], str):
                            effect_alleles = {cur['effect_allele']}
                        else:
                            effect_alleles = set(cur['effect_allele'])

                        if ref_allele in effect_alleles:
                            self.ref_effect += 1
                            # logging.info(f"Ref = effect for {snp_chrom}:{snp_pos}; ref: {ref_allele}; eff: {cur['effect_allele']}")
                    except Exception as e:
                        logging.error(f"For PGS {self.pgs_id}: Error updating reference effect for chrom {snp_chrom}, pos {snp_pos}: {e}")

                # Ambiguous strand check
                try:
                    if cur['other_allele'] is not None:
                        if (cur['effect_allele'] == 'A' and cur['other_allele'] in ['T']) or \
                        (cur['effect_allele'] == 'T' and cur['other_allele'] in ['A']) or \
                        (cur['effect_allele'] == 'G' and cur['other_allele'] in ['C']) or \
                        (cur['effect_allele'] == 'C' and cur['other_allele'] in ['G']):
                            self.ambiguous_strands += 1
                except Exception as e:
                    logging.warning(f"For PGS {self.pgs_id}: Error checking ambiguous strands for chrom {snp_chrom}, pos {snp_pos}: {e}")

            except Exception as e: 
                logging.warning(f"For PGS {self.pgs_id}, due to error {e} unable to add row {row}.")

            if haplotype_skipped:
                warning_message = f"Warning: Haplotype SNPs are skipped in PGS ID {self.pgs_id}. Using this scoring file may result in an inaccurate score."
                logging.warning(warning_message)

        return snps

    def save(self, save_path):
        with open(save_path, 'wb') as f:
            pickle.dump(self, f)
        logging.info(f"PGScore object saved to {save_path}.")

    def load(self, load_path):
        with open(load_path, 'rb') as f:
            loaded_obj = pickle.load(f)
        self.__dict__.update(loaded_obj.__dict__)
        logging.info(f"PGScore object loaded from {load_path}.")


def setup_logging(job_id):
    """Set up logging with the given job ID."""
    log_dir = f'/N/project/compgen/PGSCalc/logs/{job_id}'
    os.makedirs(log_dir, exist_ok=True)
    log_filename = f'{log_dir}/{job_id}.log'
    
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(filename=log_filename, level=logging.INFO,
                            format='%(asctime)s %(levelname)s: %(message)s')
    print(log_filename)
    return log_filename

def get_logger(name):
    """Get the logger instance for the given name."""
    return logging.getLogger(name)

def generate_pgs_id(index):
    if index < 1 or index > 5085:
        raise ValueError("Index must be between 1 and 5085.")
    return f"PGS{index:06d}"

def generate_pgs_ids_in_range(start, end):
    """Generate PGS IDs within a specified range."""
    if start < 1 or end > 5085 or start > end:
        raise ValueError("Invalid index range provided.")
    return [generate_pgs_id(i) for i in range(start, end + 1)]

def main(job_id, start_index, end_index, output_file):
    # Set up logging
    setup_logging(job_id)
    logger = get_logger(__name__)

    logger.info("Logging setup complete.")

    stats = []
    all_ids = generate_pgs_ids_in_range(start_index, end_index)
    
    for pgs_id in all_ids:
        current_score = PGScore(pgs_id=pgs_id, scoring_file_path=f'/N/project/compgen/PGSCalc/src/PGScore/scoring_files/{pgs_id}_hmPOS_GRCh38.txt.gz', check_ref='/N/project/compgen/shared/ref/GRCh38.p14.min.fa.gz', override=True)
        cur_row = {
            'pgs_id': current_score.pgs_id,
            'total_variants': sum(len(positions) for positions in current_score.snps.values()) + current_score.haplotype,
            'multiallelic': current_score.multiallelic,
            'haplotype': current_score.haplotype,
            'ambiguous_strand': current_score.ambiguous_strand,
            'ref_effect': current_score.ref_effect
        }
        stats.append(cur_row)
        logger.info(f'{pgs_id} has {current_score.ref_effect} effect alleles that are ref')
        logger.info(f"Recorded stats for {pgs_id}")

    # Check if stats has been populated
    if len(stats) == 0:
        logger.error(f"No data collected for job {job_id}, exiting.")
        return

    # Write results to the specified output file
    try:
        with open(output_file, mode='w', newline='') as file:
            fieldnames = ['pgs_id', 'total_variants', 'multiallelic', 'haplotype', 'ambiguous_strand', 'ref_effect']
            writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()  # Write the header
            writer.writerows(stats)  # Write all the rows

        logger.info(f"Results successfully written to {output_file}")

    except Exception as e:
        logger.error(f"Error writing to file {output_file}: {e}")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some inputs for PGSCalc.')
    parser.add_argument('--job_id', type=str, required=True, help='Job ID for logging.')
    parser.add_argument('--start_index', type=int, required=True, help='Start index for PGS IDs.')
    parser.add_argument('--end_index', type=int, required=True, help='End index for PGS IDs.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output file.')

    args = parser.parse_args()
    main(args.job_id, args.start_index, args.end_index, args.output_file)
    
# class PGScore:

#     def __init__(self, pgs_id, scoring_file_path):
#         self.pgs_id = pgs_id
#         self.scoring_file_path = scoring_file_path
#         self.data = self._parse_scoring()
#         self.snps = self._extract_snps()

#     def _parse_scoring(self):
#         logging.info(f"Parsing scoring file for {self.pgs_id}.")
#         try:
#             if self.scoring_file_path.endswith('.gz'):
#                 data = pd.read_csv(self.scoring_file_path, compression='gzip', sep='\t', comment='#')
#             else:
#                 data = pd.read_csv(self.scoring_file_path, sep='\t', comment='#')
#         except Exception as e:
#             logging.error(f"Error reading scoring file {self.scoring_file_path}: {e}")
#             raise
#         return data
        
#     def _extract_snps(self):
#         logging.info(f"Extracting SNPs from scoring file for {self.pgs_id}.")
#         snps = []
#         haplotype_skipped = False
#         for _, row in self.data.iterrows():
#             if 'is_haplotype' in row and row['is_haplotype']:
#                 haplotype_skipped = True
#                 continue
#             snp = {
#                 'chrom': str(row['hm_chr']),
#                 'pos': row['hm_pos'],
#                 'effect_allele': row['effect_allele'],
#                 'other_allele': tuple(row['other_allele'].split(',')),
#                 'rsID': row['hm_rsID'],
#                 'effect_weight': row['effect_weight']
#             }
#             snps.append(snp)
#         if haplotype_skipped:
#             warning_message = f"Warning: Haplotype SNPs are skipped in PGS ID {self.pgs_id}. Using this scoring file may result in an inaccurate score."
#             # print(warning_message)
#             logging.warning(warning_message)
#         return snps