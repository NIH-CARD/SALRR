#!/usr/bin/env python3
import os
import pandas as pd
import glob
import re
import argparse
import polars as pl
from gtfparse import read_gtf
import yaml
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

"""
SIRV Spike-In QC Analysis Tool

Analyzes SIRV RNA spike-in controls from ONT long-read RNA-seq data.
Compares expected vs measured concentrations and generates QC plots.

Usage:
    python sirvanalyzer.py ./config/biowulf.yaml

Requirements:
    - SIRV reference file: SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx
    - Completed StringTie assembly GTF files
    - Config YAML with base_dir, stringtie_dir, sample_file paths

Example running of the script:
    python sirvanalyzer.py ./config/biowulf.yaml
"""



###############################################################################################################################################
#
# SIRV Analysis Pipeline
###############################################################################################################################################
class DotDict(dict):
    """DotDict class allows accessing dictionary keys as attributes."""
    def __getattr__(self, attr):
        if attr in self: return self[attr]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{attr}'")
    def __setattr__(self, key, value): self[key] = value

class SIRVAnalyzer:
    def __init__(self, config_file):
        """Initialize with config and set up state"""
        self.config_file = config_file
        self.config = self._load_config()
        
        # State attributes
        self.expected_df = None
        self.count_files = []
        self.merged_df = None
        
    def _load_config(self):
        with open(self.config_file, 'r') as file:
            return DotDict(yaml.safe_load(file))

    def load_expected_values(self, sheet_name='ERCC', 
                           file_name='SIRV_Set4_Norm_sequence-design-overview_20210507a.xlsx'):
        """Loads and processes reference SIRV values"""
        data_dir = self.config.get('data_dir', '.') 
        file_path = os.path.join(data_dir, file_name)
        
        if not os.path.exists(file_path):
             raise FileNotFoundError(f"Reference file not found: {file_path}")

        print(f"Loading expected values from {sheet_name}...")
        
        # Load and flatten columns (your existing logic)
        df = pd.read_excel(file_path, sheet_name=sheet_name, header=[3,4])
        df.columns = [' '.join([str(c) for c in col if pd.notna(c)]).strip() for col in df.columns]

        # Process based on sheet type
        if sheet_name == 'ERCC':
            desired_columns = {
                'ERCC ID Unnamed: 1_level_1': 'ERCC ID',
                'GenBank* Unnamed: 2_level_1': 'GenBank*',
                'stock conc (amoles/µl)': 'conc (amoles/µl)',
                'full ERCC transcript, including poly(A) tail length (nt)': 'length (nt)',
                'full ERCC transcript, including poly(A) tail GC': 'GC'
            }
            df_subset = df[list(desired_columns.keys())].rename(columns=desired_columns)
            df_subset = df_subset[df_subset['ERCC ID'].astype(str).str.startswith("ERCC-")]
            
            # Create merge key
            df_subset['merge_key'] = (df_subset['GenBank*'].astype(str) + "_" + 
                                    df_subset['ERCC ID'].astype(str))
            
        elif sheet_name == 'longSIRV':
            desired_columns = {
                'long SIRV ID Unnamed: 1_level_1' : 'longSIRV ID',
                'conc (amoles/µl) Unnamed: 4_level_1': 'conc (amoles/µl)',
                'full transcript, including poly(A) tail length (nt)' : 'length (nt)',
                'full transcript, including poly(A) tail GC' : 'GC'
            }
            df_subset = df[list(desired_columns.keys())].rename(columns=desired_columns)
            df_subset = df_subset[df_subset['longSIRV ID'].astype(str).str.startswith('SIRV')]
            # Add merge_key logic for longSIRV if needed
        
        self.expected_df = df_subset
        return self.expected_df

    def parse_gtf_files(self, input_sample=None):
        """Parses GTF files defined in sample sheet"""
        print("Parsing GTF files...")
        if input_sample is None:
            input_sample = self.config.sample_file
            
        samples = pd.read_csv(input_sample, sep="\t", header=None, engine='python')
        
        self.count_files = {}
        for _, row in samples.iterrows():
            sample_id, flowcell_id = row[0], row[1]
            gtf_path = f"{self.config.base_dir}{self.config.stringtie_dir}/{sample_id}/sirv/"
            gtf_file = os.path.join(gtf_path, f"{sample_id}_{flowcell_id}_sirv.stringtie.gtf")
            output_file = os.path.join(gtf_path, f"{sample_id}_{flowcell_id}_sirv_counts.tsv")
            
            # Checking if GTF file exists before trying to read it
            if not os.path.exists(gtf_file):
                print(f"Warning: SIRV GTF file not found: {gtf_file}")
                continue  # Skip this file and move to the next sample

            # Use polars logic to read and filter GTF
            gtf = read_gtf(gtf_file)
            features = gtf.filter(pl.col("feature") == "transcript")[
                ["reference_id", "ref_gene_id", "FPKM"]
            ]
            features.write_csv(output_file, separator="\t")
            self.count_files[output_file] = sample_id
            
        return self.count_files

    def merge_data(self):
        """Combines expected values with sample counts"""
        if self.expected_df is None or not self.count_files:
            raise ValueError("Run load_expected_values() and parse_gtf_files() first")

        print("Merging sample data...")
        merged_list = []
        
        # 1. Prepare Expected Reference
        df_exp = self.expected_df[['merge_key', 'conc (amoles/µl)', 'length (nt)', 'GC']].copy()
        df_exp = df_exp.rename(columns={'conc (amoles/µl)': 'conc'})
        df_exp['source'] = 'Expected'
        merged_list.append(df_exp)
        
        # 2. Process Sample Files
        for f, sample_name in self.count_files.items():
            df_sample = pd.read_csv(f, sep='\t').dropna(
                subset=['reference_id', 'ref_gene_id']
            )
            
            df_sample['merge_key'] = (
                df_sample['reference_id'].astype(str) + "_" + 
                df_sample['ref_gene_id'].astype(str)
            )
            
            # Filter and Format
            subset = df_sample[df_sample['merge_key'].isin(self.expected_df['merge_key'])].copy()
            subset = subset[['merge_key', 'FPKM']].rename(columns={'FPKM': 'conc'})
            subset['source'] = sample_name
            merged_list.append(subset)
            
        self.merged_df = pd.concat(merged_list, ignore_index=True)
        
        # Save output
        out_path = os.path.join(self.config.base_dir, "merged_sirv_results.csv")
        self.merged_df.to_csv(out_path, index=False)
        return self.merged_df

    def plot_metrics(self, output_path="SIRV_analysis"):
        """Generates QC plots using merged data"""
        if self.merged_df is None:
            raise ValueError("Run merge_data() first")
            
        print("Generating plots...")
        
        # Create output directory for plots
        output_dir = f"{self.config.base_dir}{self.config.stringtie_dir}/plots/"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_path)

        merged_df = self.merged_df  # Local reference for convenience
        
        # Separate Expected values from Sample values
        df_expected = merged_df[merged_df['source'] == 'Expected'].copy()
        df_samples = merged_df[merged_df['source'] != 'Expected'].copy()
        
        # Create mapping dictionaries from Expected rows (for metadata lookup)
        expected_lookup = df_expected.set_index('merge_key')
        len_map = expected_lookup['length (nt)'].to_dict()
        gc_map = expected_lookup['GC'].to_dict()
        conc_map = expected_lookup['conc'].to_dict()
        
        # Add metadata columns to sample dataframe
        df_samples['expected_conc'] = df_samples['merge_key'].map(conc_map)
        df_samples['length'] = df_samples['merge_key'].map(len_map)
        df_samples['gc'] = df_samples['merge_key'].map(gc_map)
        
        # --- PLOT 1: Barplot - Number of Transcripts Detected per Sample ---
        sample_df = merged_df[merged_df['source'] != 'Expected']
        detected_counts = sample_df.groupby('source').apply(lambda x: (x['conc'] > 0).sum())
        average_detected = detected_counts.mean()
        
        plt.figure(figsize=(10, 6))
        ax = detected_counts.plot(kind='bar', color='steelblue', alpha=0.8)
        plt.axhline(average_detected, color='red', linestyle='--', linewidth=2, 
                   label=f"Average: {average_detected:.1f}")
        plt.xlabel("Sample", fontsize=14)
        plt.ylabel("Transcripts Detected", fontsize=14)
        plt.title("Number of ERCC Transcripts Detected per Sample", fontsize=16)
        plt.xticks(rotation=45, ha='right', fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(f"{output_path}_detection_counts.png", dpi=150)
        plt.close()
        print(f"Saved: {output_path}_detection_counts.png")

        # --- PLOT 2: Expected vs Measured - ALL SAMPLES ON ONE PLOT ---
        df_plot_valid = df_samples[(df_samples['conc'] > 0) & (df_samples['expected_conc'] > 0)].dropna(subset=['conc', 'expected_conc'])
        
        if len(df_plot_valid) > 0:
            plt.figure(figsize=(8, 8))
            sns.scatterplot(data=df_plot_valid, x='expected_conc', y='conc', 
                        hue='source', alpha=0.7, s=50)
            
            plt.xscale('log')
            plt.yscale('log')
            plt.title("Expected vs Measured Concentration (All Samples)")
            plt.xlabel("Expected Concentration (amoles/µl)")
            plt.ylabel("Measured FPKM")
            
            # Calculate overall Pearson correlation
            log_expected = np.log10(df_plot_valid['expected_conc'])
            log_measured = np.log10(df_plot_valid['conc'])
            r, p_val = pearsonr(log_expected, log_measured)
            
            plt.text(0.05, 0.95, f"Overall Pearson r = {r:.2f}\np = {p_val:.2e}",
                    transform=plt.gca().transAxes, verticalalignment='top',
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
            
            plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(f"{output_path}_all_samples_correlation.png", dpi=150, bbox_inches='tight')
            plt.close()
            print(f"Saved: {output_path}_all_samples_correlation.png")
        else:
            print("No valid data for correlation plot")
        
        # PLOT 3: Detection Rate Curve
        all_transcripts = self.expected_df['merge_key'].unique()
        all_samples = df_samples['source'].unique()
        n_samples = len(all_samples)

        # Count how many samples detected each transcript (conc > 0)
        detection_counts = df_samples[df_samples['conc'] > 0].groupby('merge_key').size()

        # Create detection rate for ALL expected transcripts (including those never detected)
        detection_rate_dict = {}
        for transcript in all_transcripts:
            if transcript in detection_counts.index:
                detection_rate_dict[transcript] = detection_counts[transcript] / n_samples
            else:
                detection_rate_dict[transcript] = 0.0  # Not detected in any sample

        detection_rate_series = pd.Series(detection_rate_dict).sort_values(ascending=False).reset_index(drop=True)

        plt.figure(figsize=(10, 4))
        plt.scatter(detection_rate_series.index, detection_rate_series.values, alpha=0.6, color='deepskyblue', s=50)
        plt.xlabel("Transcripts (sorted by detection rate)")
        plt.ylabel("Detection Rate")
        plt.title("Detection Rate per Transcript")
        plt.ylim(0, 1.05)
        plt.axhline(y=0.5, color='red', linestyle='--', linewidth=2, label='50% Detection Rate')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{output_path}_detection_rate.png", dpi=150)
        plt.close()
        print(f"Saved: {output_path}_detection_rate.png")

        # --- PLOT 4: Distribution Boxplots by Detection Group ---
        df_plot = pd.DataFrame({
            'merge_key': all_transcripts,
            'detection_rate': [detection_rate_dict[t] for t in all_transcripts]
        })

        # Map metadata back
        df_plot['length'] = df_plot['merge_key'].map(len_map)
        df_plot['gc'] = df_plot['merge_key'].map(gc_map)
        df_plot['expected_conc'] = df_plot['merge_key'].map(conc_map)

        # Create detection groups
        df_plot['Detection_Group'] = df_plot['detection_rate'].apply(
            lambda x: '>=50%' if x >= 0.50 else '<50%'
        )

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        # Concentration
        sns.boxplot(data=df_plot, x='Detection_Group', y='expected_conc', ax=axes[0], palette='Set2')
        axes[0].set_title("Concentration Distribution by Detection Group", fontsize=14)
        axes[0].set_xlabel("Detection Group", fontsize=12)
        axes[0].set_ylabel("Concentration (amoles/µl)", fontsize=12)
        axes[0].set_yscale('log')
        
        # Length
        sns.boxplot(data=df_plot, x='Detection_Group', y='length', ax=axes[1], palette='Set2')
        axes[1].set_title("Transcript Length Distribution by Detection Group", fontsize=14)
        axes[1].set_xlabel("Detection Group", fontsize=12)
        axes[1].set_ylabel("Length (nt)", fontsize=12)
        
        # GC Content
        sns.boxplot(data=df_plot, x='Detection_Group', y='gc', ax=axes[2], palette='Set2')
        axes[2].set_title("GC Content Distribution by Detection Group", fontsize=14)
        axes[2].set_xlabel("Detection Group", fontsize=12)
        axes[2].set_ylabel("GC Content", fontsize=12)
        
        plt.tight_layout()
        plt.savefig(f"{output_path}_distributions.png", dpi=150)
        plt.close()
        print(f"Saved: {output_path}_distributions.png")
        
        print(f"\n All plots generated successfully with prefix: {output_path}")

    def run(self):
        """running the full pipeline"""
        if not self.config.get('use_sirv', True):
            print("SIRV analysis disabled, skipping...")
            return

        print("Starting SIRV Analysis Pipeline...")
        self.load_expected_values()
        self.parse_gtf_files()
        self.merge_data()
        self.plot_metrics()
        print("SIRV Analysis Pipeline completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file')
    args = parser.parse_args()
    
    analyzer = SIRVAnalyzer(args.config_file)
    analyzer.run()
###########################################################################################################################################################
# End of SIRV Analysis Pipeline
###########################################################################################################################################################