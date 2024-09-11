import os
import argparse
import pandas as pd
import peptacular as pt
from matplotlib import pyplot as plt
from typing import Tuple, List


# Function to parse Casanovo .mztab files and return a dataframe along with the introductory metadata.
def casanovo_to_df(file: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    header, data, intro = None, [], []
    for line in file:
        if line == '\n':  # Skip empty lines
            continue
        elif line.startswith('PSH'):  # PSH indicates the header line
            header = line.strip().split('\t')
        elif line.startswith('PSM'):  # PSM lines contain peptide data
            data.append(line.strip().split('\t'))
        else:
            intro.append(line)  # Store the introduction section
    df = pd.DataFrame(data, columns=header)  # Create dataframe with headers
    return df, intro


# Function to fix data types in the dataframe
def fix_df(df: pd.DataFrame) -> None:
    df['exp_mass_to_charge'] = df['exp_mass_to_charge'].astype(float)
    df['calc_mass_to_charge'] = df['calc_mass_to_charge'].astype(float)
    df['search_engine_score[1]'] = df['search_engine_score[1]'].astype(float)
    df['charge'] = df['charge'].astype(float).astype(int)


# Function to modify the dataframe by converting sequences and calculating error metrics
def convert_df(df: pd.DataFrame) -> None:
    df['sequence'] = df['sequence'].apply(pt.convert_casanovo_sequence)
    df['unmodified_sequence'] = df['sequence'].apply(pt.strip_mods)  # Strip modifications from sequences
    df['ppm_error'] = df.apply(lambda x: pt.ppm_error(x['calc_mass_to_charge'], x['exp_mass_to_charge']), axis=1)
    df['dalton_error'] = df.apply(lambda x: x['calc_mass_to_charge'] - x['exp_mass_to_charge'], axis=1)


# Function to align mass-to-charge ratios and generate histograms of errors
def align_file(infile: str, aligned_dir: str, fig_dir: str, percentile: float = 0.90, min_ppm_error: float = -20,
               max_ppm_error: float = 20, keep_n_peptides: int = -1) -> None:
    print(f'Aligning {infile}')
    base_filename = os.path.basename(infile).replace('.mztab', '')

    # Set output file and figure paths
    output_base = os.path.join(aligned_dir,
                               base_filename + f'{"_" + str(keep_n_peptides) + "n_peptides" if keep_n_peptides > 0 else ""}')
    output_error_fig = os.path.join(fig_dir, base_filename + '_ppm_error_histogram.png')
    output_C_term_fig = os.path.join(fig_dir, base_filename + '_C_term_histogram.png')
    output_N_term_fig = os.path.join(fig_dir, base_filename + '_N_term_histogram.png')
    outfile = os.path.join(aligned_dir, base_filename + '_aligned.mztab')

    # Read the input file and process the data into a dataframe
    with open(infile) as f:
        df, intro = casanovo_to_df(f)
        fix_df(df)

    # If specified, keep only the best scoring peptide per sequence and charge pair
    if keep_n_peptides > 0:
        df = df.sort_values('search_engine_score[1]', ascending=False).groupby(['sequence', 'charge']).head(1)
        print(f'Keeping {keep_n_peptides} peptides')

    # Filter top peptides based on percentile score and calculate errors
    top_df = df.nlargest(int(len(df) * (1 - percentile)), 'search_engine_score[1]')
    convert_df(top_df)

    print('Percentile Score:', percentile, 'Top Score:', top_df['search_engine_score[1]'].max())

    # Filter by ppm error range
    top_df = top_df[(top_df['ppm_error'] > min_ppm_error) & (top_df['ppm_error'] < max_ppm_error)]
    before_median_ppm_error = top_df['ppm_error'].median()
    print(f'Median PPM Error (Before): {before_median_ppm_error}')

    # Align mass-to-charge values by correcting based on median ppm error
    top_df['align_exp_mass_to_charge'] = top_df['exp_mass_to_charge'] - before_median_ppm_error * top_df[
        'exp_mass_to_charge'] / 1e6
    top_df['align_dalton_error'] = top_df.apply(lambda x: x['calc_mass_to_charge'] - x['align_exp_mass_to_charge'],
                                                axis=1)
    top_df['align_ppm_error'] = top_df.apply(
        lambda x: pt.ppm_error(x['calc_mass_to_charge'], x['align_exp_mass_to_charge']), axis=1)

    after_median_ppm_error = top_df['align_ppm_error'].median()
    print(f'Median PPM Error (After): {after_median_ppm_error}')

    # Calculate the 99% confidence interval for the aligned ppm errors
    lower_bound = top_df['align_ppm_error'].quantile(0.005)
    upper_bound = top_df['align_ppm_error'].quantile(0.995)
    print(f'99% Confidence Interval: {lower_bound} to {upper_bound}')

    # Plot and save the PPM error histogram
    plt.hist(top_df['ppm_error'], bins=100, density=True, alpha=0.5, label='Before', color='b')
    plt.hist(top_df['align_ppm_error'], bins=100, density=True, alpha=0.5, label='After', color='r')
    plt.title('Histogram of PPM Error')
    plt.xlabel('PPM Error')
    plt.ylabel('Density')
    plt.axvline(before_median_ppm_error, color='b', linestyle='dashed', linewidth=1, label='Before Median')
    plt.axvline(after_median_ppm_error, color='r', linestyle='dashed', linewidth=1, label='After Median')
    plt.axvline(lower_bound, color='g', linestyle='dashed', linewidth=1, label='99% CI Lower Bound')
    plt.axvline(upper_bound, color='g', linestyle='dashed', linewidth=1, label='99% CI Upper Bound')
    plt.legend()
    plt.savefig(output_error_fig)
    plt.close()

    # Plot and save C-terminal amino acid histogram
    top_df['C_term'] = top_df['unmodified_sequence'].apply(lambda x: x[-1])
    top_df['C_term'].value_counts().plot(kind='bar')
    plt.title('Histogram of C Terminal Amino Acids')
    plt.xlabel('Amino Acid')
    plt.ylabel('Count')
    plt.savefig(output_C_term_fig)
    plt.close()

    # Plot and save N-terminal amino acid histogram
    top_df['N_term'] = top_df['unmodified_sequence'].apply(lambda x: x[0])
    top_df['N_term'].value_counts().plot(kind='bar')
    plt.title('Histogram of N Terminal Amino Acids')
    plt.xlabel('Amino Acid')
    plt.ylabel('Count')
    plt.savefig(output_N_term_fig)
    plt.close()

    # Align the entire dataframe using the before-median ppm error
    df['exp_mass_to_charge'] = df['exp_mass_to_charge'] - before_median_ppm_error * df['exp_mass_to_charge'] / 1e6

    # Write the aligned data to a new file
    with open(outfile, 'w') as f:
        f.writelines(intro)  # Write the introductory section
    df.to_csv(outfile, sep='\t', index=False, mode='a')  # Append the aligned data
    print('Wrote:', outfile)
    print()


# Function to parse command-line arguments
def get_args() -> argparse.Namespace:
    # Example command: data/casanovo_results data/align_output --percentile 0.90 --min_ppm_error -20 --max_ppm_error 20 --keep_n_peptides 1
    parser = argparse.ArgumentParser(description="Align mass-to-charge ratios in mztab files")
    parser.add_argument('input_folder', help='Folder containing mztab files to process')
    parser.add_argument('output_folder', help='Folder where aligned files and figures will be saved')
    parser.add_argument('--percentile', type=float, default=0.90, help='Percentile to use for filtering scores')
    parser.add_argument('--min_ppm_error', type=float, default=-20, help='Minimum ppm error for filtering')
    parser.add_argument('--max_ppm_error', type=float, default=20, help='Maximum ppm error for filtering')
    parser.add_argument('--keep_n_peptides', type=int, default=-1, help='Number of peptides to keep per charge pair')
    return parser.parse_args()


# Function to process all files in the input folder and align them
def run(input_folder: str, output_folder: str, percentile: float, min_ppm_error: float, max_ppm_error: float,
        keep_n_peptides: int) -> None:
    aligned_folder = os.path.join(output_folder, 'aligned_files')
    figures_folder = os.path.join(output_folder, 'figures')

    # Create output folders if they don't exist
    os.makedirs(aligned_folder, exist_ok=True)
    os.makedirs(figures_folder, exist_ok=True)

    # List all files in the input folder
    mztab_files = [os.path.join(input_folder, file) for file in os.listdir(input_folder) if file.endswith(".mztab")]

    # Align each file and generate outputs
    for file in mztab_files:
        align_file(file, aligned_folder, figures_folder, percentile, min_ppm_error, max_ppm_error, keep_n_peptides)


if __name__ == '__main__':
    # Parse arguments and run the alignment process
    args = get_args()
    run(args.input_folder, args.output_folder, args.percentile, args.min_ppm_error, args.max_ppm_error,
        args.keep_n_peptides)
