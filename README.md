# AL2 Antibody Stitching Analysis

## How to run the alignment script

First install the Python dependencies (Python 3.9+):

```bash
pip install -r requirements.txt
```

This will output the aligned files and the figures to .\data\align_output\aligned_files and .\data\align_output\figures respectively.

```bash
python .\data\align_results.py .\data\casanovo_results .\data\align_output --percentile 0.90 --min_ppm_error -20 --max_ppm_error 20 --keep_n_peptides 1 
```

## How to run our stitch analysis

1) Download Stitch from GitHub releases (I used version 1.5.0): https://github.com/snijderlab/stitch
2) Unzip the download
3) Copy the files found in this repo under the `stitch_files` directory to their respective directories from the extracted Stitch download.
4) Copy the data folder from this repo to the extracted Stitch directory
5) cd into the extracted Stitch directory
6) run `.\stitch.exe run .\batchfiles\AL2_monoclonal.txt`

Note: There are 4 provided templates in this repo. The only difference between the AL2 versions and the 
versions included in stitch is the addition of 1 additional joining region (>IGLJ2c IGLJ3c), and several additional
contaminants found in the AL2 samples.