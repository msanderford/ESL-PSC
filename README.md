# Evolutionary Sparse Learning for Convergent Traits (ESL-PSC) #

## Table of Contents ##

1. [Description](#description)
2. [Usage](#usage)
3. [Using a Configuration File with ESL-PSC Scripts](#using-a-configuration-file-with-esl-psc-scripts)
4. [Requirements](#requirements)
5. [Input Data](#input-data)
6. [Output Data](#output-data)
7. [Additional Options and Parameters](#additional-options-and-parameters)
8. [Included Data](#included-data)

## Description ##
ESL-PSC is a set of scripts for performing integrations of many species combinations using the Evolutionary Sparse Learning method. The tool is designed for analyzing genomic data to identify convergence in phenotypic traits among species. It employs a combination of Sparse Group Lasso and Contrast Tree methods to perform the analysis.

The main script, esl_multimatrix.py, takes in various input parameters and options to control the analysis process. It preprocesses input data, performs gap-cancellation, creates response matrices, and generates predictions. The output of the script includes gene ranks, species predictions, and optional plots to visualize the results.

## Usage ##
To use ESL-PSC, you will need to run the esl_multimatrix.py script with the necessary arguments and options. You can provide the input parameters and options through the command line or by creating a configuration file called esl_psc_config.txt. When using a configuration file, provide one argument per line.

Here is an example of how to run the script:

`python esl_multimatrix.py --esl_inputs_outputs_dir /path/to/inputs_outputs/dir --output_file_base_name output_file_name --species_groups_file /path/to/species_groups_file  --alignments_dir /path/to/alignments/dir --use_logspace`

To see all of the options available for any of the scripts in this directory, you can use `python [script_name].py --help`

### Using a Configuration File with ESL-PSC Scripts ###

ESL-PSC scripts can utilize a configuration file to easily manage arguments that remain constant across multiple runs. The scripts will check for the presence of an esl_psc_config.txt file in the current working directory. If it exists, the function reads the arguments from the file and combines them with any additional command-line arguments provided when running the script. This allows you to keep common arguments, such as `--esl_main_dir`, in the configuration file, while providing run-specific arguments via the command line.

To use this feature:

1. An existing esl_psc_config.txt file is included in the directory with the ESL-PSC scripts. 

2. Add any desired arguments to the esl_psc_config.txt file, using the same format as you would when providing them on the command line. Each argument should be on a new line, followed by its corresponding value. For example:

`--esl_main_dir /path/to/main/directory`

3. Save the `esl_psc_config.txt` file.

4. When running an ESL-PSC script, the scripts will automatically check for the presence of the esl_psc_config.txt file and incorporate its contents. Command-line arguments will override the values in the config file if both are provided for the same argument.

5. If the esl_psc_config.txt file is not found in the current working directory, the function will only parse command-line arguments.

It is recommended to always include the esl_main_dir argument in the configuration file, as it typically remains unchanged between runs. A configuration file with this argument is already included in the repository.

## Requirements ##

To run ESL-PSC, you will need the following software, libraries, and dependencies:

Python
ESL-PSC requires Python 3. It has not been tested with Python 2, and compatibility is not guaranteed.

Libraries
The following Python libraries are required to run ESL-PSC:

- BioPython
- NumPy
- pandas
- matplotlib
- seaborn

You can install these libraries using pip:

`pip install biopython numpy pandas matplotlib seaborn`

## Input Data ##

#### The main input files required for ESL-PSC are: ####

1. A directory of alignment files. These should be in 2-line fasta format and whose file names must have the file extension `.fas`. It is assumed that each seperate alignment file will be a different genomic component, such as a gene, a protein, an exon, a domain, etc. and each component will be treated as a "group" of sites in the analysis (see Methods in Allard et al., 2023). Use the argumemnt `--alignments_dir` and give the full absolute path to the directory.

2. A species groups file.  This is a text file that contains a comma delimited list of species on each line. In the simplest case, one species identifier can be placed on each line. The first line must contain one or more species that possess the convergent trait under analysis, and the next line must contain one or more species that can serve as trait-negative controls for the species in the first line, such that the first two lines, and each subsequent pair of lines will define a contrast pair of species to use in the analysis (see Allard et al., 2023 for details on chosing contrast pairs for ESL-PSC analysis). When more than one species is given in a line, each of those species will be used in a seperate analysis, along with all combinations of other alternative speices.  Thus, the total number of species combinations can be calculated by the product of the number of species given on each line. Use the argument `--species_groups_file` and give the full absolute path to the file.

#### Optional input files: ####

1. A species phenotype file. This is a text file which has each in the full  species name followed by a comma and then a 1 or -1 for the true phenotype class to which that species belongs. A 1 typically refers to the convergent phenotype. If this file is not provided, the ture phenotype will not be listed for each species prediction in the species_predictions output file.

2. A directory of alignments to use for preditions. By default, any species in the input alignments that are not used in building any given model will be assigned a prediction score (SPS) for that model, which will be included in the predictions output file. As an alternative, you can use a seperate directory of alignments for the predictions, however these still need to be fully aligned to any input species alignments or the predictions will be meaningless. Use the argument `--prediction_alignments_dir` and give the full absolute path to the directory.

3. Canceled alignments directory. Full path to the new alignments directory. Gap-canceled alignments for each species combo will be placed here. This may also be an existing folder of gap-canceled alignments for multimatrix ESL. Use the argument `--canceled_alignments_dir` and give the full absolute path to the directory.

4. Limited genes list. If you want to use a subset of the alignment files for model building without having to remove files from your alignments directory, you can submit a limited genes list file, which is a text file containing one alignment file name on each line. Note that these names must exactly match the ones in the alignments directory, and must end in `.fas` like they do. Use the argument `--limited_genes_list` and give the full absolute path to the file.

## Output Data ##

ESL-PSC generates two main types of output files: a Predictions File and a Gene Ranks File.

#### Predictions File ####
The predictions file contains every prediction made by every model generated using every species combination in the analysis. Each line in the file lists the following information:

1. Species combination
2. Lambda1 (first sparsity hyperparameter)
3. Lambda2 (second sparsity hyperparameter)
4. Penalty term
5. Number of genes
6. Input Root Mean Squared Error (RMSE; this is referred to as the Model Fit Score by Allard et al. (2023))
7. Species being predicted
8. Species-Phenotype Score (SPS)
9. True phenotype for the species (taken from the species_pheno_file if provided)

#### Gene Ranks File ####
The gene ranks file lists the genes (or proteins or other genomic components) used in the analysis, along with information about their rankings based on their model contributions. Each line in the file includes the following information:

1. Gene name (taken from the alignment file)
2. Number of species combinations in which the gene is ranked (i.e. number of combinations for which it recieved a non-zero GSS as part of any model)
3. Number of species combinations in which the gene is ranked among the top contributors (the percentage of genes to consider "top genes" by GSS in any model can be set using the `--top_rank_frac` argument.)
4. Highest ever Gene Selection Score (GSS)
5. Best ever rank (the best ever rank, 1 being the best possible, recieved in any model)

## Additional Options and Parameters ##

The following additional options and parameters can be specified when running ESL-PSC to fine-tune the analysis and control various aspects of the process. These options can be added as command line arguments or specified in the config file.

##### Hyperparameters:
* `--initial_lambda1`: Initial lambda 1 value (position sparsity parameter). Default = .01.
* `--final_lambda1`: Final lambda 1 value (position sparsity parameter). Default = .99.
* `--initial_lambda2`: Initial lambda 2 value (group sparsity parameter). Default = .01.
* `--final_lambda2`: Final lambda 2 value (group sparsity parameter). Default = .99.
* `--lambda_step`: The increment to increase the lambda values with each step. It is recommended to use a logspace (see options below) but in a linear gridsearch of sparsity hyperparameters, this controls the step between values.
* `--group_penalty_type`: Group penalty calculation type ("sqrt", "median", "linear", or "default"). Median will be used by default (see Methods in Allard et al. 2023)
* `--initial_gp_value`: Group penalty constant term initial value. If a linear group lenalty type is selected, the group penalties for each gene will be equal to the number of variable sites in the gene's alignment plus a constant term that is the same across all genes. By default, this will be 1 for all genes, but it is also possible to use a range of different constant terms and repeat all model ensembles for each group penalty term. In order to do this, the initial, final and step can be set using this and the following two arguments.
* `--final_gp_value`: Group penalty constant term final value. See initial_gp_value above for explanation.
* `--gp_step`: Group penalty constant term increment. The default is 6. See initial_gp_value above for explanation.
* `--num_log_points`: The number of values per sparsity hyperparameter (lambda1 and lambda2) in a logspace of values to test.
* `--pheno_names`: The names of the two phenotypes separated by a space, with the convergent phenotype coming first. by default "1" and "-1" will be used
* `--min_genes`: Minimum number of genes a model must have in order for that model to be included in the prediction scores plots. Default = 0.

##### Options:
* `--use_logspace`: *Recommended* Use a log space of points for lambda values instead of initial and final lambda values with a lambda step.
* `--use_existing_preprocess`: Use existing preprocess folder and skip running the preprocess step.
* `--use_default_gp`: Don't replace group penalties (automatically set to True if the group_penalty_type is "default").
* `--keep_raw_output`: Don't delete the raw model output files for each run.
* `--show_selected_sites`: Print a dictionary of all selected sites with their highest model score for every gene in the gene_ranks output file.
* `--no_genes_output`: Don't output a gene ranks file. If only predictions output is desired, including the option will speed up the analysis.
* `--no_pred_output`: Don't output a species predictions file. If only gene ranks output is desired, including the option will significantly speed up the analysis.
* `--make_sps_plot`: Make a violin plot showing SPS density for each true phenotype (SPS of > 1 or < -1 as 1 and -1 by default).
* `--make_sps_kde_plot`: Make a KDE plot showing SPS density for each true phenotype.

##### Deletion Canceler Options:
* `--nix_full_deletions`: Don't create files for fully canceled genes.
* `--cancel_only_partner`: Only cancel partner of any gap species at the site.
* `--min_pairs MIN_PAIRS`: The minimum number of pairs that must not have gaps or the whole site will be canceled.
* `--limited_genes_list`: Use only genes in this list. One file per line.

##### Multimatrix-specific Optional Arguments:
* `--top_rank_frac`: Fraction of genes to count as "top genes."  The default is .01 (1%)
* `--response_dir`: Folder with response matrices. Any txt file in this folder is assumed to be a response matrix file.
* `--use_uncanceled_alignments`: Use the alignments_dir alignments for all matrices without doing gap canceling (not recommended).
* `--use_existing_alignments`: Use existing files in canceled_alignments_dir.
* `--delete_preprocess`: Clear preprocess folders after each matrix run.
* `--make_null_models`: Make null response-flipped ESL-CT models. Must have an even number of pairs. All balanced flippings of the response values will be generated for each combo and all will be run and aggregated to maximally decouple true convergences (see Methods in Allard et al. 2023). 
* `--make_pair_randomized_null_models`: Make null pair randomized ESL-CT models. A copy of input deletion-canceled alignment will, for each variable site, be randomized such that the residues of each contrast pair will be either flipped or not and the ESL integration will be repeated for each one. The results are then aggregated for all (see Methods in Allard et al. 2023).
* `--num_randomized_alignments`: Number of pair-randomized alignments to make. Default is 10.

## Included Data ##

#### We have included two sample species_groups files for use in ESL-PSC alignments ####
1. photo_single_LC_matrix_species_groups.txt (the grass species with the closest contrast partners with the longest sequences (fewest gaps used for photosynthesis analyses by Allard et al. (2023))
2. orthomam_echo_species_groups.txt (this can be used to reproduce the echolocation analyses using all 16 species combinations (Allard et al. 2023) 

#### We have included the protein sequence alignments used for ESL-PSC analyses by Allard et al. (2023). If you use these data, please cite these sources: ####


##### Grass chloroplast alignments which were used by Allard et al. (2023) were derived from:

Casola C, Li J. 2022. Beyond RuBisCO: convergent molecular evolution of multiple chloroplast genes in C4 plants. PeerJ 10:e12791 https://doi.org/10.7717/peerj.12791
More information regarding these alignments can be found in the supplemental information kindly provided online by these authors.

##### Mammalian protein sequence alignments for echolocators and their control species were derived from the OrthoMaM database:
https://orthomam.mbb.cnrs.fr/#

OrthoMaM v10: Scaling-Up Orthologous Coding Sequence and Exon Alignments with More than One Hundred Mammalian Genomes Celine Scornavacca, Khalid Belkhir, Jimmy Lopez, Rémy Dernat, Frédéric Delsuc, Emmanuel J P Douzery, Vincent Ranwez Molecular Biology and Evolution, Volume 36, Issue 4, April 2019, Pages 861–862

## Citation ##
If you use this software in your research, please cite our paper:

Allard, J. B., Sharma, S., Patel, R., Sanderford, M., Tamura, K., Vucetic, S., Gerhard, G. S., & Kumar, S. (2023). Evolutionary sparse learning reveals the shared genetic basis of convergent traits. Nature Ecology and Evolution. Submitted.


