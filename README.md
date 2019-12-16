### Tool for protein-domain-composition-based-gene-family scoring

This container contains application for scoring gene-families using their Pfam protein domain composition. Cosine function is used to compare the protein domains of the sequences within gene-families.

 1. #### Downloading the container
  ```
  docker pull akshayayadav/genefamily-domain-composition-cosine-scoring
  ```
 2. #### Preparing the data
  Create a data directory *<my_directory>* with a user-defined name. This directory **MUST** contain a directory called *family_fasta* containing fasta files for the desired gene families and directory named *pfam_database* containing the Pfam database HMM file. The HMM file must be *hmmpressed* using the hmmpress tool from the HMMER package.

 3. #### Running the analysis
  Execute the analysis on the prepared data with *n* cores using the following command
  ```
  docker run -v <absolute_path_to_the_data_directory>:/data akshayayadav/genefamily-domain-composition-cosine-scoring run_analysis.sh -c <n>
  ```
