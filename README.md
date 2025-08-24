# IEDBTool2
Tool to match epitope predictions using IEDB Data and user provided data.

This tool has been developed using a package called EpitopeFinder (Gong and Bloom, PLoS Genetics, 2014)

Their publication can be found [here](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004328), and their Github Repository can be found [here](https://github.com/jbloom/epitopefinder).

This pipeline aimed to simplify and reduce the tedious task of manual epitope selection for IFN-γ ELISPOT assay testing by performing a two-step process and is also used in [Temporal dynamics of viral fitness and the adaptive immune response in HCV infection](https://doi.org/10.7554/eLife.102232.2.sa3) by Walker & Leung et al [^1].

1) *iedb_predictionParser.py* - Predicted epitopes generated from IEDB [prediction algorithm](http://tools.iedb.org/mhci/) downloaded in plain text format is to be parsed into this script to prepare data for the next script. The options available include boundary conditions such as *-lower N* and *-upper N* to indicate the percentile ranking (indicated by N. For instance: -lower 76.0 -upper 5.0) range for extracting relevant subject information such as subject HLA genotypes (the ones chosen for prediction), predicted epitope sequence and the start and end positions of the epitope sequence. *Note* - Percentile ranking in the prediction algorithm puts 100.0 as lower (bad) ranking and numbers closer to 0 are better ranking. 

```
  -p IEDB_PREDICTION, --IEDB_prediction IEDB_PREDICTION
                        Filename the text file from IEDB MHC Class-I
                        Prediction tool.
  
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
                        Output name for the text file.
  
  -lower RANKING_THRESHOLD_LOWER, --ranking_threshold_lower RANKING_THRESHOLD_LOWER
                        Cut off of epitopes at a defined percentile ranking
                        (Range 0-100). Default = 10.0
  
  -upper RANKING_THRESHOLD_UPPER, --ranking_threshold_upper RANKING_THRESHOLD_UPPER
                        Cut off of epitopes at a defined percentile ranking
                        (Range 0-100). Default = 0.0
```


2) *iedb_tool2.X.py* - The parsed predicted data is to be used together with the experimentally validated epitopes (explained below) which was also retrieved from IEDB as input for *iedb_tool.py* . Here a categorised list of ranked epitopes for which predicted epitopes were matched with experimentally validated epitopes are compiled. HLA-I typing for each subject was also matched to the experimentally validated epitopes if the HLA type information was available from IEDB (some entries in IEDB did not specify the HLA types used in experiments, and these were labelled as "Undeteremined" on IEDB). 

Parameters for use:

  ```
  -csv IEDB_CSV, --IEDB_CSV IEDB_CSV
                        Filename of compact file from IEDB.
  
  -pfile IEDB_PRED, --IEDB_PRED IEDB_PRED
                        Parsed prediction file from iedbPredictionParser.
  
  -hfile HLA_LIST, --HLA_LIST HLA_LIST
                        A file containing subject HLA information. See
                        template => SAMPLE-HLA_Set.txt
  
  -res R_LIST, --R_LIST R_LIST
                        File containing a list of HLA that can be synethesized
                        into Dextramers
  
  -hv H_VALUE, --H_VALUE H_VALUE
                        Percentage of homology between IEDB Known epitope and
                        IEDB Prediction. Default = 0.8
  
  -o OUT_DIR, --OUT_DIR OUT_DIR
                        Output directory. Otherwise files are written to
                        current directory.
  
  -extra, --EXTRA_FILES
                        Prevents removal of potentially useful output files
                        that iedb_tool.py makes.
  
  -mfile IEDB_MPRED, --IEDB_MPRED IEDB_MPRED
                        Parsed prediction file of mutated epitopes from
                        iedbPredictionParser.
  
  -mref MUT_SEQ, --MUT_SEQ MUT_SEQ
                        A fasta sequence file in amino acids with mutations
                        subbed in. Output of fixationSubber.py
  
  -n OUT_NAME, --OUT_NAME OUT_NAME
                        Prefix to add to output file name. Otherwise names
                        files with default prefix "Result"
```

---

[^1]:Melanie R Walker, Preston Leung, Elizabeth Keoshkerian, Mehdi R Pirozyan, Andrew R Lloyd, Fabio Luciani and Rowena A Bull. **Temporal dynamics of viral fitness and the adaptive immune response in HCV infection.** *eLife, 2025.* DOI:10.7554/eLife.102232.2


