These scripts allow the user to infer a distribution of fitness effects from a set of somatic mutations.

Dependencies include:
dadi pipeline https://github.com/dportik/dadi_pipeline 
fitdadi https://github.com/LohmuellerLab/fitdadi
the dadi package https://bitbucket.org/gutenkunstlab/dadi/src

While these scripts have many options and are not canonized into a pipeline, the workflow roughly follows the following logic.
1. Infer a set of betabinomial parameters for the distribution of body frequencies from the pileup with "predict_freq.py".
2. Create a likelihood-annotated vcf from the pileup with "pileup_to_annotated_vcf.py" and the betabinomial parameters. 
3. Create a "lookup" format file, a custom text format used for manual annotation of mutation types from a GTF file for your species of interest using "gtf_to_lookup.py".
4. Pass both the likelihood-annotated vcf and the lookup file to "make_mutation_frame.py", which will create and save a filtered pandas dataframe containing your final set of mutations. Any number of annotated VCFs can be passed at once to create a single large frame. Note that this expects a specific format for vcf file names- a 3 digit signifier at the front of a file name (strain-stage-identifier) followed by an underscore.
5. Pass the frame file to "calculate_dfe_from_frame.py", which will optimize gamma DFE parameters and print them out along with a lot of other information.

The "subset_pileup.py" script and "filter_somatic.py" scripts are additional scripts that may be useful for removing mutations from the pileup or the vcf step, respectively, but are not required.

These scripts are fairly specific to my current project at the moment, but much of their content can be adapted for DFE inference in other species. If your goal is just to annotate the effects of your mutations, please use the pipeline in the subdirectory "synonymity_analysis", which should serve the needs of most species without any special setup.
