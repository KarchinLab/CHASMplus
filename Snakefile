from os.path import join

# configuration file 
configfile: "chasm2/data/config.yaml"
configfile: join(config['twentyTwentyPlus'], 'config.yaml')

# include 20/20+ snakefile
include: join(config['twentyTwentyPlus'], 'Snakefile')

# parameters from command line
output_dir=config['output_dir']
mutsigcv_dir=config['mutsigcv_dir']
mutations_hg38=config['mutations_hg38']
mutations=config['mutations']

# data files
#output_dir=config["output_dir"]
snvGet="/mnt/disk003/projects/CVS-dev/SNVBox/snvGetGenomic"
data_dir=config['data']
bed=join(data_dir, config["bed"])
fasta=join(data_dir, config["fasta"])
featureList=join(data_dir, config["feature_list"])
liftoverChain=join(data_dir, config['liftoverChain'])

# parameters
hotmaps1d_windows=[0, 5, 10]
folds=range(1, 11) 
iters=range(1, 11)

#rule all:
    #input:
        #expand(join(output_dir, 'cv_trained_model/train{fold}.Rdata'), fold=folds)
        #

rule chasm2:
    input:
        null=join(output_dir, "chasm2_null_distribution.txt"),
        chasm=join(output_dir, 'chasm2_result.txt'), 
        ttplus=join(output_dir, 'output/results/r_random_forest_prediction.txt')
    output:
        join(output_dir, 'chasm2_result_final.txt')
    shell:
        "python chasm2/console/chasm.py "
        "   -c {input.chasm} "
        "   -t {input.ttplus} "
        "   -nd {input.null} "
        "   -o {output}"

#########################
# Run CHASM2 on the observed data
#########################
# convert MAF mutation file into snvbox input
rule prepSnvboxInput:
    input:
        mutations=mutations_hg38,
    params:
        mutsigcv=mutsigcv_dir
    output:
        driver=join(output_dir, "driver.snvbox_input.txt"),
        passenger=join(output_dir, "passenger.snvbox_input.txt"),
        geneFile=join(output_dir, "id2gene.txt")
    shell:
        "python chasm2/console/chasm.py prepSnvboxInput "
        "   -i {input.mutations} "
        "   -m {params.mutsigcv} "
        "   -g {output.geneFile} "
        "   -m {output.passenger} "
        "   -od {output.driver}"

# Fetch snvbox features
rule snvGetGenomic:
    input:
        driver=join(output_dir, "driver.snvbox_input.txt"),
        passenger=join(output_dir, "passenger.snvbox_input.txt"),
        geneFile=join(output_dir, "id2gene.txt")
    params:
        featureList=featureList
    output:
        join(output_dir, "snvbox_features.txt")
    shell:
        "{snvGet} "
        "   --pickone -c -f {{params.featureList}} "
        "   --o {{output}} "
        "   driver {{input.driver}} "
        "   passenger {{input.passenger}} ".format(snvGet=snvGet)

# merge additional features outside snvbox
rule mergeAdditionalFeatures:
    input:
        expand(join(output_dir, 'hotmaps1d/window{win}/result.txt'), win=hotmaps1d_windows),
        mutations=mutations_hg38,
        windowZero=join(output_dir, 'hotmaps1d/window0/result.txt'),
        windowFive=join(output_dir, 'hotmaps1d/window5/result.txt'),
        windowTen=join(output_dir, 'hotmaps1d/window10/result.txt'),
        snvbox=join(output_dir, "snvbox_features.txt"),
        geneFile=join(output_dir, "id2gene.txt")
    output:
        join(output_dir, "snvbox_features_merged.txt")
    shell:
        "python chasm2/console/chasm.py mergeFeatures "
        "   -i {input.mutations} "
        "   -ig {input.geneFile} "
        "   -s {input.snvbox} "
        "   -hm {input.windowZero},{input.windowFive},{input.windowTen} " 
        "   -o {output}"

# run hotmaps1d algorithm
rule hotmaps:
    input:
        mutations=mutations,
        bed=bed,
        fasta=fasta
    threads: 10
    params:
        window='{win}'
    output:
        result=join(output_dir, "hotmaps1d/window{win}/result.txt"),
        null=join(output_dir, "hotmaps1d/window{win}/null_distribution")
    shell:
        "probabilistic2020 hotmaps1d "
        "   -i {input.fasta} "
        "   -b {input.bed} "
        "   -m {input.mutations} "
        "   -w {params.window} "
        "   -p {threads} "
        "   -n 10000 "
        "   --report-index "
        "   -o {output.result} "
        "   -nd {output.null}"

# train a gene-hold-out cross-validated model
rule cv_train:
    input:
        features=join(output_dir, 'snvbox_features_merged.txt')
    params:
        outdir=join(output_dir, 'cv_trained_model')
    output:
        expand(join(output_dir, 'cv_trained_model/train{fold}.Rdata'), fold=folds)
    shell:
        "Rscript chasm2/r/cv_train.R -i {input.features} -t {params.outdir}"

# score mutations based on a gene-hold-out cross-validation model
# previously saved.
rule cv_test:
    input:
        expand(join(output_dir, 'cv_trained_model/train{fold}.Rdata'), fold=folds),
        features=join(output_dir, 'snvbox_features_merged.txt')
    params:
        model_dir=join(output_dir, 'cv_trained_model')
    output:
        join(output_dir, 'chasm2_result.txt') 
    shell:
        "Rscript chasm2/r/cv_test.R "
        "   -i {input.features} "
        "   -t {params.model_dir} "
        "   -o {output}"


#########################
# Handle simuated mutations
#########################
# prepare snvbox input for simulated mutations
rule simPrepSnvboxInput:
    input:
        mutations=join(output_dir, "simulated_summary/chasm_sim_maf{iter}.hg38.txt")
    output:
        driver=join(output_dir, "simulated_summary/driver_{iter,[0-9]+}.snvbox_input.txt"),
        passenger=join(output_dir, "simulated_summary/passenger_{iter,[0-9]+}.snvbox_input.txt"),
        geneFile=join(output_dir, "simulated_summary/id2gene_{iter,[0-9]+}.txt")
    shell:
        "python chasm2/console/chasm.py prepSnvboxInput "
        "   -i {input.mutations} "
        "   -g {output.geneFile} "
        "   -m {output.passenger} "
        "   -od {output.driver}"

# Fetch snvbox features for simulations
rule simSnvGetGenomic:
    input:
        driver=join(output_dir, "simulated_summary/driver_{iter}.snvbox_input.txt"),
        passenger=join(output_dir, "simulated_summary/passenger_{iter}.snvbox_input.txt"),
        geneFile=join(output_dir, "simulated_summary/id2gene_{iter}.txt")
    params:
        featureList=featureList
    output:
        join(output_dir, "simulated_summary/snvbox_features_{iter,[0-9]+}.txt")
    shell:
        "{snvGet} "
        "   --pickone -c -f {{params.featureList}} "
        "   --o {{output}} "
        "   driver {{input.driver}} "
        "   passenger {{input.passenger}} ".format(snvGet=snvGet)

# merge additional features outside snvbox
rule simMergeAdditionalFeatures:
    input:
        expand(join(output_dir, 'simulated_summary/hotmaps1d/result_window{win_sim}_{{iter}}.txt'), win_sim=hotmaps1d_windows),
        mutations=mutations_hg38,
        windowZero=join(output_dir, 'simulated_summary/hotmaps1d/result_window0_{iter}.txt'),
        windowFive=join(output_dir, 'simulated_summary/hotmaps1d/result_window5_{iter}.txt'),
        windowTen=join(output_dir, 'simulated_summary/hotmaps1d/result_window10_{iter}.txt'),
        snvbox=join(output_dir, "simulated_summary/snvbox_features_{iter}.txt"),
        geneFile=join(output_dir, "simulated_summary/id2gene_{iter}.txt")
    output:
        join(output_dir, "simulated_summary/snvbox_features_merged_{iter,[0-9]+}.txt")
    shell:
        "python chasm2/console/chasm.py mergeFeatures "
        "   -i {input.mutations} "
        "   -ig {input.geneFile} "
        "   -s {input.snvbox} "
        "   -hm {input.windowZero},{input.windowFive},{input.windowTen} " 
        "   -o {output}"

# run hotmaps1d on the simulated mutations
rule simHotmaps:
    input:
        mutations=join(output_dir, 'simulated_summary/chasm_sim_maf{iter}.txt'),
        bed=bed,
        fasta=fasta
    params:
        window='{win_sim}'
    output:
        result=join(output_dir, 'simulated_summary/hotmaps1d/result_window{win_sim,[0-9]+}_{iter,[0-9]+}.txt')
    shell:
        "probabilistic2020 hotmaps1d "
        "   -i {input.fasta} "
        "   -b {input.bed} "
        "   -m {input.mutations} "
        "   -w {params.window} "
        "   -p {threads} "
        "   -n 10000 "
        "   --report-index "
        "   -o {output.result} "

rule simCvTest:
    input: 
        expand(join(output_dir, 'cv_trained_model/train{fold}.Rdata'), fold=folds),
        features=join(output_dir, 'simulated_summary/snvbox_features_merged_{iter}.txt')
    params:
        model_dir=join(output_dir, 'cv_trained_model')
    output:
        join(output_dir, 'simulated_summary/chasm2_result{iter,[0-9]+}.txt') 
    shell:
        "Rscript chasm2/r/cv_test.R "
        "   -i {input.features} "
        "   -t {params.model_dir} "
        "   -o {output}"


############################
# exclusively only predict 
# on indicated mutations
# for 20/20+
############################
rule predict_ttplus_only:
    input:
        trained_classifier=join(output_dir, "trained.Rdata"),
        sim_features=join(output_dir, "simulated_summary/simulated_features{iter}.txt"),
    params:
        outdir=join(output_dir, "simulated_summary/2020plus/sim{iter,[0-9]+}")
    output: 
        join(output_dir, "simulated_summary/2020plus/sim{iter,[0-9]+}/results/r_random_forest_prediction.txt")
    shell:
        "python `which 2020plus.py` --log-level=INFO "
        "   --out-dir {params.outdir} "
        "   classify "
        "   --trained-classifier {input.trained_classifier} "
        "   --features {input.sim_features}"


##########################
# Create null distribution
##########################
rule nullDistribution:
    input:
        expand(join(output_dir, "simulated_summary/2020plus/sim{iter}/results/r_random_forest_prediction.txt"), iter=iters),
        expand(join(output_dir, "simulated_summary/chasm2_result{iter}.txt"), iter=iters)
    params:
        sim_dir=join(output_dir, "simulated_summary")
    output:
        join(output_dir, "chasm2_null_distribution.txt")
    shell:
        "python chasm2/console/chasm.py nullDistribution "
        "   -s {params.sim_dir} "
        "   -o {output}"


#########################
# Handle conversion to hg38 in simulations
#########################
# create input for liftover
rule simMafToBed:
    input:
        mutations=join(output_dir, 'simulated_summary/chasm_sim_maf{iter}.txt'),
    output:
        bed_out=join(output_dir, 'simulated_summary/chasm_sim_maf{iter,[0-9]+}.bed')
    shell:
        "python scripts/maf_to_bed.py "
        "   -i {input.mutations} "
        "   -o {output.bed_out}"

# liftover coordinates of simulations
rule simLiftover:
    input:
        bed_in=join(output_dir, 'simulated_summary/chasm_sim_maf{iter}.bed'),
        chain=liftoverChain
    output:
        bed_hg38=join(output_dir, 'simulated_summary/chasm_sim_maf_hg38_{iter,[0-9]+}.bed'),
        unmapped=join(output_dir, 'simulated_summary/chasm_sim_maf{iter,[0-9]+}.unmapped.bed')
    shell:
        "liftOver {input.bed_in} {input.chain} {output.bed_hg38} {output.unmapped}"

# Create simulated MAF file with hg38 coordinates
rule simLiftoverMaf:
    input:
        bed_hg38=join(output_dir, 'simulated_summary/chasm_sim_maf_hg38_{iter}.bed'),
        mutations=join(output_dir, 'simulated_summary/chasm_sim_maf{iter}.txt')
    output:
        mutations_hg38=join(output_dir, 'simulated_summary/chasm_sim_maf{iter,[0-9]+}.hg38.txt')
    shell:
        "python scripts/create_hg38_maf.py "
        "   -i {input.bed_hg38} "
        "   -m {input.mutations} "
        "   -o {output.mutations_hg38}"
