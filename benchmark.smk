from os.path import join
configfile: "chasm2/data/config.yaml"

# important variables
mutations=config['mutations']
benchmk=config['benchmark']
snvbox=config['snvbox']
hotmaps_null=config['hotmaps_null']
benchmark_dir=config['benchmark_dir']
trained_chasm2=config['trained_model']

# snvget 
prot_to_gen="/mnt/disk003/projects/CVS-dev/SNVBox/proteinToGenomic.py"

#include: 'Snakefile'

# doesn't compute features based on provided data
rule chasm2_benchmark:
    input:
        features=join(benchmark_dir, 'snvbox_features_merged.txt')
    params:
        model_dir=trained_chasm2
    output:
        join(benchmark_dir, 'chasm2_result_pretrained.txt') 
    shell:
        "Rscript chasm2/r/cv_test.R "
        "   -i {input.features} "
        "   -t {params.model_dir} "
        "   -o {output}"

# merge additional features outside snvbox
# for benchmark assesment
rule mergeAdditionalFeaturesBenchmark:
    input:
        hotmaps=hotmaps_null,
        snvbox=snvbox,
        benchmark=benchmk,
        mutations=mutations
    output:
        join(benchmark_dir, 'snvbox_features_merged.txt')
    shell:
        "python chasm2/console/chasm.py mergeBenchmarkFeatures "
        "   -m {input.mutations} "
        "   -s {input.snvbox} "
        "   -w 0,5,10 "
        "   -b {input.benchmark} "
        "   -n {input.hotmaps} " 
        "   -o {output}"
