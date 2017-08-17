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

rule prepSnvboxTxInput:
    input:
        # mutation files
        berger=join(benchmark_dir, 'berger_et_al.txt'),
        berger_egfr=join(benchmark_dir, 'berger_et_al_egfr.txt'),
        kim=join(benchmark_dir, 'kim_et_al.txt'),
        iarc_tp53=join(benchmark_dir, 'iarc_tp53.txt'),
        # preferred transcript files
        berger_tx=join(benchmark_dir, 'preferred_tx/berger_et_al.preferred_tx.txt'),
        berger_egfr_tx=join(benchmark_dir, 'preferred_tx/berger_et_al_egfr.preferred_tx.txt'),
        kim_tx=join(benchmark_dir, 'preferred_tx/kim_et_al.preferred_tx.txt'),
        iarc_tp53_tx=join(benchmark_dir, 'preferred_tx/iarc_tp53.preferred_tx.txt')
    output:
        berger_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al.snvbox_tx.txt'),
        berger_egfr_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.snvbox_tx.txt'),
        kim_sbox=join(benchmark_dir, 'snvbox_input/kim_et_al.snvbox_tx.txt'),
        iarc_tp53_sbox=join(benchmark_dir, 'snvbox_input/iarc_tp53.snvbox_tx.txt')
    shell:
        """
        python scripts/benchmark/create_kim_snvbox_input.py -i {input.kim} -p {input.kim_tx} -o {output.kim_sbox}
        python scripts/benchmark/create_berger_snvbox_input.py -i {input.berger} -p {input.berger_tx} -o {output.berger_sbox}
        python scripts/benchmark/create_berger_snvbox_input.py -i {input.berger_egfr} -p {input.berger_egfr_tx} -o {output.berger_egfr_sbox}
        python scripts/benchmark/create_tp53_snvbox_input.py -i {input.iarc_tp53} -p {input.iarc_tp53_tx} -o {output.iarc_tp53_sbox}
        """
