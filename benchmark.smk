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

rule snvboxTxToGenomic:
    input:
        berger_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al.snvbox_tx.txt'),
        berger_egfr_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.snvbox_tx.txt'),
        kim_sbox=join(benchmark_dir, 'snvbox_input/kim_et_al.snvbox_tx.txt'),
        iarc_tp53_sbox=join(benchmark_dir, 'snvbox_input/iarc_tp53.snvbox_tx.txt')
    params:
        txToGenome='python2.7 {0}'.format(prot_to_gen)
    output:
        berger_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al.snvbox_genomic.txt'),
        berger_egfr_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.snvbox_genomic.txt'),
        kim_sbox=join(benchmark_dir, 'snvbox_input/kim_et_al.snvbox_genomic.txt'),
        iarc_tp53_sbox=join(benchmark_dir, 'snvbox_input/iarc_tp53.snvbox_genomic.txt')
    shell:
        """
        {params.txToGenome} -i {input.berger_sbox} -o {output.berger_sbox}
        {params.txToGenome} -i {input.berger_egfr_sbox} -o {output.berger_egfr_sbox}
        {params.txToGenome} -i {input.kim_sbox} -o {output.kim_sbox}
        {params.txToGenome} -i {input.iarc_tp53_sbox} -o {output.iarc_tp53_sbox}
        """

# create BED file for liftover from hg38 to hg19
rule snvboxToBed:
    input:
        berger_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al.snvbox_genomic.txt'),
        berger_egfr_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.snvbox_genomic.txt'),
        kim_sbox=join(benchmark_dir, 'snvbox_input/kim_et_al.snvbox_genomic.txt'),
        iarc_tp53_sbox=join(benchmark_dir, 'snvbox_input/iarc_tp53.snvbox_genomic.txt')
    output:
        berger_bed=join(benchmark_dir, 'snvbox_input/berger_et_al.hg38.bed'),
        berger_egfr_bed=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.hg38.bed'),
        kim_bed=join(benchmark_dir, 'snvbox_input/kim_et_al.hg38.bed'),
        iarc_tp53_bed=join(benchmark_dir, 'snvbox_input/iarc_tp53.hg38.bed')
    shell:
        """
        python scripts/benchmark/snvbox_to_bed.py -i {input.berger_sbox} -o {output.berger_bed}
        python scripts/benchmark/snvbox_to_bed.py -i {input.berger_egfr_sbox} -o {output.berger_egfr_bed}
        python scripts/benchmark/snvbox_to_bed.py -i {input.kim_sbox} -o {output.kim_bed}
        python scripts/benchmark/snvbox_to_bed.py -i {input.iarc_tp53_sbox} -o {output.iarc_tp53_bed}
        """

# liftover from hg38 bed file to hg19
rule liftOver:
    input:
        berger_bed=join(benchmark_dir, 'snvbox_input/berger_et_al.hg38.bed'),
        berger_egfr_bed=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.hg38.bed'),
        kim_bed=join(benchmark_dir, 'snvbox_input/kim_et_al.hg38.bed'),
        iarc_tp53_bed=join(benchmark_dir, 'snvbox_input/iarc_tp53.hg38.bed')
    params:
        chain='chasm2/data/hg38ToHg19.over.chain.gz'
    output:
        berger_bed=join(benchmark_dir, 'snvbox_input/berger_et_al.hg19.bed'),
        berger_egfr_bed=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.hg19.bed'),
        kim_bed=join(benchmark_dir, 'snvbox_input/kim_et_al.hg19.bed'),
        iarc_tp53_bed=join(benchmark_dir, 'snvbox_input/iarc_tp53.hg19.bed'),
        berger_unmapped=join(benchmark_dir, 'snvbox_input/berger_et_al.unmapped.bed'),
        berger_egfr_unmapped=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.unmapped.bed'),
        kim_unmapped=join(benchmark_dir, 'snvbox_input/kim_et_al.unmapped.bed'),
        iarc_tp53_unmapped=join(benchmark_dir, 'snvbox_input/iarc_tp53.unmapped.bed')
    shell:
       """
       liftOver {input.berger_bed} {params.chain} {output.berger_bed} {output.berger_unmapped}
       liftOver {input.berger_egfr_bed} {params.chain} {output.berger_egfr_bed} {output.berger_egfr_unmapped}
       liftOver {input.kim_bed} {params.chain} {output.kim_bed} {output.kim_unmapped}
       liftOver {input.iarc_tp53_bed} {params.chain} {output.iarc_tp53_bed} {output.iarc_tp53_unmapped}
       """

rule bedToSnvbox:
    input:
        berger_bed=join(benchmark_dir, 'snvbox_input/berger_et_al.hg19.bed'),
        berger_egfr_bed=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.hg19.bed'),
        kim_bed=join(benchmark_dir, 'snvbox_input/kim_et_al.hg19.bed'),
        iarc_tp53_bed=join(benchmark_dir, 'snvbox_input/iarc_tp53.hg19.bed'),
        berger_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al.snvbox_genomic.txt'),
        berger_egfr_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.snvbox_genomic.txt'),
        kim_sbox=join(benchmark_dir, 'snvbox_input/kim_et_al.snvbox_genomic.txt'),
        iarc_tp53_sbox=join(benchmark_dir, 'snvbox_input/iarc_tp53.snvbox_genomic.txt')
    output:
        berger_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al.snvbox_genomic.hg19.txt'),
        berger_egfr_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.snvbox_genomic.hg19.txt'),
        kim_sbox=join(benchmark_dir, 'snvbox_input/kim_et_al.snvbox_genomic.hg19.txt'),
        iarc_tp53_sbox=join(benchmark_dir, 'snvbox_input/iarc_tp53.snvbox_genomic.hg19.txt')
    shell:
        """
        python scripts/benchmark/bed_to_snvbox.py -i {input.berger_sbox} -b {input.berger_bed} -o {output.berger_sbox}
        python scripts/benchmark/bed_to_snvbox.py -i {input.berger_egfr_sbox} -b {input.berger_egfr_bed} -o {output.berger_egfr_sbox}
        python scripts/benchmark/bed_to_snvbox.py -i {input.kim_sbox} -b {input.kim_bed} -o {output.kim_sbox}
        python scripts/benchmark/bed_to_snvbox.py -i {input.iarc_tp53_sbox} -b {input.iarc_tp53_bed} -o {output.iarc_tp53_sbox}
        """
