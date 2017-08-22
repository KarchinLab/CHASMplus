from os.path import join, abspath
configfile: "chasm2/data/config.yaml"

# important variables
mutations=config['mutations']
#benchmk=config['benchmark']
#snvbox=config['snvbox']
hotmaps_null=config['hotmaps_null']
benchmark_dir=config['benchmark_dir']
trained_chasm2=config['trained_model']

# snvget 
prot_to_gen="/mnt/disk003/projects/CVS-dev/SNVBox/proteinToGenomic.py"
snvgetgenomic="/mnt/disk003/projects/CVS-dev/SNVBox/snvGetGenomic"
snvgettranscript="/mnt/disk003/projects/CVS-dev/SNVBox/snvGetTranscript"

# original chasm
orig_chasm_dir="/mnt/disk003/projects/CVS-dev/CHASM/"

# other methods
candra="methods/CanDrA.v1.0/open_candra.pl"
annovar="methods/annovar/table_annovar.pl"
transfic="methods/transfic/bin/transf_scores.pl"
parssnp="methods/ParsSNP/ParsSNP_application.r"

# outher config
data_dir=config['chasm2_data']
featureList=join(data_dir, config["feature_list"])

# list of benchmarks to run
mybenchmarks=['berger_et_al', 'berger_et_al_egfr', 'kim_et_al', 'iarc_tp53']
mybenchmarks_maf=['patrick_et_al', 'msk_impact', 'mc3']

rule perform_benchmark:
    input:
        # regular benchmarks
        expand(join(benchmark_dir, 'methods/output/{benchmark}.chasm2_output.txt'), benchmark=mybenchmarks),
        expand(join(benchmark_dir, 'methods/input/{benchmark}.fathmm_input.txt'), benchmark=mybenchmarks),
        expand(join(benchmark_dir, "methods/output/{benchmark}.candra_output.txt"), benchmark=mybenchmarks),
        expand(join(benchmark_dir, 'methods/output/{benchmark}.annovar_output.hg19_multianno.txt'), benchmark=mybenchmarks),
        expand(join(benchmark_dir, 'methods/output/{benchmark}.transfic_output.txt'), benchmark=mybenchmarks),
        expand(abspath(join(benchmark_dir, 'methods/output/ParsSNP.output.{benchmark}.annovar_output.hg19_multianno.txt')), benchmark=mybenchmarks),
        expand(join(benchmark_dir, 'methods/output/{benchmark}.chasm_output.txt'), benchmark=mybenchmarks)

rule perform_benchmark_maflike:
    input:
        # benchmarks from MAF like files
        expand(join(benchmark_dir, 'methods/output/{benchmark_maf}.chasm2_output.txt'), benchmark_maf=mybenchmarks_maf[:-1]),
        expand(join(benchmark_dir, "methods/output/{benchmark_maf}.candra_output.txt"), benchmark_maf=mybenchmarks_maf),
        expand(abspath(join(benchmark_dir, 'methods/output/{benchmark_maf}.annovar_output.hg19_multianno.txt')), benchmark_maf=mybenchmarks_maf),
        expand(join(benchmark_dir, 'methods/output/{benchmark_maf}.transfic_output.txt'), benchmark_maf=mybenchmarks_maf),
        expand(abspath(join(benchmark_dir, 'methods/output/ParsSNP.output.{benchmark_maf}.annovar_output.hg19_multianno.txt')), benchmark_maf=mybenchmarks_maf),
        expand(join(benchmark_dir, 'methods/output/{benchmark_maf}.chasm_output.txt'), benchmark_maf=mybenchmarks_maf)

####################
# CHASM2
####################
# doesn't compute features based on provided data
rule chasm2_benchmark:
    input:
        features=join(benchmark_dir, 'snvbox_output/features_{benchmark}_merged.txt')
    params:
        model_dir=trained_chasm2
    output:
        join(benchmark_dir, 'methods/output/{benchmark}.chasm2_output.txt')
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
        myinput=join(benchmark_dir, '{benchmark}.txt'),
        snvbox=join(benchmark_dir, 'snvbox_output/features_{benchmark}.txt'),
        mutations=mutations
    output:
        join(benchmark_dir, 'snvbox_output/features_{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}_merged.txt')
    shell:
        "python chasm2/console/chasm.py mergeBenchmarkFeatures "
        "   -m {input.mutations} "
        "   -s {input.snvbox} "
        "   -w 0,5,10 "
        "   -b {input.myinput} "
        "   -n {input.hotmaps} " 
        "   -o {output}"

rule mergeAdditionalFeaturesBenchmarkMaf:
    input:
        hotmaps=hotmaps_null,
        myinput=join(benchmark_dir, '{benchmark_maf}.maf'),
        snvbox=join(benchmark_dir, 'snvbox_output/features_{benchmark_maf}.txt'),
        mutations=mutations
    output:
        join(benchmark_dir, 'snvbox_output/features_{benchmark_maf,patrick_et_al|msk_impact}_merged.txt')
    shell:
        "python chasm2/console/chasm.py mergeBenchmarkFeatures "
        "   -m {input.mutations} "
        "   -s {input.snvbox} "
        "   -w 0,5,10 "
        "   -b {input.myinput} "
        "   -n {input.hotmaps} " 
        "   -o {output}"

########################
# Convert maf to hg38
########################
# create input for liftover for observed data
rule mafToBed:
    input:
        mutations=join(benchmark_dir, '{benchmark_maf}.maf')
    output:
        bed_out=join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact}.bed')
    shell:
        "python scripts/maf_to_bed.py "
        "   -i {input.mutations} "
        "   -o {output.bed_out}"

# liftover coordinates of observed data
rule liftover:
    input:
        bed_in=join(benchmark_dir, 'snvbox_input/{benchmark_maf}.bed'),
        chain='chasm2/data/hg19ToHg38.over.chain.gz'
    output:
        bed_hg38=join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact}_hg38.bed'),
        unmapped=join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact}.unmapped.bed')
    shell:
        "liftOver {input.bed_in} {input.chain} {output.bed_hg38} {output.unmapped}"

# Create simulated MAF file with hg38 coordinates
rule liftoverMaf:
    input:
        bed_hg38=join(benchmark_dir, 'snvbox_input/{benchmark_maf}_hg38.bed'),
        mutations=join(benchmark_dir, '{benchmark_maf}.maf')
    output:
        mutations_hg38=join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact}.hg38.maf')
    shell:
        "python scripts/create_hg38_maf.py "
        "   -i {input.bed_hg38} "
        "   -m {input.mutations} "
        "   -o {output.mutations_hg38}"

#################
# Query SNVBox with genomic from MAF
#################
# convert MAF mutation file into snvbox input
rule prepSnvboxInput:
    input:
        mutations=join(benchmark_dir, 'snvbox_input/{benchmark_maf}.hg38.maf'),
    output:
        driver=join(benchmark_dir, "snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_driver.snvbox_input.txt"),
        passenger=join(benchmark_dir, "snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_passenger.snvbox_input.txt"),
        geneFile=join(benchmark_dir, "snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_id2gene.txt")
    shell:
        "python chasm2/console/chasm.py prepSnvboxInput "
        "   -i {input.mutations} "
        "   -g {output.geneFile} "
        "   -op {output.passenger} "
        "   -od {output.driver}"

# Fetch snvbox features
rule snvGetGenomic:
    input:
        driver=join(benchmark_dir, "snvbox_input/{benchmark_maf}_driver.snvbox_input.txt"),
        passenger=join(benchmark_dir, "snvbox_input/{benchmark_maf}_passenger.snvbox_input.txt"),
        geneFile=join(benchmark_dir, "snvbox_input/{benchmark_maf}_id2gene.txt")
    params:
        featureList=featureList
    output:
        join(benchmark_dir, "snvbox_output/features_{benchmark_maf,patrick_et_al|msk_impact|mc3}.txt")
    shell:
        "python2.7 {snvgetgenomic} "
        "   --pickone -c -f {{params.featureList}} "
        "   --o {{output}} "
        "   driver {{input.driver}} "
        "   passenger {{input.passenger}} ".format(snvgetgenomic=snvgetgenomic)

#########################
# Handle SNVBox
#########################
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

# get features from snvbox
rule snvGetTranscript:
    input:
        berger_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al.snvbox_tx.txt'),
        berger_egfr_sbox=join(benchmark_dir, 'snvbox_input/berger_et_al_egfr.snvbox_tx.txt'),
        kim_sbox=join(benchmark_dir, 'snvbox_input/kim_et_al.snvbox_tx.txt'),
        iarc_tp53_sbox=join(benchmark_dir, 'snvbox_input/iarc_tp53.snvbox_tx.txt')
    params:
        snvget="python2.7 {0}".format(snvgettranscript),
        featlist="chasm2/data/Features.list"
    output:
        berger_feat=join(benchmark_dir, 'snvbox_output/features_berger_et_al.txt'),
        berger_egfr_feat=join(benchmark_dir, 'snvbox_output/features_berger_et_al_egfr.txt'),
        kim_feat=join(benchmark_dir, 'snvbox_output/features_kim_et_al.txt'),
        iarc_tp53_feat=join(benchmark_dir, 'snvbox_output/features_iarc_tp53.txt'),
    shell:
        """
        {params.snvget} -c -f {params.featlist} -o {output.berger_feat} NA {input.berger_sbox}
        {params.snvget} -c -f {params.featlist} -o {output.berger_egfr_feat} NA {input.berger_egfr_sbox}
        {params.snvget} -c -f {params.featlist} -o {output.kim_feat} NA {input.kim_sbox}
        {params.snvget} -c -f {params.featlist} -o {output.iarc_tp53_feat} NA {input.iarc_tp53_sbox}
        """


###################
# candra
###################
# prep input files
rule prepCandraInput:
    input:
        join(benchmark_dir, 'snvbox_input/{benchmark}.snvbox_genomic.hg19.txt')
    output:
        join(benchmark_dir, 'methods/input/{benchmark}.candra_input.txt')
    shell:
        "python scripts/benchmark/snvbox2candra.py -i {input} -o {output}"

rule prepCandraInputMaf:
    input:
        join(benchmark_dir, '{benchmark_maf}.maf')
    output:
        join(benchmark_dir, 'methods/input/{benchmark_maf,patrick_et_al|msk_impact|mc3}.candra_input.txt')
    shell:
        "python scripts/benchmark/maf_to_candra.py -i {input} -o {output}"

# run candra
rule runCandra:
    input:
        join(benchmark_dir, "methods/input/{benchmark}.candra_input.txt")
    params:
        candra="perl {0}".format(candra)
    output:
        join(benchmark_dir, "methods/output/{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}.candra_output.txt")
    shell:
        "{params.candra} OVC {input} > {output}"

rule runCandraMaf:
    input:
        join(benchmark_dir, "methods/input/{benchmark_maf}.candra_input.txt")
    params:
        candra="perl {0}".format(candra)
    output:
        join(benchmark_dir, "methods/output/{benchmark_maf,patrick_et_al|msk_impact|mc3}.candra_output.txt")
    shell:
        "{params.candra} OVC {input} > {output}"

##################
# annovar
##################
rule prepAnnovarInput:
    input:
        join(benchmark_dir, 'snvbox_input/{benchmark}.snvbox_genomic.hg19.txt')
    output:
        join(benchmark_dir, 'methods/input/{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}.annovar_input.txt')
    shell:
        "python scripts/benchmark/snvbox2annovar.py "
        "   -i {input} "
        "   -o {output} "

rule prepAnnovarInputMaf:
    input:
        join(benchmark_dir, '{benchmark_maf}.maf')
    output:
        join(benchmark_dir, 'methods/input/{benchmark_maf,patrick_et_al|msk_impact|mc3}.annovar_input.txt')
    shell:
        "python scripts/benchmark/maf_to_annovar.py "
        "   -i {input} "
        "   -o {output} "

rule runAnnovar:
    input:
        join(benchmark_dir, 'methods/input/{benchmark}.annovar_input.txt')
    params:
        annovar='perl {0}'.format(annovar),
        prefix=join(benchmark_dir, 'methods/output/{benchmark}.annovar_output')
    output:
        join(benchmark_dir, 'methods/output/{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}.annovar_output.hg19_multianno.txt')
    shell:
        "{params.annovar} {input} -buildver hg19 methods/annovar/humandb -out {params.prefix} -remove -protocol refGene,ensGene,ljb26_all,revel,mcap -operation g,g,f,f,f -nastring NA"

rule runAnnovarMaf:
    input:
        join(benchmark_dir, 'methods/input/{benchmark_maf}.annovar_input.txt')
    params:
        annovar='perl {0}'.format(annovar),
        prefix=join(benchmark_dir, 'methods/output/{benchmark_maf}.annovar_output')
    output:
        abspath(join(benchmark_dir, 'methods/output/{benchmark_maf,patrick_et_al|msk_impact|mc3}.annovar_output.hg19_multianno.txt'))
    shell:
        "{params.annovar} {input} -buildver hg19 methods/annovar/humandb -out {params.prefix} -remove -protocol refGene,ensGene,ljb26_all,revel,mcap -operation g,g,f,f,f -nastring NA"

##################
# fathmm
##################
rule prepFathmmInput:
    input:
        infile=join(benchmark_dir, 'snvbox_input/{benchmark}.snvbox_tx.txt')
    params:
        user=config['mysql_user'],
        passwd=config['mysql_passwd']
    output:
        join(benchmark_dir, 'methods/input/{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}.fathmm_input.txt')
    shell:
        "python scripts/benchmark/snvbox2fathmm.py -i {input} -o {output}"
        "   -i {input.infile} "
        "   --mysql-user {params.user} "
        "   --mysql-passwd {params.passwd} "
        "   -o {output} "

#####################
# transfic
#####################
rule prepTransficInput:
    input:
        join(benchmark_dir, 'methods/output/{benchmark}.annovar_output.hg19_multianno.txt')
    output:
        join(benchmark_dir, 'methods/input/{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}.transfic_input.txt')
    shell:
        "python scripts/benchmark/annovar2transfic.py "
        "   -i {input} "
        "   -o {output} "

rule prepTransficInputMaf:
    input:
        abspath(join(benchmark_dir, 'methods/output/{benchmark_maf}.annovar_output.hg19_multianno.txt'))
    output:
        join(benchmark_dir, 'methods/input/{benchmark_maf,patrick_et_al|msk_impact|mc3}.transfic_input.txt')
    shell:
        "python scripts/benchmark/annovar2transfic.py "
        "   -i {input} "
        "   -o {output} "

# run transfic
rule runTransfic:
    input:
        join(benchmark_dir, 'methods/input/{benchmark}.transfic_input.txt')
    params:
        transfic='perl {0}'.format(transfic)
    output:
        join(benchmark_dir, 'methods/output/{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}.transfic_output.txt')
    shell:
        "{params.transfic} gosmf {input} {output}"

rule runTransficMaf:
    input:
        join(benchmark_dir, 'methods/input/{benchmark_maf}.transfic_input.txt')
    params:
        transfic='perl {0}'.format(transfic)
    output:
        join(benchmark_dir, 'methods/output/{benchmark_maf,patrick_et_al|msk_impact|mc3}.transfic_output.txt')
    shell:
        "{params.transfic} gosmf {input} {output}"

######################
# ParsSNP
######################
rule runParssnp:
    input:
        abspath(join(benchmark_dir, 'methods/output/{benchmark}.annovar_output.hg19_multianno.txt'))
    params:
        parssnp='Rscript {0}'.format(parssnp)
    output:
        abspath(join(benchmark_dir, 'methods/output/ParsSNP.output.{benchmark,^(?!.*(patrick_et_al|msk_impact|mc3)).*$}.annovar_output.hg19_multianno.txt'))
    shell:
        "{params.parssnp} {input}"

rule runParssnpMaf:
    input:
        abspath(join(benchmark_dir, 'methods/output/{benchmark_maf}.annovar_output.hg19_multianno.txt'))
    params:
        parssnp='Rscript {0}'.format(parssnp)
    output:
        abspath(join(benchmark_dir, 'methods/output/ParsSNP.output.{benchmark_maf,patrick_et_al|msk_impact|mc3}.annovar_output.hg19_multianno.txt'))
    shell:
        "{params.parssnp} {input}"

######################
# Original CHASM
######################

## CHASM for transcript format ##
rule runChasmTranscript:
    input:
        expand(abspath(join(benchmark_dir, 'snvbox_input/{{benchmark}}_chasm_input{citer}driver{citer}.output')), citer=range(10)),
        abspath(join(benchmark_dir, 'snvbox_input/{benchmark}_chasm_input_passengerdriver0.output'))
    params:
        input_prefix=join(benchmark_dir, 'snvbox_input/{benchmark}')
    output:
        join(benchmark_dir, 'methods/output/{benchmark}.chasm_output.txt')
    shell:
        "python scripts/benchmark/merge_chasm_result.py "
        "   -i {params.input_prefix} "
        "   -o {output}"

rule runChasmTranscriptDriver:
    input:
        abspath(join(benchmark_dir, 'snvbox_input/{benchmark}_chasm_input{citer}.txt')),
    params:
        work_dir=join(benchmark_dir, 'snvbox_input'),
        orig_chasm_dir=orig_chasm_dir,
        chasm_model='driver{citer}'
    output:
        output_file=abspath(join(benchmark_dir, 'snvbox_input/{benchmark}_chasm_input{citer}driver{citer}.output'))
    shell:
        """
        cd {params.orig_chasm_dir}
        python2.7 RunChasm {params.chasm_model} {input} 
        cd -
        """

rule runChasmTranscriptPassenger:
    input:
        abspath(join(benchmark_dir, 'snvbox_input/{benchmark}_chasm_input_passenger.txt')),
    params:
        work_dir=join(benchmark_dir, 'snvbox_input'),
        orig_chasm_dir=orig_chasm_dir,
        chasm_model='driver0'
    output:
        output_file=abspath(join(benchmark_dir, 'snvbox_input/{benchmark}_chasm_input_passengerdriver0.output'))
    shell:
        """
        cd {params.orig_chasm_dir}
        python2.7 RunChasm {params.chasm_model} {input} 
        """

rule prepChasmTxInput:
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
    params:
        outdir=join(benchmark_dir, 'snvbox_input'),
        chasm_dir=orig_chasm_dir
    output:
        expand(abspath(join(benchmark_dir, 'snvbox_input/berger_et_al_chasm_input{num}.txt')),num=range(10)),
        expand(abspath(join(benchmark_dir, 'snvbox_input/berger_et_al_egfr_chasm_input{num}.txt')),num=range(10)),
        expand(abspath(join(benchmark_dir, 'snvbox_input/kim_et_al_chasm_input{num}.txt')),num=range(10)),
        expand(abspath(join(benchmark_dir, 'snvbox_input/iarc_tp53_chasm_input{num}.txt')),num=range(10)),
        abspath(join(benchmark_dir, 'snvbox_input/berger_et_al_chasm_input_passenger.txt')),
        abspath(join(benchmark_dir, 'snvbox_input/berger_et_al_egfr_chasm_input_passenger.txt')),
        abspath(join(benchmark_dir, 'snvbox_input/kim_et_al_chasm_input_passenger.txt')),
        abspath(join(benchmark_dir, 'snvbox_input/iarc_tp53_chasm_input_passenger.txt'))
    shell:
        """
        python scripts/benchmark/create_kim_snvbox_input.py -i {input.kim} -p {input.kim_tx} -pf kim_et_al -c {params.chasm_dir}/drivers_gene_crossval -o {params.outdir}
        python scripts/benchmark/create_berger_snvbox_input.py -i {input.berger} -p {input.berger_tx} -pf berger_et_al -c {params.chasm_dir}/drivers_gene_crossval -o {params.outdir}
        python scripts/benchmark/create_berger_snvbox_input.py -i {input.berger_egfr} -p {input.berger_egfr_tx} -pf berger_et_al_egfr -c {params.chasm_dir}/drivers_gene_crossval -o {params.outdir}
        python scripts/benchmark/create_tp53_snvbox_input.py -i {input.iarc_tp53} -p {input.iarc_tp53_tx} -pf iarc_tp53 -c {params.chasm_dir}/drivers_gene_crossval -o {params.outdir}
        """

## CHASM for MAF files ##
rule runChasmMaf:
    input:
        expand(abspath(join(benchmark_dir, 'snvbox_input/{{benchmark_maf}}_chasm_input{chasm_iter}driver{chasm_iter}.output')), chasm_iter=range(10)),
        abspath(join(benchmark_dir, 'snvbox_input/{benchmark_maf}_chasm_input_passengerdriver0.output'))
    params:
        input_prefix=join(benchmark_dir, 'snvbox_input/{benchmark_maf}')
    output:
        join(benchmark_dir, 'methods/output/{benchmark_maf,patrick_et_al|msk_impact|mc3}.chasm_output.txt')
    shell:
        "python scripts/benchmark/merge_chasm_result.py "
        "   -i {params.input_prefix} "
        "   --use-maf "
        "   -o {output}"

rule runChasmMafDriver:
    input:
        join(benchmark_dir, 'snvbox_input/{benchmark_maf}_split_driver{chasm_iter}.txt'),
    params:
        work_dir=join(benchmark_dir, 'snvbox_input'),
        orig_chasm_dir=orig_chasm_dir,
        chasm_model='driver{chasm_iter}'
    output:
        idfile=join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_id2gene_driver{chasm_iter}.txt'),
        snvbox_input=abspath(join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_chasm_input{chasm_iter}.txt')),
        output_file=abspath(join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_chasm_input{chasm_iter}driver{chasm_iter}.output'))
    shell:
        """
        python chasm2/python/snvbox/make_snvbox_input.py \
            -i {input} \
            -onco chasm2/data/oncogenes.txt -tsg chasm2/data/tsgs.txt  \
            -g {output.idfile} \
            -od {output.snvbox_input}
        cd {params.orig_chasm_dir} 
        python2.7 RunChasm {params.chasm_model} {output.snvbox_input} -g --pickone
        cd -
        """

rule runChasmMafPassenger:
    input:
        join(benchmark_dir, 'snvbox_input/{benchmark_maf}_split_passenger.txt'),
    params:
        work_dir=join(benchmark_dir, 'snvbox_input'),
        orig_chasm_dir=orig_chasm_dir,
        chasm_model='driver0'
    output:
        idfile=join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_id2gene_passenger.txt'),
        snvbox_input=abspath(join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_chasm_input_passenger.txt')),
        output_file=abspath(join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_chasm_input_passengerdriver0.output'))
    shell:
        """
        python chasm2/python/snvbox/make_snvbox_input.py \
            -i {input} \
            -onco chasm2/data/oncogenes.txt -tsg chasm2/data/tsgs.txt  \
            -g {output.idfile} \
            -od {output.snvbox_input}
        cd {params.orig_chasm_dir}
        python2.7 RunChasm driver0 {output.snvbox_input} -g --pickone
        cd -
        """

rule splitMafChasm:
    input:
        join(benchmark_dir, 'snvbox_input/{benchmark_maf}.hg38.maf')
    params:
        chasm_dir=orig_chasm_dir,
        out_prefix=join(benchmark_dir, 'snvbox_input/{benchmark_maf}')
    output:
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver0.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver1.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver2.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver3.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver4.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver5.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver6.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver7.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver8.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_driver9.txt'),
        join(benchmark_dir, 'snvbox_input/{benchmark_maf,patrick_et_al|msk_impact|mc3}_split_passenger.txt'),
    shell:
        "python scripts/benchmark/split_maf.py "
        "   -m {input} "
        "   -c {params.chasm_dir}/drivers_gene_crossval/ "
        "   -o {params.out_prefix}"
