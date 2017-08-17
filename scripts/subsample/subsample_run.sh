mydir=$1
n=30
for i in {0..29} ; do
    snakemake -s Snakefile chasm2_pretrained -p -j 999 -w 100 --max-jobs-per-second 1 --config  mutations="$mydir/sample$i/mutations.maf" output_dir="$mydir/sample$i" mutsigcv_dir="tmp2/mutsigcv_formatted" trained_model="/projects/CHASM2_results/output/snvbox_hg38_6_13_2017/cv_trained_no_gtex_no_exac_hotmaps_multip" trained_classifier="output/cancer_types/2020plus_10k.Rdata" --cluster-config cluster.yaml --cluster "qsub -cwd -pe smp {threads} -l mem_free={cluster.mem},h_vmem={cluster.vmem} -v PATH=$PATH -o $mydir/sample$i/sge_log -e $mydir/sample$i/sge_log"
    sleep 5
done
