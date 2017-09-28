mydir=$1
n=30
for i in {0..29} ; do
    snakemake -s Snakefile chasm2_pretrained -p -j 999 -w 100 --max-jobs-per-second 1 --config  mutations="$mydir/sample$i/mutations.maf" output_dir="$mydir/sample$i" mutsigcv_dir="output/mutsigcv_formatted" trained_model="/projects/CHASM2/output/pancancer_v8_regular/cv_trained_model" trained_classifier="output/pancancer_v8_regular/2020plus.Rdata" --cluster-config cluster.yaml --cluster "qsub -cwd -pe smp {threads} -l mem_free={cluster.mem},h_vmem={cluster.vmem} -v PATH=$PATH -o $mydir/sample$i/sge_log -e $mydir/sample$i/sge_log"
    sleep 5
done
