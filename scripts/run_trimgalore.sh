for F in $(cat index.txt) ; do
        trim_galore --paired --nextera --length 70 --fastqc --fastqc_args "-t 20 -o qc/ny/fastqc.trimgalore.afterqc/" -o fastq/ny/trimgalore.afterqc $F\_1.good.fq $F\_2.good.fq
done
