for F in $(cat index.txt) ; do
    dsrc d /bi/ala/copy/2017-CAMDA/MetaSUB/Sacramento/data/$F.fastq.dsrc $F.fastq
    cat $F | paste - - - - | sort | tr '\t' '\n' > $F.fastq.out
    cat $F.fastq.out | paste - - - - | grep ' 1:' | tr '\t' '\n' > $F\_R1.fastq
    cat $F.fastq.out | paste - - - - | grep ' 2:' | tr '\t' '\n' > $F\_R2.fastq
    rm $F.fastq $F.fastq.out
done

