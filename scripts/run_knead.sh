for F in $(cat index.txt); do
    kneaddata --input $F\_1.good_val_1.fq --input $F\_2.good_val_2.fq -db HS.removed/db/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens --output HS.removed/ny/ --bypass-trim -t 20
done
