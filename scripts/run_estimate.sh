function new_name {
        sed -e 's/^\([^_]*\).*\.csv$/\1_'"$1"'.csv/'
}

for F in $(cat $1) ; do
  echo $F
  ./estimate_abundance.sh -F $F -D clark/db --highconfidence -a 0.01 --krona > `echo $F | new_name clark_OTU`
  mv results.krn $F.krn
done
