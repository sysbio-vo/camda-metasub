function new_name {
	sed -e 's/^\([^_]*\).*\.csv$/\1_'"$1"'.csv/'
}

for file in *.csv; do
	mv $file `echo $file | new_name $1`
done
