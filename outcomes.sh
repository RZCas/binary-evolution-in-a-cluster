shopt -s extglob

root_folder="storage/m1=m2=10"
for folder in $root_folder/*/
do
	outcomes_file="$root_folder/outcomes-$(basename $folder).txt"
	# touch $outcomes_file
	for file in $folder?(?|??|???|????).txt
	do
		echo "${file##*/}" >> "$outcomes_file"
		head -n3 $file | tail -n1 >> "$outcomes_file"
		tail -n2 $file >> "$outcomes_file"
		echo "" >> "$outcomes_file"
	done
done