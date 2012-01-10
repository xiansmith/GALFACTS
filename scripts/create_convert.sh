#!/bin/bash
array=(`ls -d 5*`)
echo ${array[*]}
len=${#array[*]}
count=0 
while [ $count -lt $len ];
do
	echo "#!/bin/bash" >> ${array[$count]}/convert.sh
	echo "beams=\"0 1 2 3 4 5 6\"" >> ${array[$count]}/convert.sh
	echo "for beam in \$beams" >> ${array[$count]}/convert.sh
	echo "do" >> ${array[$count]}/convert.sh
	day=(`ls -1 ${array[$count]}/*b0s0*.lis | cut -d '_' -f 2`) 
	echo "/n/ras/data/FIELD2/dtm/bin/dtm2spec a2130 ${array[$count]} 0 \${beam} . fitsfiles_b\${beam}s0.lis" >>  ${array[$count]}/convert.sh
	echo "done" >> ${array[$count]}/convert.sh
	count=$(( $count + 1 ))	
done
