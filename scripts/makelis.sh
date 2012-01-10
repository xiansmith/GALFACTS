#!/bin/bash
array=(`ls -d 5*/`)
echo ${array[*]}
len=${#array[*]}
count=0
while [ $count -lt $len ];
do

  beams="0 1 2 3 4 5 6"

  cd ${array[$count]}

  pwd 

  echo "*** band0 ***"

  for beam in $beams
  do

    echo "*** beam${beam} ***"

    ls -1 *b${beam}s0*.fits > fitsfiles_b${beam}s0.lis

    #/n/ras1/people/ricci/galfacts/bin/dtm2spec a2130 54785 0 ${beam} . fitsfiles_081114_b${beam}s0.lis

  done

  echo "*** band1 ***"

  for beam in $beams
  do

    echo "*** beam${beam} ***"

    ls -1 *b${beam}s1*.fits > fitsfiles_b${beam}s1.lis

    #/n/ras1/people/ricci/galfacts/bin/dtm2spec a2130 54785 1 ${beam} . fitsfiles_081114_b${beam}s1.lis

  done

  cd ..
  count=$(( $count + 1 ))
done
