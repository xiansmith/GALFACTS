#!/bin/bash
#
# For each calibration run, we get files named this: 
# a2130_dtm.20091218.b2s0g0.00100.fits
#
# What we need is a directory structure SOURCE/BAND?/MJDDATE
# and the data sorted into those directories
# 
# For each run add the sources, and dates to the following variables
#
# Christian Smith

beams="0 1 2 3 4 5 6"

SOURCES="S0957+161 S1009+140 S1026+064 S1041+027 S1054+032 S1106-008 S1123+055 S1135-003"

DATESCAL=(20091218 20091219 20091220 20091221 20091222 20091223 20091224 20091226 20091227 20091228 20091229)

DATESMJD=(55183 55184 55185 55186 55187 55188 55189 55191 55192 55193 55194)

basedir=`pwd`
numdays=${#DATESMJD[*]}

echo BASE DIR is: $basedir
echo Number of Days = $numdays

for s in $SOURCES
  do
  cd $basedir
  cd $s
  
  for ((i=0; i < $numdays; i++ )) 
    do
    
#for j in ${DATESMJD[$i]}
      #do
  
  j=${DATESMJD[$i]}
  mkdir -p band0/$j
  mkdir -p band1/$j

  find . -name a2130_dtm.${DATESCAL[$i]}*s0*.fits -exec mv {} band0/$j \;
  find . -name a2130_dtm.${DATESCAL[$i]}*s1*.fits -exec mv {} band1/$j \;

  cd band0/$j
      
  # don't bother with empty dirs
  if [ `ls a2130* | wc -l` -gt "1" ];
      then
      for beam in $beams
	do
		# script to convert and run in parallel
	ls -1 *b${beam}s0*.fits > fitsfiles-$j-b${beam}s0.lis
	echo "#!/bin/bash" > convert.sh
	echo "beams=\"0 1 2 3 4 5 6\"" >> convert.sh
	echo "for beam in \$beams" >> convert.sh
	echo "do" >> convert.sh
	echo "/n/ras/processed/FIELD1/dtm2spec a2130 $j 0 \${beam} . fitsfiles-$j-b\${beam}s0.lis" >>  convert.sh
	echo "done" >> convert.sh
	chmod a+x convert.sh
      done
  fi
   
  cd ../../band1/$j
  
  # don't bother with empty dirs
  if [ `ls a2130* | wc -l` -gt "1" ];
      then
      for beam in $beams
	do
      	# script to convert and run in parallel
	ls -1 *b${beam}s1*.fits > fitsfiles-$j-b${beam}s1.lis
	echo "#!/bin/bash" > convert.sh
	echo "beams=\"0 1 2 3 4 5 6\"" >> convert.sh
	echo "for beam in \$beams" >> convert.sh
	echo "do" >> convert.sh
	echo "/n/ras/processed/FIELD1/dtm2spec a2130 $j 1 \${beam} . fitsfiles-$j-b\${beam}s1.lis" >>  convert.sh
	echo "done" >> convert.sh
	chmod a+x convert.sh
      done
  fi
  
  cd ../..
  done

done 

