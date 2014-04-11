#!/bin/bash
#By Samuel George; samuel@ras.ucalgary.ca
#Original: 27/02/2010
#Mod: 01/03/2010 - added support for more machines.
#Purpose: start lots of jobs. Script to be place in base dir of all files to be processed.
#just works on one machine at the moment.
#some ideas: http://pebblesinthesand.wordpress.com/2008/05/22/a-srcipt-for-running-processes-in-parallel-in-bash/
#	     ##http://ubuntuforums.org/showthread.php?t=382330
#run in correct dir, i.e. above the directory of sources like gal_parra.sh &
#uses passwordless ssh to make life easier, see passwordless_login.txt

date
dir_use=$PWD
array=(`ls -d band?/55*/`) #need to ensure no other directories.
echo ${array[*]}  
#echo ${array[*]}  >> "logfile"
len=${#array[*]}
maxjobs=8 #cant be higher as we are always overwriting the same file for each source / beam #has to equal number of sources note starts at zero
jobnumber=0
ddmmyy=(`date`)
#echo ${ddmmyy[*]} | mail -s "Started conversion jobs..." gu
#loop over the maximum number of jobs based on the number of files in array
while [ $jobnumber -lt $len ]; do
        jobsrunning=0
	while [ $jobsrunning -le $maxjobs -a $jobnumber -lt $len ]
                #start jobs up till maximum and then wait for them to finish before continuing.
		do
			if [ `ssh crow00 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow00..."
				ssh  crow00 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow01 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow01..."
				ssh  crow01 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow02 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow02..."
				ssh  crow02 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow03 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow03..."
				ssh  crow03 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow04 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow04..."
				ssh  crow04 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow05 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow05..."
				ssh  crow05 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow06 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow06..."
				ssh  crow06 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow07 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow07..."
				ssh  crow07 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			elif [ `ssh crow08 "ps -e | grep convert " | wc -l` == 0 ]
			then
				echo "Submitting ${array[$jobnumber]} to crow08..."
				ssh  crow08 "cd $dir_use/${array[$jobnumber]}; ./convert.sh &" &
	 			jobsrunning=$(( $jobsrunning + 1 ))
 				jobnumber=$(( $jobnumber + 1 ))
			fi
			ddmmyy=(`date`)
			#echo ${ddmmyy[*]} | mail -s "Processed $jobnumber  days so far..." gurm.sukhpreet@gmail.com
		done 
	wait
        echo $jobnumber
done

ddmmyy=(`date`)
#echo ${ddmmyy[*]} | mail -s "Processed all of the conversion..." $email
