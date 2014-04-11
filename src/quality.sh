#!/bin/bash
#change the 12 items below suitably
#executables directory
BINDIR=/export/cheetah2/people/sguram/wasim_data/galfacts/src/bin
#location of the data files
OBSDIR=/export/cheetah2/people/sguram/wasim_data/data
SMOOTH=0 #parameter to do hanning smoothing (1) or not (0)i
NUMSIGMA=5 #sigma value for rfi detection
NUMSIGMATHRESH=100 #sgima threshold for rfi detection
IGNOREA_LOW=0 #lower channel number for the spectral line to ignore during rfi detection
IGNOREA_HIGH=0 #higher channel number for the spectral line to ignore during rfi detection
PROJ="A2303" #filename prefix for the .spec files
MAX_CHANNELS=2048 #number of channels in the data
days="54773 54774 54778" #list of days
beams="0 1 2 3 4 5 6" #list of beams
RFI_XRANGE="[1365:1465]" #the band

PROG=$BINDIR/spec2fits #executable name for producing text tables from the .spec files
PLOT=.plot #temp file name

#set these to 1 or 0 to turn on and off functionality
PROCESS_SPEC=1
POINTING_PLOTS=1
AGGREGATE_PLOTS=1
BANDAVG_PLOTS=1
CALAVG_PLOTS=1
TIMEAVG_PLOTS=1
RFI_PLOTS=1
HTML=1

#create the local directories
for dir in $days 
do
	mkdir -p $dir
done


for beam in $beams
do
echo "*** beam${beam} ***"


#run the processing code
if [ "$PROCESS_SPEC" = 1 ]; then
echo "Processing files"
for dir in $days 
do
	$PROG $OBSDIR/$dir/$PROJ*.beam${beam}*.$dir.spec $SMOOTH $NUMSIGMA $NUMSIGMATHRESH $IGNOREA_LOW $IGNOREA_HIGH 
done
fi

#determine plotting ranges from the data 
RA_MIN=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 1 | head -1 | awk '{print $1}'`
RA_MAX=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 1 | tail -1 | awk '{print $1}'`
RA_RANGE="[$RA_MIN:$RA_MAX]"
echo RA Range: $RA_RANGE

DEC_MIN=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 2 | head -1 | awk '{print $2}'`
DEC_MAX=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 2 | tail -1 | awk '{print $2}'`
DEC_RANGE="[$DEC_MIN:$DEC_MAX]"
echo DEC Range: $DEC_RANGE

BANDAVG_YRANGE="[0.6:2.0]"
TIMEAVG_YRANGE="[0:2.0]"
CALAVG_YRANGE="[-0.1:0.1]"

#create the pointing plots
if [ "$POINTING_PLOTS" = 1 ]; then
echo "Creating pointing plots"
for dir in $days
do
	echo "set term post color" > $PLOT
	echo "set output '${dir}/beam${beam}_pointing.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} Pointing Positions'" >> $PLOT
	echo "set xlabel 'RA (hours)'" >> $PLOT
	echo "set xrange $RA_RANGE" >> $PLOT
	echo "set ylabel 'DEC (degrees)'" >> $PLOT
	echo "set yrange $DEC_RANGE" >> $PLOT
	echo "plot '< cat *.beam${beam}.$dir.spec_pointing.dat' using 1:2 with lines notitle" >> $PLOT
	gnuplot $PLOT
	convert "${dir}/beam${beam}_pointing.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_pointing_small.png"
	convert "${dir}/beam${beam}_pointing.ps" -rotate 90 "${dir}/beam${beam}_pointing_large.png"
done
fi

if [ "$AGGREGATE_PLOTS" = 1 ]; then
echo "Creating aggregate plot"

#Create the aggregate pointing plot
echo "set term post color" > $PLOT
echo "set output 'beam${beam}_pointing.ps'" >> $PLOT
echo "set title 'Beam${beam} Aggregate Pointing Positions'" >> $PLOT
echo "set xlabel 'RA (hours)'" >> $PLOT
echo "set xrange $RA_RANGE" >> $PLOT
echo "set ylabel 'DEC (degrees)'" >> $PLOT
echo "set yrange $DEC_RANGE" >> $PLOT
echo "plot '< cat *.beam${beam}.*.spec_pointing.dat' using 1:2 with lines notitle" >> $PLOT
gnuplot $PLOT
convert "beam${beam}_pointing.ps" -rotate 90 "beam${beam}_pointing_large.png"
fi


#create the band average plots
if [ "$BANDAVG_PLOTS" = 1 ]; then
echo "Creating bandavg plots"
for dir in $days 
do
	echo "set term post color" > $PLOT
	echo "set output '${dir}/beam${beam}_bandavg.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} Band Average'" >> $PLOT
	echo "set xlabel 'RA (Hours)'" >> $PLOT
	echo "set xrange $RA_RANGE" >> $PLOT
	echo "set ylabel 'Power'" >> $PLOT
	echo "set yrange $BANDAVG_YRANGE" >> $PLOT
	echo "plot \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_bandavg.dat' using 1:4 with dots title 'CalON XX', \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_bandavg.dat' using 1:8 with dots title 'CalOFF XX', \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_bandavg.dat' using 1:7 with dots title 'CalON YY', \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_bandavg.dat' using 1:11 with dots title 'CalOFF YY'" >> $PLOT
	gnuplot $PLOT
	convert "${dir}/beam${beam}_bandavg.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_bandavg_small.png"
	convert "${dir}/beam${beam}_bandavg.ps" -rotate 90 "${dir}/beam${beam}_bandavg_large.png"
done
fi

#create the time average plots
if [ "$TIMEAVG_PLOTS" = 1 ]; then
echo "Creating timeavg plots"
for dir in $days 
do
	echo "set term post color" > $PLOT
	echo "set output '${dir}/beam${beam}_timeavg.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} beam${beam}_Time Average'" >> $PLOT
	echo "set xrange [0:$MAX_CHANNELS]" >> $PLOT
	echo "set xlabel 'Channel" >> $PLOT
	echo "set ylabel 'Power'" >> $PLOT
	echo "set yrange $TIMEAVG_YRANGE" >> $PLOT
	echo "plot \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_timeavg.dat' using 1 with dots title 'CalON XX', \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_timeavg.dat' using 5 with dots title 'CalOFF XX', \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_timeavg.dat' using 4 with dots title 'CalON YY', \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_timeavg.dat' using 8 with dots title 'CalOFF YY'" >> $PLOT
	gnuplot $PLOT
	convert "${dir}/beam${beam}_timeavg.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_timeavg_small.png"
	convert "${dir}/beam${beam}_timeavg.ps" -rotate 90 "${dir}/beam${beam}_timeavg_large.png"
done
fi

#create the CAL  plots
if [ "$CALAVG_PLOTS" = 1 ]; then
echo "Creating cal plots"
for dir in $days 
do
	echo "set term post color" > $PLOT
	echo "set output '${dir}/beam${beam}_calavg.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} Cal Average'" >> $PLOT
	echo "set xlabel 'RA (Hours)'" >> $PLOT
	echo "set xrange $RA_RANGE" >> $PLOT
	echo "set ylabel 'Power'" >> $PLOT
	echo "set yrange $CALAVG_YRANGE" >> $PLOT
	echo "plot \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_bandavg.dat' using 1:(\$4-\$8) with dots title 'Cal XX', \\" >> $PLOT
	echo "'< cat *.beam${beam}.$dir.spec_bandavg.dat' using 1:(\$7-\$11) with dots title 'Cal YY'" >> $PLOT
	gnuplot $PLOT
	convert "${dir}/beam${beam}_calavg.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_calavg_small.png"
	convert "${dir}/beam${beam}_calavg.ps" -rotate 90 "${dir}/beam${beam}_calavg_large.png"
done
fi



#create the RFI plots
if [ "$RFI_PLOTS" = 1 ]; then
echo "Creating rfi plots"
for dir in $days 
do
	echo "set term post color" > $PLOT
	echo "set output '${dir}/beam${beam}_rfi.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} RFI'" >> $PLOT
	echo "set xlabel 'Frequency (MHz)'" >> $PLOT
	echo "set xrange $RFI_XRANGE" >> $PLOT
	echo "set ylabel 'Time (seconds)'" >> $PLOT
	echo "plot '< cat *.beam${beam}.$dir.spec_rfi.dat' using 2:5 with dots notitle" >> $PLOT
	gnuplot $PLOT
	convert "${dir}/beam${beam}_rfi.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_rfi_small.png"
	convert "${dir}/beam${beam}_rfi.ps" -rotate 90 "${dir}/beam${beam}_rfi_large.png"
done
fi

if [ "$HTML" = 1 ]; then

#create the per beam html files
echo "Creating per beam html"
html_file="beam${beam}.html"
echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">" > $html_file
echo "<html>" >> $html_file
echo "<head>" >> $html_file
echo "<title>${PROJ} WAPPS beam${beam} Data Quality Heuristics</title>" >> $html_file
echo "<meta name=author content=\"Jeff Dever\">" >> $html_file
echo "</head>" >> $html_file
echo "<body>" >> $html_file
echo "<h1 align=center>${PROJ} beam${beam} WAPPS Data Quality Heuristics</h1>" >> $html_file

#link navigation bar
echo "<table align=center><tr><td align=center>" >> $html_file
echo "<a href=index.html>main</a>" >> $html_file
echo "</td></tr><tr><td align=center>" >> $html_file
for beamid in $beams; do
echo "<a href=beam${beamid}.html>beam${beamid}</a> " >> $html_file
done
echo "</td></tr><tr><td align=center>" >> $html_file
for dir in $days; do
echo "<a href=${dir}.html>${dir}</a> " >> $html_file
done
echo "</td></tr></table>" >> $html_file


echo "<table border=1>" >> $html_file
echo "<tr>" >> $html_file
echo "<th></th>" >> $html_file
echo "<th>Pointing Plot</th>" >> $html_file
echo "<th>Band Average</th>" >> $html_file
echo "<th>Cal Average</th>" >> $html_file
echo "<th>Time Average</th>" >> $html_file
echo "<th>RFI Detection</th>" >> $html_file
echo "</tr>" >> $html_file

for dir in $days 
do
echo "<tr>" >> $html_file

echo "<th align=center>" >> $html_file
echo "$dir<br>beam${beam}" >> $html_file
echo "</th>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_pointing_large.png><img src=\"${dir}/beam${beam}_pointing_small.png\" alt=\"Pointing Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_pointing.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_bandavg_large.png><img src=\"${dir}/beam${beam}_bandavg_small.png\" alt=\"Band Average Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_bandavg.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_calavg_large.png><img src=\"${dir}/beam${beam}_calavg_small.png\" alt=\"CAL Average Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_calavg.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_timeavg_large.png><img src=\"${dir}/beam${beam}_timeavg_small.png\" alt=\"Time Average Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_timeavg.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_rfi_large.png><img src=\"${dir}/beam${beam}_rfi_small.png\" alt=\"RFI Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_rfi.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file
echo "</tr>" >> $html_file

done
echo "</table>" >> $html_file

echo "<hr>" >> $html_file

echo "<table width=\"100%\">" >> $html_file
echo "<tr>" >> $html_file
date=`date`
echo "<td align=left width=\"20%\">Jeff Dever<br>Centre for Radio Astronomy<br>University of Calgary</td>" >> $html_file
echo "<td align=center>Generated: $date</td>" >> $html_file
echo "<td align=right width=\"20%\"><a href=\"http://validator.w3.org/check?uri=referer\"><img src=\"http://www.w3.org/Icons/valid-html401\" alt=\"Valid HTML 4.01 Transitional\" height=31 width=88></a></td>" >> $html_file
echo "</tr>" >> $html_file
echo "</table>" >> $html_file

echo "</body>" >> $html_file
echo "</html>" >> $html_file

fi

done

if [ "$HTML" = 1 ]; then

#create the per day html files
echo "Creating per day html"
for day in $days 
do
html_file="${day}.html"
echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">" > $html_file
echo "<html>" >> $html_file
echo "<head>" >> $html_file
echo "<title>${PROJ} WAPPS ${day} Data Quality Heuristics</title>" >> $html_file
echo "<meta name=author content=\"Jeff Dever\">" >> $html_file
echo "</head>" >> $html_file
echo "<body>" >> $html_file
echo "<h1 align=center>${PROJ} ${day} WAPPS Data Quality Heuristics</h1>" >> $html_file

#link navigation bar
#TODO: use conditional to drop the link for the current page
echo "<table align=center><tr><td align=center>" >> $html_file
echo "<a href=index.html>main</a>" >> $html_file
echo "</td></tr><tr><td align=center>" >> $html_file
for beamid in $beams; do
echo "<a href=beam${beamid}.html>beam${beamid}</a> " >> $html_file
done
echo "</td></tr><tr><td align=center>" >> $html_file
for dir in $days; do
echo "<a href=${dir}.html>${dir}</a> " >> $html_file
done
echo "</td></tr></table>" >> $html_file

echo "<h2 align=center>Pointing Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_pointing_large.png><img src=\"${day}/beam3_pointing_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_pointing_large.png><img src=\"${day}/beam4_pointing_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_pointing_large.png><img src=\"${day}/beam2_pointing_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_pointing_large.png><img src=\"${day}/beam0_pointing_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_pointing_large.png><img src=\"${day}/beam5_pointing_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_pointing_large.png><img src=\"${day}/beam1_pointing_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_pointing_large.png><img src=\"${day}/beam6_pointing_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file


echo "<h2 align=center>Band Average Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_bandavg_large.png><img src=\"${day}/beam3_bandavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_bandavg_large.png><img src=\"${day}/beam4_bandavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_bandavg_large.png><img src=\"${day}/beam2_bandavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_bandavg_large.png><img src=\"${day}/beam0_bandavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_bandavg_large.png><img src=\"${day}/beam5_bandavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_bandavg_large.png><img src=\"${day}/beam1_bandavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_bandavg_large.png><img src=\"${day}/beam6_bandavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file

echo "<h2 align=center>CAL Average Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_calavg_large.png><img src=\"${day}/beam3_calavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_calavg_large.png><img src=\"${day}/beam4_calavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_calavg_large.png><img src=\"${day}/beam2_calavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_calavg_large.png><img src=\"${day}/beam0_calavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_calavg_large.png><img src=\"${day}/beam5_calavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_calavg_large.png><img src=\"${day}/beam1_calavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_calavg_large.png><img src=\"${day}/beam6_calavg_small.png\" alt=\"Band Average Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file

echo "<h2 align=center>Time Average Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_timeavg_large.png><img src=\"${day}/beam3_timeavg_small.png\" alt=\"Time Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_timeavg_large.png><img src=\"${day}/beam4_timeavg_small.png\" alt=\"Time Average Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_timeavg_large.png><img src=\"${day}/beam2_timeavg_small.png\" alt=\"Time Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_timeavg_large.png><img src=\"${day}/beam0_timeavg_small.png\" alt=\"Time Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_timeavg_large.png><img src=\"${day}/beam5_timeavg_small.png\" alt=\"Time Average Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_timeavg_large.png><img src=\"${day}/beam1_timeavg_small.png\" alt=\"Time Average Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_timeavg_large.png><img src=\"${day}/beam6_timeavg_small.png\" alt=\"Time Average Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file


echo "<h2 align=center>RFI Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_rfi_large.png><img src=\"${day}/beam3_rfi_small.png\" alt=\"RFI Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_rfi_large.png><img src=\"${day}/beam4_rfi_small.png\" alt=\"RFI Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_rfi_large.png><img src=\"${day}/beam2_rfi_small.png\" alt=\"RFI Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_rfi_large.png><img src=\"${day}/beam0_rfi_small.png\" alt=\"RFI Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_rfi_large.png><img src=\"${day}/beam5_rfi_small.png\" alt=\"RFI Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_rfi_large.png><img src=\"${day}/beam1_rfi_small.png\" alt=\"RFI Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_rfi_large.png><img src=\"${day}/beam6_rfi_small.png\" alt=\"RFI Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file


echo "<hr>" >> $html_file

echo "<table width=\"100%\">" >> $html_file
echo "<tr>" >> $html_file
date=`date`
echo "<td align=left width=\"20%\">Jeff Dever<br>University of Calgary<br>Radio Astronomy Lab<br></td>" >> $html_file
echo "<td align=center>Generated: $date</td>" >> $html_file
echo "<td align=right width=\"20%\"><a href=\"http://validator.w3.org/check?uri=referer\"><img src=\"http://www.w3.org/Icons/valid-html401\" alt=\"Valid HTML 4.01 Transitional\" height=31 width=88></a></td>" >> $html_file
echo "</tr>" >> $html_file
echo "</table>" >> $html_file

echo "</body>" >> $html_file
echo "</html>" >> $html_file

done


#generate index.html
echo "Creating index html"

html_file="index.html"
echo "<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">" > $html_file
echo "<html>" >> $html_file
echo "<head>" >> $html_file
echo "<title>${PROJ} Data Quality Heuristics</title>" >> $html_file
echo "<meta name=author content="Jeff Dever">" >> $html_file
echo "</head>" >> $html_file
echo "<body>" >> $html_file
echo "<h1 align=center>${PROJ} Data Quality Heuristics</h1>" >> $html_file
echo "<table align=center>" >> $html_file
echo "<tr><td align=center>" >> $html_file
for beam in $beams
do
echo "<a href=beam${beam}.html>beam${beam}</a> " >> $html_file
done
echo "</td></tr><tr><td align=center>" >> $html_file
for dir in $days
do
echo "<a href=$dir.html>$dir</a> " >> $html_file
done
echo "</td></tr></table>" >> $html_file
echo "<br>" >> $html_file
echo "<center><a href=beam0_pointing.ps><img src=beam0_pointing_large.png alt="Pointing Plot"></a></center>" >> $html_file
echo "<br>" >> $html_file
echo "<hr>" >> $html_file
echo "<table width="100%">" >> $html_file
echo "<tr>" >> $html_file
echo "<td align=left width=20%>Jeff Dever<br>University of Calgary<br>Radio Astronomy Lab</td>" >> $html_file
echo "<td align=center></td>" >> $html_file
echo "<td align=right width=20%><a href="http://validator.w3.org/check?uri=referer"><img src="http://www.w3.org/Icons/valid-html401" alt="Valid HTML 4.01 Transitional" height=31 width=88></a></td>" >> $html_file
echo "</tr>" >> $html_file
echo "</table>" >> $html_file
echo "</body>" >> $html_file
echo "</html>" >> $html_file

fi

#rm .plot

exit 0
