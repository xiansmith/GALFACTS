#!/bin/bash
#change the 12 items below suitably
#executables directory
#BINDIR=/export/cheetah2/people/sguram/wasim_data/galfacts/src/bin
#location of the data files
OBSDIR=/n/fox/processed/cal
days="54785 54786 54787 54788 54789 54790 54793" #list of days
beams="0 1 2 3 4 5 6" #list of beams
#beams="0" #list of beams

PLOT=.plot #temp file name

#set these to 1 or 0 to turn on and off functionality
PROCESS_SPEC=0
FITPLOTS=1
POINTING_PLOTS=0
AGGREGATE_PLOTS=0
BANDAVG_PLOTS=0
CALAVG_PLOTS=0
TIMEAVG_PLOTS=0
RFI_PLOTS=0
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
DP_MIN=`grep -v '^#' $OBSDIR/$dir/beam$beam/rawbandcal.dat | sort -n -k 1 | head -1 | awk '{print $1}'`
DP_MAX=`grep -v '^#' $OBSDIR/$dir/beam$beam/rawbandcal.dat | sort -n -k 1 | tail -1 | awk '{print $1}'`
DP_RANGE="[$DP_MIN:$DP_MAX]"
echo Datapoint Range: $DP_RANGE

#determine plotting ranges from the data 
#RA_MIN=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 1 | head -1 | awk '{print $1}'`
#RA_MAX=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 1 | tail -1 | awk '{print $1}'`
#RA_RANGE="[$RA_MIN:$RA_MAX]"
#echo RA Range: $RA_RANGE

#DEC_MIN=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 2 | head -1 | awk '{print $2}'`
#DEC_MAX=`grep -v '^#' *beam${beam}*.$dir.spec_pointing.dat | sort -n -k 2 | tail -1 | awk '{print $2}'`
#DEC_RANGE="[$DEC_MIN:$DEC_MAX]"
#echo DEC Range: $DEC_RANGE

BANDAVG_YRANGE="[0.6:2.0]"
TIMEAVG_YRANGE="[0:2.0]"
CALAVG_YRANGE="[-0.1:0.1]"

#create the calfit plots
if [ "$FITPLOTS" = 1 ]; then
echo "Creating calfit plots"
for dir in $days
do
	echo "set term post color" > $PLOT
	#echo "set output '${dir}/beam${beam}_calfit.ps'" >> $PLOT
	echo "set output '${dir}/beam${beam}_calxxfit.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} Cal XX fits'" >> $PLOT
	echo "set xlabel 'n'" >> $PLOT
	echo "set xrange $DP_RANGE" >> $PLOT
	echo "set ylabel 'cal'" >> $PLOT
	echo "plot '< cat $OBSDIR/$dir/beam$beam/rawbandcal.dat' using 1:2 with dots notitle,'<cat $OBSDIR/$dir/beam$beam/smoothcal.dat' using 1:2 with lines notitle" >> $PLOT
	echo "set output '${dir}/beam${beam}_calyyfit.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} Cal YY fits'" >> $PLOT
	echo "set xlabel 'n'" >> $PLOT
	echo "set xrange $DP_RANGE" >> $PLOT
	echo "set ylabel 'cal'" >> $PLOT
	echo "plot '< cat $OBSDIR/$dir/beam$beam/rawbandcal.dat' using 1:3 with dots notitle ,'<cat $OBSDIR/$dir/beam$beam/smoothcal.dat' using 1:3 with lines notitle" >> $PLOT
	echo "set output '${dir}/beam${beam}_calxyfit.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} Cal XY fits'" >> $PLOT
	echo "set xlabel 'n'" >> $PLOT
	echo "set xrange $DP_RANGE" >> $PLOT
	echo "set ylabel 'cal'" >> $PLOT
	echo "plot '< cat $OBSDIR/$dir/beam$beam/rawbandcal.dat' using 1:4 with dots notitle ,'<cat $OBSDIR/$dir/beam$beam/smoothcal.dat' using 1:4 with lines notitle" >> $PLOT
	echo "set output '${dir}/beam${beam}_calyxfit.ps'" >> $PLOT
	echo "set title '${dir} beam${beam} Cal YX fits'" >> $PLOT
	echo "set xlabel 'n'" >> $PLOT
	echo "set xrange $DP_RANGE" >> $PLOT
	echo "set ylabel 'cal'" >> $PLOT
	echo "plot '< cat $OBSDIR/$dir/beam$beam/rawbandcal.dat' using 1:5 with dots notitle ,'<cat $OBSDIR/$dir/beam$beam/smoothcal.dat' using 1:5 with lines notitle" >> $PLOT
	gnuplot $PLOT
	convert "${dir}/beam${beam}_calxxfit.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_calxx_small.png"
	convert "${dir}/beam${beam}_calxxfit.ps" -rotate 90 "${dir}/beam${beam}_calxx_large.png"
	convert "${dir}/beam${beam}_calyyfit.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_calyy_small.png"
	convert "${dir}/beam${beam}_calyyfit.ps" -rotate 90 "${dir}/beam${beam}_calyy_large.png"
	convert "${dir}/beam${beam}_calxyfit.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_calxy_small.png"
	convert "${dir}/beam${beam}_calxyfit.ps" -rotate 90 "${dir}/beam${beam}_calxy_large.png"
	convert "${dir}/beam${beam}_calyxfit.ps" -rotate 90 -thumbnail 300x210 "${dir}/beam${beam}_calyx_small.png"
	convert "${dir}/beam${beam}_calyxfit.ps" -rotate 90 "${dir}/beam${beam}_calyx_large.png"
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
echo "<th>Cal XX Plot</th>" >> $html_file
echo "<th>Cal YY Plot</th>" >> $html_file
echo "<th>Cal XY Plot</th>" >> $html_file
echo "<th>Cal YX Plot</th>" >> $html_file
echo "</tr>" >> $html_file

for dir in $days 
do
echo "<tr>" >> $html_file

echo "<th align=center>" >> $html_file
echo "$dir<br>beam${beam}" >> $html_file
echo "</th>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_calxx_large.png><img src=\"${dir}/beam${beam}_calxx_small.png\" alt=\"XX Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_calxxfit.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_calyy_large.png><img src=\"${dir}/beam${beam}_calyy_small.png\" alt=\"YY Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_calyyfit.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_calxy_large.png><img src=\"${dir}/beam${beam}_calxy_small.png\" alt=\"XY Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_calxyfit.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

echo "<td align=center>" >> $html_file
echo "<a href=${dir}/beam${beam}_calyx_large.png><img src=\"${dir}/beam${beam}_calyx_small.png\" alt=\"YX Plot\"></a><br>" >> $html_file
echo "<a href=${dir}/beam${beam}_calyxfit.ps>Postscript Plot</a><br>" >> $html_file
echo "</td>" >> $html_file

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

echo "<h2 align=center>Cal XX Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_calxx_large.png><img src=\"${day}/beam3_calxx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_calxx_large.png><img src=\"${day}/beam4_calxx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_calxx_large.png><img src=\"${day}/beam2_calxx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_calxx_large.png><img src=\"${day}/beam0_calxx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_calxx_large.png><img src=\"${day}/beam5_calxx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_calxx_large.png><img src=\"${day}/beam1_calxx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_calxx_large.png><img src=\"${day}/beam6_calxx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file

echo "<h2 align=center>Cal YY Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_calyy_large.png><img src=\"${day}/beam3_calyy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_calyy_large.png><img src=\"${day}/beam4_calyy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_calyy_large.png><img src=\"${day}/beam2_calyy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_calyy_large.png><img src=\"${day}/beam0_calyy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_calyy_large.png><img src=\"${day}/beam5_calyy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_calyy_large.png><img src=\"${day}/beam1_calyy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_calyy_large.png><img src=\"${day}/beam6_calyy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file

echo "<h2 align=center>Cal XY Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_calxy_large.png><img src=\"${day}/beam3_calxy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_calxy_large.png><img src=\"${day}/beam4_calxy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_calxy_large.png><img src=\"${day}/beam2_calxy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_calxy_large.png><img src=\"${day}/beam0_calxy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_calxy_large.png><img src=\"${day}/beam5_calxy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_calxy_large.png><img src=\"${day}/beam1_calxy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_calxy_large.png><img src=\"${day}/beam6_calxy_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "</table>" >> $html_file

echo "<h2 align=center>Cal YX Plots</h2>" >> $html_file
echo "<table align=center border=1>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam3_calyx_large.png><img src=\"${day}/beam3_calyx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam4_calyx_large.png><img src=\"${day}/beam4_calyx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr>" >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam2_calyx_large.png><img src=\"${day}/beam2_calyx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam0_calyx_large.png><img src=\"${day}/beam0_calyx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam5_calyx_large.png><img src=\"${day}/beam5_calyx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "</td> </tr> " >> $html_file
echo "<tr> <td align=center>" >> $html_file
echo "<a href=${day}/beam1_calyx_large.png><img src=\"${day}/beam1_calyx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
echo "<a href=${day}/beam6_calyx_large.png><img src=\"${day}/beam6_calyx_small.png\" alt=\"Pointing Plot\"></a>" >> $html_file
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
echo "<title>${PROJ} Cal fit quality</title>" >> $html_file
echo "<meta name=author content="Sukhpreet Guram">" >> $html_file
echo "</head>" >> $html_file
echo "<body>" >> $html_file
echo "<h1 align=center>${PROJ} Cal fit Quality</h1>" >> $html_file
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
