#!/bin/bash

tools="plotCRTmod_batch.m CRTomo epstopdf matlab mogrify";

echo "checking for tools"
forcecheck(){
	Ret_val=`which $1`
	if [ -z "$Ret_val" ];then
		echo "$1 not in search path.. searching $HOME..."
		Ret_val=`find $HOME -name "$1" -print`
		if [ -z "$arg" ];then
			exit
		fi
	fi
	return
}

i=1;
for x in $tools;do
	forcecheck $x
	case "$i" in
		1)
			plotmod=$Ret_val
		;;
		2)
			crt=$Ret_val
		;;
		3)
			epsp=$Ret_val
		;;
		4)
			mtlb=$Ret_val
		;;
		5)
			mogri=$Ret_val
		;;
		*)
			echo "not possible"
			exit
		;;
	esac
	i=$[$i+1];
done

awk '{if (NR==2){print}}' crtomo.cfg > tmp.meshname
awk '{if (NR==3){print}}' crtomo.cfg > tmp.elecname
echo $1 > tmp.fenster

c1=`pwd|awk '/2_decades/{printf("%.1f\t%.1f",1.0,3.0)}'` 
c2=`pwd|awk '/3_decades/{printf("%.1f\t%.1f",0.5,3.5)}'` 
c3=`pwd|awk '/4_decades/{printf("%.1f\t%.1f",0.0,4.0)}'` 
echo range check: $c1 $c2 $c3 $2
echo $c1 $c2 $c3 $2 > tmp.range

wdir=`pwd|xargs basename`

if [ -n $1 ];then
    wdir=$1
fi

cp $crt crtomo.$wdir

./crtomo.$wdir >& $wdir'.crtrun'

$mtlb < $plotmod >& /dev/null

tmp=`cat tmp.lastmod`
myeps=`basename $tmp`'.eps'
mypdf=`basename $myeps '.eps'`'.pdf'
echo "making image conversion..."
$epsp $myeps -o=$1$mypdf >& /dev/null
$mogri -format jpg $myeps
echo "done $1"

#pdfcr=`which pdfcrop`
#forcecheck $pdfcr
#$pdfcr --margins 1 --hires $mypdf cropped.pdf >& /dev/null
