#!/bin/bash

echo running $0 with PID $$ at `uname -n`
date
modes="Gau"
regud="true smo exp gau sph"
regu=(0 0 15 15 15)
rsto=(1 1 1 11 21)

var='1'

cur=`pwd`

skel=$cur/skel

for mode in $modes;do
	for x in *_$mode;do
		ix=`echo $x|tr '_' ' '|awk '{print $2}'`
		iy=`echo $x|tr '_' ' '|awk '{print $4}'`
		echo ix:$ix iy:$iy
		cd $x
		crm_skel=crm_skel
		if [ -d "$crm_skel" ];then
			echo $crm_skel exists
		else
			mkdir $crm_skel
			mv config $crm_skel
			mv exe $crm_skel
			mv grid $crm_skel
			mv mod $crm_skel
			mv rho $crm_skel
		fi
		pub=`pwd`/pub
		if [ -d $pub ];then
			echo $pub already there
		else
			mkdir $pub
		fi
		let i=0 # counts the different regus and rstos
		for z in $regud;do
			echo $z $i ${regu[i]} ${rsto[i]}
			if [ -z "$var" ];then
				exit
			if [ -d "$z" ];then
				rm -fR $z;
				rm -fR $pub/$z;
			fi
			mkdir $z
			mkdir $pub/$z
			cd $z
			cp -R $cur/$x/$crm_skel/* .
			cd exe
			date
			pwd
			cp $skel/exe/crtomo.cfg .
			mybin='CRTomo_'$x'_reg_'$z
			cp ~/bin/CRTomo_`uname -n` ./$mybin
			sed -in "s/INVDR/\.\.\/inv/g" crtomo.cfg
			sed -in "s/RSTO/${rsto[i]}/g" crtomo.cfg
			sed -in "s/REGU/${regu[i]}/g" crtomo.cfg
			sed -in "s/IX/$ix/g" crtomo.cfg
			sed -in "s/IY/$iy/g" crtomo.cfg
			sed -in "s/ANO/F/g" crtomo.cfg
			awk '{if(NR==12){print 0}else{print}}' crtomo.cfg > tmp.cfg
			cp tmp.cfg crtomo.cfg
			if [ "$z" == "true" ]; then
			  awk '{if(NR==8){printf("../rho/rho.dat\n")}else{print}}' crtomo.cfg > tmp.cfg
			  awk '{if(NR==15){print 0}else{print}}' tmp.cfg > crtomo.cfg
			fi
			if [ ${rsto[i]} -ge 40 ];then
				echo ''|./$mybin >& `echo $mybin|tr 'CRT' 'crt'`'.out'
			else
				./$mybin >& `echo $mybin|tr 'CRT' 'crt'`'.out'
			fi
			plt_all_crt $z
			else
			cd $z/exe
			fi
			cp inv.variogram* $pub/$z
			cp rho??_mag.pdf $pub/$z/model.pdf
			cp inv.stats_it $pub/$z
			cp inv_stat*.pdf $pub/$z
			let i=i+1
			cd $cur/$x
		done
		cd $pub
		variogram.sh $ix $iy $mode
		cd $cur
	done
done
