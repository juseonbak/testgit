#!/bin/bash

type=ML
# current directory
pwd=`pwd` 

# Whehter to use qsub or not (1: use qsub disabled 0: command line, else: generate L1L2_fnames.inp)
use_qsub=0
use_omicldo2=1
ver=3

# Set up output and input directory (RSL)
lv2dir=/home/JDB/data/paper2/sonpxl/ 
lv2dir=/home/JDB/data/paper6/omison/sonpxl2/
#lv2dir=/JDB/data/gems/sonpxl/

listdir=/home/JDB/Data/SONDE/collocations/
listdir=/home/JDB/Data/SONDE/omison/


if [ ! -e $lv2dir ]
then
    mkdir $lv2dir
fi

lv2prefix=OMIO3PROF_$type'-'
solprefix=OMI-Aura_L1-OML1BIRR_
radprefix=OMI-Aura_L1-OML1BRUG_
soldir=/omi/live/col3/OML1BIRR/
cldprefix=OMI-Aura_L2-OMCLDO2_

nx=0  # number of pixels xtrack      (+/-)
ny=0  # number of pixels along track (+/-)

grplist=`cat /home/JDB/Data/SONDE/grplist.txt`
grplist=`cat ./grplist.txt`

for grp in $grplist
do
    listfname=$grp'_exact_colist.txt'
    listfname=$grp'_UV1_colist_2004-2008.txt'
   
    sta=${listfname:0:${#grp}}      
    subdir=$lv2dir$sta/
    if [ ! -e $subdir ]
    then
        mkdir $subdir
   fi
        
    cat $listdir$listfname | while read aline
    do
        #echo $aline
        aline=($aline)        # turn line into an array
        date=${aline[0]}
        orb=${aline[1]}
        sline1=${aline[2]}
        pix=${aline[3]}
       # sta=${aline[4]}
       # utc=${aline[5]}
       # dis=${aline[6]}
        date=${date:0:9}
        year=${date:0:4}
        mon=${date:5:2}
        day=${date:7:2}

        raddir=/home/Data/OMI/1_OML1BRUG/$year/
        clddir=/home/Data/OMI/2_OML2CLDO2/$year/
        
     
#	subdir=$lv2dir$sta'/'$type'/'
	if [ ! -e $subdir ]
	then
	     mkdir $subdir
  fi
        
        # set up range along the track
       # sline1=`expr $sline1 + 1`  # input line number starts from 1 in my algorithm (04/19/2008) when use xiong's collocation
        sline2=`expr $sline1 + $ny`
        sline1=`expr $sline1 - $ny`
        
        # set up range across the track
        if [ `expr $pix % 2` -eq 0 ]
        then
            pix=`expr $pix - 1`
        fi

        let dx=$nx*2
        pix1=`expr $pix - $dx`
        pix2=`expr $pix + $dx + 1`
        if [[ "$pix1" -lt 1 ]]
        then
            pix1=1
        fi
        if [[ "$pix2" -gt 60 ]]
        then
            pix2=60
        fi
        #echo $pix1 $pix2
        
        pixsel="$sline1 $sline2 $pix1 $pix2"
        timesel=$date't????-o*'$orb
        #echo $timesel $pixsel
        
        # initialized the orbits to be arrays
        # note the index starts from zero
        orbs[0]=''
        norb=0
        
        for radfname in `ls $raddir$radprefix$timesel*he4`
        do
            orb=${radfname#*-o}  
            orb=${orb:0:5}
            #echo $orb
            mondayutc=${radfname#*$radprefix}
            mondayutc=${mondayutc:0:14}
            monday=${mondayutc:0:9}
            #echo $mondayutc
            
            #solfname=`ls $soldir$year/$mon/$day/$solprefix$monday*he4`
            cldfname=`ls $clddir$cldprefix*o$orb*he5`
            lv2fname=$subdir$lv2prefix'o'$orb
            # Solar radiance could not be found, make up a name in case use backup radiance
            #if [ -z $solfname ]
            #then 
                #continue
                solfname=$soldir$solprefix$mondayutc'-o'$orb'_v00'$ver'.he4'
            #fi
            
            # fake a cloud file name (not used anyway for the moment)
            if [ -z $cldfname ]
                then
                cldfname=$clddir$cldprefix$mondayutc'-o'$orb'_v00'$ver'.he5'
            fi
            
            orbs[$norb]=$orb
            norb=`expr $norb + 1`
                    
            # Write input file INP/L1L2_fnames.inp
            # use echo -e to tell that there are control variable new line (\x0a )
            echo -e "$solfname \x0a$radfname \x0a$cldfname \x0a$lv2fname \x0a$pixsel" > INP/L1L2_fnames.inp
            if [ $use_qsub -eq 0 ] 
                then  
                echo $radfname
                echo $lv2fname          
                SOMIPROF'.exe' 
            #    sleep 30
            fi
        done  
    done
done

