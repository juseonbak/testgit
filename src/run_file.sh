# format of inputfile
# timesel pixsel

#!/bin/bash


omdir=/home/Data/OMI
type=ML29
# current directory
pwd=`pwd` 

# Whehter to use qsub or not (1: use qsub disabled 0: command line, else: generate L1L2_fnames.inp)
use_qsub=2
use_omicldo2=1
select_lonlat=0
# set up output and input directory
listdir='./'
listfname='filelist.txt'
lv2dir='./'
solprefix=OMI-Aura_L1-OML1BIRR_
radprefix=OMI-Aura_L1-OML1BRUG_
lv2prefix=OMIO3PROF_$type'-'
cldprefix=OMI-Aura_L2-OMCLDO2_
soldir=/pool/cluster1/xliu/OML1BIRR/
raddir=$omdir/1_OML1BRUG/
clddir=$omdir/2_OML2CLDO2/

nby=1         # number of coadded pixels along the track
subdir=HRES

if [ ! -e $lv2dir ]
then
    mkdir $lv2dir
fi


#for grp in $grplist
#do 
 
cat $listdir$listfname | while read aline
do
    
    aline=($aline)        # turn line into an array
    timesel=${aline[0]}
    pixsel=${aline[1]}' '${aline[2]}' '${aline[3]}' '${aline[4]}
    yr=${timesel:0:4}
    echo $timesel $pixsel

for radfname in `ls $raddir$yr'/'$radprefix$timesel*he4`
do
    orb=${radfname#*-o}  
    orb=${orb:0:5}
    mondayutc=${radfname#*$radprefix}
    mondayutc=${mondayutc:0:14}
    monday=${mondayutc:0:9}
    #echo $mondayutc

   # solfname=`ls $soldir$solprefix$monday*he4`
    cldfname=`ls $clddir$yr'/'$cldprefix*o$orb*he5`
    lv2fname=$lv2dir$lv2prefix'o'$orb

    # Solar radiance could not be found, make up a name in case use backup radiance
    if [ -z $solfname ]
    then 
        #continue
        solfname=$soldir$solprefix$mondayutc'-o'$orb'_v00'$ver'.he4'
    fi

    # fake a cloud file name (not used anyway for the moment)
    if [ -z $cldfname ]
        then
        cldfname=$clddir$cldprefix$mondayutc'-o'$orb'_v00'$ver'.he5'
    fi

    ## cloud data could not be found, skip this orbit
    #if [ -z $cldfname ]
    #then 
    #    continue
    #fi


    orbs[$norb]=$orb
    norb=`expr $norb + 1`

    # generate subdirectory for each orbit
    if [ $use_qsub -eq 1 ] 
        then
        # used by qsub
        echo $orb >> orbs.list
        
        orbdir='o'$orb'_run'
        mkdir -p $orbdir
        cd $orbdir
        
        # cp/link fixed input files
        ln -sf ../SOMIPROF.exe ./SOMIPROF.exe
        cp ../SOMIPROF.pcf ./SOMIPROF.pcf
        mkdir -p INP
        cp -f ../INP/* INP/
        #cp -f ../INP/SOMIPROF.inp INP/SOMIPROF.inp
        #cp -f ../INP/ozprof.inp INP/ozprof.inp
        #cp -f ../INP/o3prof_control.inp INP/o3prof_control.inp
    fi
                   
    # Write input file INP/L1L2_fnames.inp
    # use echo -e to tell that there are control variable new line (\x0a )
    echo -e "$solfname \x0a$radfname \x0a$cldfname \x0a$lv2fname \x0a$pixsel" > INP/L1L2_fnames.inp
    echo $mondayutc $lv2fname
    if [ $select_lonlat -eq 1 ]
    then
        echo -e "T \x0a$latlonsel" >> INP/L1L2_fnames.inp
    fi

    # exit orbit directory
    cd $pwd
            
    if [ $use_qsub -eq 0 ] 
        then            
        SOMIPROF.exe  > $type$orb'.dat' 
        sleep 15
    fi
done


done
#done

