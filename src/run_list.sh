#!/bin/bash
homedir=/JDB
omdir=/OMI 
type=ML29
# current directory
pwd=`pwd` 

# Whehter to use qsub or not (1: use qsub disabled 0: command line, else: generate L1L2_fnames.inp)
use_qsub=2
use_omicldo2=1

# set up output and input directory
listdir='./'
listfname='fail.txt'
#listfname='fail.txt'
solprefix=OMI-Aura_L1-OML1BIRR_
radprefix=OMI-Aura_L1-OML1BRUG_
lv2prefix=OMIO3PROF_$type'-'
soldir=/pool/cluster1/xliu/OML1BIRR/
raddir=$omdir/1_OML1BRUG/
clddir=$omdir/2_OML2CLDO2/
cldprefix=OMI-Aura_L2-OMCLDO2_
lv2dir='./'
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
    date=${aline[0]}
    orb=${aline[1]}
    sline1=${aline[2]}
    pix=${aline[3]}
    sline2=`expr $sline1 + $nby - 1`
    #sline1=$sline2
    pix1=`expr $pix`
    pix2=`expr $pix + 1 `
    pixsel="$sline1 $sline2 $pix1 $pix2"
    timesel=$date'*-o*'$orb
    year=${date:0:4}
    #echo $timesel
    lv2dirsub=$lv2dir/
    raddirsub=$raddir$year/
    clddirsub=$clddir$year/
    # initialized the orbits to be arrays
    # note the index starts from zero
    orbs[0]=''
    norb=0

    for radfname in `ls $raddirsub$radprefix$timesel*he4`
    do
        echo $grp $radfname $sline1 $pix1
        orb=${radfname#*-o}  
        orb=${orb:0:5}
        #echo $orb
        mondayutc=${radfname#*$radprefix}
        mondayutc=${mondayutc:0:14}
        monday=${mondayutc:0:9}
        #echo $mondayutc
        
    #    solfname=`ls $soldir$solprefix$monday*he4`
        cldfname=`ls $clddirsub$cldprefix*o$orb*he5`
        lv2fname=$lv2dirsub$lv2prefix'o'$orb
        
        # Solar radiance could not be found, make up a name in case use backup radiance
    #    if [ -z $solfname ]
    #    then 
    #        #continue
            solfname=$soldir$solprefix$mondayutc'-o'$orb'_v00'$ver'.he4'
    #    fi

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
            SOMIPROF'.exe'  > $type'.dat' 
            cat $type'.dat'
        fi
    done  
done
#done

