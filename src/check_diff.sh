dir1=./o3prof/
dir2=/home/JDB/OzoneFit/OZBOREAS-OMI/src-201605/o3prof/
dir2=/data/tempo1/Shared/jbak/OzoneFit/OZBOREAS-OMI-zcai-cld/src/o3prof/

files1=$(ls $dir1*)

savefile=$(pwd)/diff.txt
rm -rf $savefile

for file1 in $files1
do
   pos1=$(echo $dir1 | wc -m)-1
   totaln=$(echo $file1 | wc -m)
   pos2=$[$totaln-$pos1+1]
   filename=${file1:$pos1:$pos2}
   file2=$dir2$filename

   echo '============================' >>$savefile
   echo $filename >> $savefile
   diff -d $file1 $file2 >> $savefile

done

cat $savefile
