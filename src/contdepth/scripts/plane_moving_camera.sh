#
# Script to compute the results for the cloth-rot sequence.
#
dir=/tmp/plane_4
dmin=3.5
dmax=5
dimg=1.0
fimg=2.0
vdisp=2
vdisp_const=0

if [ $# == "0" ];
then
    echo "Give directory.";
else
    dir=$1
fi;

echo "Directory is ${dir}";

mkdir -p ${dir}/results ${dir}/debug


fimg=1.2

echo vdisp: $vdisp fimg: $fimg;

rdir=${dir}/results
mkdir -p ${rdir}

# Uncomment these if you want to run everything again.

./discrete --dmin=${dmin} --dmax=${dmax}  \
   --odepth=${dir}/results/discrete-depth-%d.rfi -d 2 --ndisp=100 --diters=2  --nclosest=-1 \
    ${dir}/seq/{0,1,2}/image-%04d.png

./main --depth=${dir}/results/discrete-depth-%d.rfi \
    --odepth=${dir}/results/vari-depth-%d.rfi -d 4 --vdisp=${vdisp} \
    --vflow=0 --dimg=${dimg} --root=${dir}/debug --nclosest=-1   \
    ${dir}/seq/{0,1,2}/image-%04d.png --vdisp_dsmax=0

./main --depth=${dir}/results/vari-depth-%d.rfi \
    --odepth=${dir}/results/sf-depth-%d.rfi \
    --oflow=${dir}/results/sf-flow-%d.rfi -d 3 \
    --vdisp=0  --vflow=1 --fimg=${fimg} \
    --nclosest=-1 --root=${dir}/debug  --vdisp_const=${vdisp_const} \
    ${dir}/seq/{0,1,2}/image-%04d.png --vdisp_dsmax=0
