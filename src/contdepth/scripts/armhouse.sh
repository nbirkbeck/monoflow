#
# Script to compute the results for the cloth-rot sequence.
#
dir=${HOME}/data/contdepth/armhouse5
dmin=2
dmax=4.5
dimg=0.3


# Good
#dmin=4
#dmax=6
#dimg=0.3
#fimg=0.1
#fbeta=1
#vflow=1
vdisp_const=0.0025

fimg=0.125
fbeta=0.6
vdisp=3
vflow=10

mkdir -p ${dir}/results ${dir}/debug

./discrete --dmin=${dmin} --dmax=${dmax}  \
    --odepth=${dir}/results/discrete-depth-%d.rfi -d 3 --ndisp=120 --diters=1 \
    ${dir}/seq/{0,1,2}/image-%04d.tga


./main --depth=${dir}/results/discrete-depth-%d.rfi --odepth=${dir}/results/vari-depth-%d.rfi -d 4 --vdisp=${vdisp} --vdisp_dsmax=2 --nclosest=1  --vflow=0 --dimg=${dimg} --root=${dir}/debug ${dir}/seq/{0,1,2}/image-%04d.tga

./main --depth=${dir}/results/vari-depth-%d.rfi \
    --odepth=${dir}/results/sf-depth-%d.rfi \
    --oflow=${dir}/results/sf-flow-%d.rfi -d 4 \
    --vdisp=0  --vflow=${vflow} --fimg=${fimg} --fbeta=${fbeta} \
    --nclosest=1 --root=${dir}/debug  --vdisp_const=${vdisp_const} --vdisp_dsmax=0 \
    ${dir}/seq/{0,1,2}/image-%04d.tga

