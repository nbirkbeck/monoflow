#
# Script to compute the results for the cloth-rot sequence.
#
dir=${data_dir}/kpop
dmin=1.5
dmax=3
dimg=0.8

vdisp=1
# Good
#dmin=4
#dmax=6
#dimg=0.3
#fimg=0.1
#fbeta=1
#vflow=1
#vdisp_const=0.0025

#0.1, 1.0, 0.0025
fimg=0.16
fbeta=0.9
vdisp_const=0.03
vflow=20

mkdir -p ${dir}/results ${dir}/debug

./discrete --dmin=${dmin} --dmax=${dmax}  \
    --odepth=${dir}/results/discrete-depth-%d.rfi -d 2 --ndisp=120 --diters=2 \
    ${dir}/seq/{0,1,2}/image-%04d.tga

./main --depth=${dir}/results/discrete-depth-%d.rfi \
    --odepth=${dir}/results/vari-depth-%d.rfi -d 4 --vdisp=${vdisp} \
    --vflow=0 --dimg=${dimg} --root=${dir}/debug  \
    ${dir}/seq/{0,1,2}/image-%04d.tga --vdisp_dsmax=2

./main --depth=${dir}/results/vari-depth-%d.rfi \
    --odepth=${dir}/results/sf-depth-%d.rfi \
    --oflow=${dir}/results/sf-flow-%d.rfi -d 4 \
    --vdisp=0  --vflow=${vflow} --fimg=${fimg} --fbeta=${fbeta} \
    --nclosest=-1 --root=${dir}/debug  --vdisp_const=${vdisp_const} --vdisp_dsmax=0 \
    ${dir}/seq/{0,1,2}/image-%04d.tga
