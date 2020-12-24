#
# Script to compute the results for the rotating sphere
#
dir=/home/osboxes/data/contdepth/rsphererend
dmin=4
dmax=7
dimg=0.25
fimg=0.25
fbeta=0.5
vdisp=3
vdisp_const=0
mkdir -p ${dir}/results ${dir}/debug

./discrete --dmin=${dmin} --dmax=${dmax}  \
   --odepth=${dir}/results/discrete-depth-%d.rfi -d 2 --ndisp=100 --diters=2  --nclosest=2 \
    ${dir}/seq/{0,1,2}/image-%04d.tga

./main --depth=${dir}/results/discrete-depth-%d.rfi \
    --odepth=${dir}/results/vari-depth-%d.rfi -d 4 --vdisp=${vdisp} \
    --vflow=0 --dimg=${dimg} --root=${dir}/debug --nclosest=2   \
    ${dir}/seq/{0,1,2}/image-%04d.tga --vdisp_dsmax=0

./main --depth=${dir}/results/vari-depth-%d.rfi \
    --odepth=${dir}/results/sf-depth-%d.rfi \
    --oflow=${dir}/results/sf-flow-%d.rfi -d 4 \
    --vdisp=0  --vflow=4 --fimg=${fimg} --fbeta=${fbeta} \
    --nclosest=2 --root=${dir}/debug  --vdisp_const=${vdisp_const} --vdisp_dsmax=0 \
    ${dir}/seq/{0,1,2}/image-%04d.tga --vdisp_dsmax=0

