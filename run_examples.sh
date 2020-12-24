
# Delete old results with
# find ~/data/contdepth -name results -exec rm -rfv {} \;
# find ~/data/contdepth -name debug -exec rm -rfv {} \;

d=$(pwd)

script_dir=$(pwd)/contdepth/scripts
cd $(pwd)/bazel-bin/contdepth

sh ${script_dir}/rsphererend.sh
sh ${script_dir}/kpop.sh
sh ${script_dir}/plane_moving_camera.sh
sh ${script_dir}/armhouse.sh

cd $d
