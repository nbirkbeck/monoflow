
# Delete old results with
# find monoflow-depth -name results -exec rm -rfv {} \;
# find monoflow-depth -name debug -exec rm -rfv {} \;

d=$(pwd)

data_dir=$(pwd)/monoflow-data
script_dir=$(pwd)/src/contdepth/scripts

# The data is too big to put in github. Grab if it doesn't exist.
if [ ! -f monoflow-data.tar.gz ]; then
    echo Data file doesnt exist!
    echo You can download and install with:
    echo curl http://neilbirkbeck.com/files/monoflow-data.tar.gz -o monoflow-data.tar.gz
    echo tar -xzf monoflow-data.tar.gz
    exit 1
fi

get_varying_baseline() {
    for plane in 2 3 4 6 10; do 
        sh ${script_dir}/plane_moving_camera.sh ${data_dir}/plane_moving_camera/plane_${plane}
    done
}

cd $(pwd)/src/bazel-bin/contdepth

data_dir=${data_dir} sh ${script_dir}/rsphererend.sh
data_dir=${data_dir} sh ${script_dir}/kpop.sh
data_dir=${data_dir} sh ${script_dir}/armhouse.sh

get_varying_baseline

cd $d
