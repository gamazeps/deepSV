cd dev/deepSV/to_tensor
source to_tensor/bin/activate
python merge_hdf5.py cluster_merge_conf.json ../data/targets_$(hostname) ../logs/LOL_$(hostname) &> ../logs/fuck_$(hostname)
