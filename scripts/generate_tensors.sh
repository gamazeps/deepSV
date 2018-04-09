cd dev/deepSV/to_tensor
source to_tensor/bin/activate
python generate_tensors.py cluster_conf.json ../data/targets_$(hostname) ../logs/LOL_$(hostname) &> ../logs/fuck_$(hostname)
