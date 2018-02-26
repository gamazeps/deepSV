set -x
cargo run --release "cluster_conf.json" "../data/targets_$(hostname)" > "../logs/logs-prepoc-$(hostname)-$(date +%F-%H-%M)"
