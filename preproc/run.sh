set -x
cargo run --release "../data/targets_$(hostname)" > "../logs/logs-prepoc-$(hostname)-$(date +%F-%H-%M)"
