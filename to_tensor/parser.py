import argparse

def main():
    parser = argparse.ArgumentParser(description='Tool for generating TFRecords from bam files')
    parser.add_argument("bar", type=int)
    parser.add_argument("--foo", type=int)
    args= parser.parse_args()
    print(args.bar)


if __name__ == "__main__":
    main()
