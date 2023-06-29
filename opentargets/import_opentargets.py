import argparse
import hail as hl
from hail.utils.java import Env
from ukbb_pan_ancestry import otg_release

bucket = "ukbb-diverse-open-targets-genetics-releases"

# export release=22.09
# gsutil -u ukbb-diversepops-neale -m rsync -rd gs://open-targets-genetics-releases/$release/lut/study-index gs://ukbb-diverse-open-targets-genetics-releases/releases/$release/lut/study-index
# gsutil -u ukbb-diversepops-neale -m rsync -rd gs://open-targets-genetics-releases/$release/lut/variant-index gs://ukbb-diverse-open-targets-genetics-releases/releases/$release/lut/variant-index
# gsutil -u ukbb-diversepops-neale -m rsync -rd gs://open-targets-genetics-releases/$release/v2d gs://ukbb-diverse-open-targets-genetics-releases/releases/$release/v2d
# gsutil -u ukbb-diversepops-neale -m rsync -rd gs://open-targets-genetics-releases/$release/v2d_credset gs://ukbb-diverse-open-targets-genetics-releases/releases/$release/v2d_credset
# gsutil -u ukbb-diversepops-neale -m rsync -rd gs://open-targets-genetics-releases/$release/variant-index gs://ukbb-diverse-open-targets-genetics-releases/releases/$release/variant-index


def main(args):
    spark = Env.spark_session()
    datasets = ["lut/study-index", "lut/variant-index", "v2d", "v2d_credset", "variant-index"]

    for dataset in datasets:
        # up to 22.02.01
        # _reader = spark.read.parquet if dataset == "variant-index" else spark.read.json
        # from 22.09
        _reader = spark.read.parquet
        df = _reader(f"gs://{bucket}/releases/{args.release}/{dataset}/")
        ht = hl.Table.from_spark(df)
        ht = ht.repartition(1000)
        ht = ht.checkpoint(f"gs://{bucket}/releases/{args.release}/{dataset}.ht", overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=str, default=otg_release)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    main(args)
