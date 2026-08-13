import argparse

from campaign_reader import build_feature_table


DEFAULT_ACX = (
    "/global/homes/b/bhowmic/campaign-store/IFE/rhino.acx"
)

DEFAULT_QUERY_DIR = (
    "/global/u1/b/bhowmic/Projects/IFE_AmSC/RHINO/campaign/queries"
)

DEFAULT_CAMPAIGN_STORE = (
    "/global/homes/b/bhowmic/campaign-store"
)

DEFAULT_ARCHIVE_NAME = "%"


def main():
    parser = argparse.ArgumentParser(
        description="Build RHINO feature table."
    )

    parser.add_argument(
        "--output",
        default="rhino_features.csv",
        help="Output CSV file",
    )

    args = parser.parse_args()

    df = build_feature_table(
        acx_path=DEFAULT_ACX,
        query_dir=DEFAULT_QUERY_DIR,
        campaign_store=DEFAULT_CAMPAIGN_STORE,
        archive_name=DEFAULT_ARCHIVE_NAME,
    )

    print(df.head())
    print(df.shape)

    print(df.isna().sum())
    print(df.describe())

    df.to_csv(args.output, index=False)

    print(f"\nSaved feature table to: {args.output}")


if __name__ == "__main__":
    main()