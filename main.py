"""Main execution script for MicroPhenoDB data pipeline
Entry point for running the complete data processing and analysis pipeline.
"""

from utils.record_manager import RecordCacheManager
from utils.record_stats_report import generate_record_statistics_report


def main():
    """Main function to run the data pipeline."""
    print("\nStarting MicroPhenoDB data pipeline...")

    records, jsonl_path = RecordCacheManager.run_data_pipeline()

    print("\nPipeline completed successfully!")
    print(f"-> Processed {len(records)} unique records")
    print(f"-> JSONL export: {jsonl_path}")

    try:
        generate_record_statistics_report()
    except ImportError as e:
        print(f"!!! Skipping statistics report: {e}")


if __name__ == "__main__":
    main()
