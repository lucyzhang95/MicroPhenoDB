"""Main execution script for MicroPhenoDB data pipeline
Entry point for running the complete data processing and analysis pipeline.
"""

from utils.record_manager import RecordCacheManager


def main():
    """Main execution function."""
    print("Starting MicroPhenoDB data pipeline...")

    # run the pipeline
    records = RecordCacheManager.run_data_pipeline()

    print("🎉 Pipeline completed successfully!")
    print(f"-> Processed {len(records)} unique records")

    try:
        from utils.record_stats_report import generate_record_statistics_report

        generate_record_statistics_report()
    except ImportError as e:
        print(f"‼️ Could not import record_stats_report: {e}")
        print("⏭️ Skipping statistics report generation")


if __name__ == "__main__":
    main()
