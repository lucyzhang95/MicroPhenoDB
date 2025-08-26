import csv

import chardet


class FileReader:
    """Handles reading and decoding of input files."""

    def _detect_encoding(self, in_file_path):
        """Detects the character encoding of a file."""
        with open(in_file_path, "rb") as f:
            raw = f.read(5000)
        return chardet.detect(raw)["encoding"]

    def read_file(self, in_file_path, has_header=True):
        """Reads a file and yields its lines, handling encoding issues.
        Does not handle an ill-formatted file with an inconsistent number of lines/columns."""
        encoding = self._detect_encoding(in_file_path)
        if encoding == "ascii":
            encoding = "utf-8"
        try:
            with open(in_file_path, "r", encoding=encoding, errors="ignore") as in_f:
                reader = csv.reader(in_f, delimiter="\t")
                if has_header:
                    try:
                        next(reader)
                    except StopIteration:
                        return
                for line in reader:
                    yield line
        except UnicodeDecodeError as e:
            print(f"❗️ Unicode error with {encoding} on file {in_file_path}: {e}")
        except FileNotFoundError:
            print(f"‼️ Error: File not found at {in_file_path}")
