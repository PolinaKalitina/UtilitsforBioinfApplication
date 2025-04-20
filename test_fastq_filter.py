import pytest
import os
from Bio import SeqIO
from UtilitsforBioinfApplication_OOP_version import (
    quality_check,
    length_check,
    gc_check,
    filter_fastq
)


class TestQualityChecks:
    def test_quality_pass(self):
        assert quality_check([30, 30, 30], 20) is True

    def test_quality_fail(self):
        assert quality_check([15, 30, 25], 20) is False

    def test_empty_quality_scores(self):
        assert quality_check([], 20) is False


class TestGCChecks:
    def test_gc_content_edge_cases(self):
        assert gc_check("", (0, 100)) is False
        assert gc_check("GGCC", (100, 100)) is True
        assert gc_check("ATTA", (0, 0)) is True
        assert gc_check("ACGT", (50, 50)) is True


class TestLengthChecks:
    def test_length_filtering(self, tmp_path):
        input_f = tmp_path / "length_test.fastq"
        output_f = tmp_path / "length_out.fastq"

        with open(input_f, 'w') as f:
            f.write("@short\nA\n+\n!\n"
                    "@medium\nACGTACGT\n+\n!!!!!!!!\n"
                    "@long\n" + "A"*100 + "\n+\n" + "!"*100 + "\n")

        import sys
        sys.argv = [
            'script.py',
            '-i', str(input_f),
            '-o', str(output_f),
            '--length-bounds', '4', '10',
            '-q', '0'
        ]

        count = filter_fastq()
        assert count == 1

        records = list(SeqIO.parse(output_f, "fastq"))
        assert len(records) == 1
        assert records[0].id == "medium"


class TestFileOperations:
    @pytest.fixture
    def temp_files(self, tmp_path):
        input_f = tmp_path / "test.fastq"
        output_f = tmp_path / "out.fastq"
        with open(input_f, 'w') as f:
            f.write("@test\nACGT\n+\n!!!!\n")
        yield input_f, output_f
        if output_f.exists():
            os.remove(output_f)

    def test_file_processing(self, temp_files, monkeypatch):
        input_f, output_f = temp_files

        with monkeypatch.context() as m:
            m.setattr('sys.argv', [
                'fastq_filter.py',
                '-i', str(input_f),
                '-o', str(output_f),
                '--gc-bounds', '0', '100',
                '--length-bounds', '1', '10',
                '-q', '0' 
            ])

            count = filter_fastq()

        assert count == 1
        assert os.path.exists(output_f)

    def test_invalid_fastq_handling(self, tmp_path):
        invalid_f = tmp_path / "invalid.fastq"

        with open(invalid_f, 'w') as f:
            f.write("@test\nACGT\nQUAL\n!!!!\n")

        with pytest.raises(ValueError, match="not in FASTQ format"):
            import sys
            sys.argv = [
                'script.py',
                '-i', str(invalid_f),
                '-o', str(tmp_path / "output.fastq"),
                '-q', '0'
            ]
            filter_fastq()


class TestErrorHandling:
    def test_invalid_gc_bounds(self):
        with pytest.raises(ValueError):
            gc_check("ACGT", (110, 120))
