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

    def test_file_processing(self, temp_files):
        input_f, output_f = temp_files
        count = filter_fastq(input_f, output_f)
        assert count == 1
        assert os.path.exists(output_f)


class TestErrorHandling:
    def test_invalid_gc_bounds(self):
        with pytest.raises(ValueError):
            gc_check("ACGT", (110, 120))
