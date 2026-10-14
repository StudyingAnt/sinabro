import csv
import json
from pathlib import Path
import tempfile
import unittest
import contextlib
import io
from sinabro.cli import main
ROOT = Path(__file__).resolve().parents[1]

class CliTests(unittest.TestCase):

    def run_cli(self, output, extra=None):
        args = ['--input', str(ROOT / 'examples/cdsplus.fa'), '--output', str(output), '--signatures', 'SBS2', 'SBS13', '--n-trajectories', '3', '--seed', '42']
        with contextlib.redirect_stdout(io.StringIO()):
            return main(args + (extra or []))

    def test_smoke_csv_and_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'result.csv'
            self.assertEqual(self.run_cli(output), 0)
            with output.open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 12)
            self.assertTrue(all((r['terminal_reason'] == 'nonsynonymous' for r in rows)))
            manifest = json.loads(output.with_suffix('.manifest.json').read_text())
            self.assertEqual(manifest['rows'], 12)
            self.assertEqual(manifest['incomplete'], 0)
            self.assertEqual(len(manifest['input_sha256']), 64)

    def test_repeatability(self):
        with tempfile.TemporaryDirectory() as directory:
            a = Path(directory) / 'a.csv'
            b = Path(directory) / 'b.csv'
            self.run_cli(a)
            self.run_cli(b)
            self.assertEqual(a.read_bytes(), b.read_bytes())

    def test_blosum_csv_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'result.csv'
            self.run_cli(output, ['--evaluation', 'blosum', '--threshold', '20'])
            with output.open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertTrue(all((r['evaluation'] == 'blosum_distance' for r in rows)))

    def test_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'result.csv'
            output.write_text('keep')
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.run_cli(output)
            self.assertEqual(output.read_text(), 'keep')

    def test_invalid_cli_count(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'result.csv'
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.run_cli(output, ['--n-trajectories', '0'])
            self.assertFalse(output.exists())
if __name__ == '__main__':
    unittest.main()
