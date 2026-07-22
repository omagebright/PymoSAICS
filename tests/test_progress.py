import unittest

from pymosaics.core.progress import RunProgressTracker, format_duration, total_steps_from_input


class RunProgressTests(unittest.TestCase):
    def test_streamed_report_steps_produce_percentage_and_eta(self):
        tracker = RunProgressTracker(1000, started_at=10.0)
        tracker.ingest(">>Report energies at st", now=20.0)
        tracker.ingest("ep 250<<\n", now=20.0)
        progress = tracker.snapshot(now=30.0)
        self.assertEqual(progress.current_step, 250)
        self.assertAlmostEqual(progress.fraction, 0.25)
        self.assertAlmostEqual(progress.steps_per_second, 25.0)
        self.assertAlmostEqual(progress.remaining_seconds, 20.0)

    def test_progress_is_monotonic_and_clamped_to_total(self):
        tracker = RunProgressTracker(100, started_at=0.0)
        tracker.ingest("Report energies at step 80", now=8.0)
        tracker.ingest("Report energies at step 20", now=9.0)
        tracker.ingest("Report energies at step 120", now=10.0)
        self.assertEqual(tracker.snapshot(now=10.0).current_step, 100)

    def test_total_steps_and_duration_are_reviewable(self):
        self.assertEqual(
            total_steps_from_input(
                "# \\total_step_mc{9}\n\\total_step_mc{100000}\n"
            ),
            100000,
        )
        self.assertEqual(format_duration(3661.2), "01:01:01")
        self.assertEqual(format_duration(None), "--:--:--")


if __name__ == "__main__":
    unittest.main()
