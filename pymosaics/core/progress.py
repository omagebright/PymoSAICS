"""Streaming MOSAICS progress and ETA estimation."""

import math
import re
from dataclasses import dataclass
from typing import Optional


STEP_PATTERN = re.compile(r"Report\s+energies\s+at\s+step\s+([0-9]+)", re.IGNORECASE)
TOTAL_STEP_PATTERN = re.compile(r"\\total_step_mc\{\s*([0-9]+)\s*\}")


@dataclass(frozen=True)
class RunProgress:
    current_step: int
    total_steps: int
    fraction: float
    elapsed_seconds: float
    steps_per_second: Optional[float]
    remaining_seconds: Optional[float]


def total_steps_from_input(text: str) -> int:
    matches = [
        match
        for raw in text.splitlines()
        if not raw.lstrip().startswith(("#", "!"))
        for match in TOTAL_STEP_PATTERN.finditer(raw)
    ]
    if not matches:
        raise ValueError("MOSAICS input does not define total_step_mc")
    return int(matches[-1].group(1))


def format_duration(seconds: Optional[float]) -> str:
    if seconds is None or not math.isfinite(seconds) or seconds < 0:
        return "--:--:--"
    rounded = max(0, int(seconds))
    hours, remainder = divmod(rounded, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02d}:{:02d}:{:02d}".format(hours, minutes, seconds)


class RunProgressTracker:
    """Consume arbitrary stdout chunks and expose a monotonic run estimate."""

    def __init__(self, total_steps: int, started_at: float):
        if total_steps < 0:
            raise ValueError("total_steps must be non-negative")
        self.total_steps = int(total_steps)
        self.started_at = float(started_at)
        self.current_step = 0
        self._tail = ""
        self._rate = None
        self._estimated_finish = None

    def ingest(self, text: str, now: float) -> RunProgress:
        combined = self._tail + text
        for match in STEP_PATTERN.finditer(combined):
            observed = int(match.group(1))
            if self.total_steps:
                observed = min(observed, self.total_steps)
            if observed > self.current_step:
                self.current_step = observed
                elapsed = max(0.0, float(now) - self.started_at)
                if elapsed > 0:
                    self._rate = self.current_step / elapsed
                    if self.total_steps > 0 and self._rate > 0:
                        self._estimated_finish = float(now) + (
                            self.total_steps - self.current_step
                        ) / self._rate
        self._tail = combined[-160:]
        return self.snapshot(now)

    def snapshot(self, now: float) -> RunProgress:
        elapsed = max(0.0, float(now) - self.started_at)
        fraction = (
            min(1.0, self.current_step / self.total_steps)
            if self.total_steps > 0
            else (1.0 if self.current_step > 0 else 0.0)
        )
        remaining = (
            max(0.0, self._estimated_finish - float(now))
            if self._estimated_finish is not None
            else None
        )
        return RunProgress(
            self.current_step,
            self.total_steps,
            fraction,
            elapsed,
            self._rate,
            remaining,
        )
