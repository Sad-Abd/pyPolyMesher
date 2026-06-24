from __future__ import annotations

import os
import sys
from pathlib import Path

# Ensure tests run headless and never open GUI windows before matplotlib import.
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg", force=True)

# Add the src layout to the import path for tests.
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


@pytest.fixture(autouse=True)
def _deterministic_seed():
    """Reset NumPy's legacy RNG before every test for deterministic behavior."""
    np.random.seed(12345)
    yield


@pytest.fixture
def rng():
    """Return a reproducible NumPy Generator for tests that prefer local RNG state."""
    return np.random.default_rng(12345)


@pytest.fixture
def np_helpers():
    """Expose a small namespace of NumPy helpers for readability in tests."""
    return {
        "array": np.array,
        "asarray": np.asarray,
        "column_stack": np.column_stack,
        "linspace": np.linspace,
        "meshgrid": np.meshgrid,
    }


@pytest.fixture
def temp_dxf_path(tmp_path):
    """Provide a temporary directory for DXF-related tests."""
    return tmp_path


@pytest.fixture(autouse=True)
def _no_gui_backend(monkeypatch):
    """Prevent plot windows from opening during tests."""
    monkeypatch.setenv("MPLBACKEND", "Agg")
    matplotlib.use("Agg", force=True)
    yield
