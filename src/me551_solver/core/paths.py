"""Path resolution utilities.

Works transparently in both development and PyInstaller --onedir bundles.
"""

import os
import sys


def get_data_dir() -> str:
    """Return the absolute path to the me551_solver/data/ directory.

    When running as a PyInstaller bundle, data files are extracted to
    sys._MEIPASS at startup.  In a normal Python environment the data
    directory sits two levels above this file (core/ -> me551_solver/ -> data/).
    """
    if getattr(sys, 'frozen', False):
        # Running as PyInstaller bundle — data is under _MEIPASS
        base_path = sys._MEIPASS  # type: ignore[attr-defined]
    else:
        # Running in development — walk up: core/ -> me551_solver/
        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, 'data')


def get_data_file(filename: str) -> str:
    """Return the full path to *filename* inside the data directory."""
    return os.path.join(get_data_dir(), filename)
