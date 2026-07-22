"""PymoSAICS: a transparent MOSAICS workbench inside PyMOL."""

__version__ = "0.3.0"


def __init_plugin__(app=None):
    """Register the plugin in current Qt-based PyMOL releases."""

    from pymol.plugins import addmenuitemqt

    from .plugin import show_dialog

    addmenuitemqt("PymoSAICS", show_dialog)
