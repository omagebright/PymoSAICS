"""PyMOL-facing window lifecycle."""

_dialog = None


def show_dialog():
    """Create or raise the single PymoSAICS window."""

    global _dialog
    from pymol.Qt import QtWidgets

    from .gui import PymoSAICSDialog

    if _dialog is None:
        _dialog = PymoSAICSDialog(parent=QtWidgets.QApplication.activeWindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
