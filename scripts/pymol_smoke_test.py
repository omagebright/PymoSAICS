#!/usr/bin/env python3
"""Smoke-test the Qt window and process path inside a PyMOL Python runtime."""

import sys
import tempfile
from pathlib import Path

from pymol.Qt import QtCore, QtWidgets

from pymosaics.core.analysis import latest_log
from pymosaics.core.config import ConfigStore
from pymosaics.gui import PymoSAICSDialog


def main():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    with tempfile.TemporaryDirectory(prefix="PymoSAICS smoke ") as temporary:
        root = Path(temporary)
        forcefields = root / "forcefield"
        forcefields.mkdir()
        (forcefields / "test.rtf").write_text("test\n", encoding="utf-8")
        (root / "start.pdb").write_text("END\n", encoding="utf-8")
        parameter_input = root / "parameters.input"
        parameter_input.write_text(
            "# \\mol_parm_file{forcefield/test.rtf}\n"
            "# \\pos_init_file{start.pdb}\n"
            "print('PYMOSAICS_QPROCESS_SMOKE')\n",
            encoding="utf-8",
        )

        dialog = PymoSAICSDialog(config_store=ConfigStore(root / "config.json"))
        dialog.runtime_combo.setCurrentIndex(dialog.runtime_combo.findData("custom"))
        dialog.custom_executable_edit.setText(sys.executable)
        dialog.project_edit.setText(str(root))
        dialog.input_combo.setEditText(str(parameter_input))
        dialog.auto_load.setChecked(False)
        dialog.process.finished.connect(application.quit)
        QtCore.QTimer.singleShot(15000, application.quit)
        dialog._start_process()
        application.exec_()
        output = dialog.log_output.toPlainText()
        if "PYMOSAICS_QPROCESS_SMOKE" not in output or "exit code 0 (success)" not in output:
            raise SystemExit("PymoSAICS smoke test failed:\n{}".format(output))
        log = latest_log(root)
        if log is None or "exit code 0 (success)" not in log.read_text(encoding="utf-8"):
            raise SystemExit("PymoSAICS did not persist the process completion status")
        dialog.close()
        print("PASS: Qt dialog and shell-free QProcess execution")


if __name__ == "__main__":
    main()
