#!/usr/bin/env python3
"""Smoke-test the Qt window and process path inside a PyMOL Python runtime."""

import sys
import tempfile
from pathlib import Path

from pymol.Qt import QtCore, QtWidgets

from pymosaics.core.config import ConfigStore
from pymosaics.gui import PymoSAICSDialog


def main():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    with tempfile.TemporaryDirectory(prefix="PymoSAICS smoke ") as temporary:
        root = Path(temporary)
        forcefields = root / "force fields"
        (forcefields / "top_database").mkdir(parents=True)
        (forcefields / "pot_database").mkdir()
        (forcefields / "top_database" / "test.rtf").write_text("test\n", encoding="utf-8")
        (root / "start.pdb").write_text("END\n", encoding="utf-8")
        parameter_input = root / "parameters.input"
        parameter_input.write_text(
            "# \\mol_parm_file{${PYMOSAICS_FORCEFIELD_DIR}/top_database/test.rtf}\n"
            "# \\pos_init_file{${PROJECT_DIR}/start.pdb}\n"
            "print('PYMOSAICS_QPROCESS_SMOKE')\n",
            encoding="utf-8",
        )

        dialog = PymoSAICSDialog(config_store=ConfigStore(root / "config.json"))
        dialog.executable_edit.setText(sys.executable)
        dialog.forcefield_edit.setText(str(forcefields))
        dialog.input_edit.setText(str(parameter_input))
        dialog.auto_load.setChecked(False)
        dialog.process.finished.connect(application.quit)
        QtCore.QTimer.singleShot(15000, application.quit)
        dialog._start_process()
        application.exec_()
        output = dialog.log_output.toPlainText()
        if "PYMOSAICS_QPROCESS_SMOKE" not in output or "exit code 0 (success)" not in output:
            raise SystemExit("PymoSAICS smoke test failed:\n{}".format(output))
        dialog.close()
        print("PASS: Qt dialog and shell-free QProcess execution")


if __name__ == "__main__":
    main()
