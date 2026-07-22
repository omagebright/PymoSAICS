#!/usr/bin/env python3
"""Smoke-test the Qt window and process path inside a PyMOL Python runtime."""

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pymol.Qt import QtCore, QtWidgets

from pymosaics.core.analysis import latest_log
from pymosaics.core.catalog import FORCE_FIELD_PROFILES
from pymosaics.core.config import ConfigStore
from pymosaics.gui import PymoSAICSDialog, RegionEditor


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
        dialog.resize(900, 700)
        dialog.show()
        application.processEvents()
        force_field_index = dialog.setup_forcefield_combo.findData("ff14sb-protein")
        dialog.setup_forcefield_combo.setCurrentIndex(force_field_index)
        if dialog.forcefield_combo.currentData() != "ff14sb-protein":
            raise SystemExit("Setup and Build force-field selectors are not synchronized")
        if dialog.setup_forcefield_combo.count() != len(FORCE_FIELD_PROFILES):
            raise SystemExit("The Setup selector does not expose every force-field profile")
        for profile_file in FORCE_FIELD_PROFILES[force_field_index].all_paths():
            if profile_file.name not in dialog.input_preview.toPlainText():
                raise SystemExit("The selected profile did not populate {}".format(profile_file.name))
        scroll = dialog.findChild(QtWidgets.QScrollArea, "buildScroll")
        content = dialog.findChild(QtWidgets.QWidget, "buildScrollContent")
        if content.width() > scroll.viewport().width():
            raise SystemExit("The Build form is wider than its viewport")
        if content.palette().window().color().name().lower() != "#0d252f":
            raise SystemExit("The Build form did not receive the deterministic dark palette")
        for combo in dialog.findChildren(QtWidgets.QComboBox):
            popup = combo.view()
            if not isinstance(popup, QtWidgets.QListView):
                raise SystemExit("A combo box is still using a host-native popup")
            if popup.palette().base().color().name().lower() != "#102c36":
                raise SystemExit("A combo popup inherited a light host palette")
            if popup.palette().text().color().name().lower() != "#edf7f8":
                raise SystemExit("A combo popup does not use readable foreground text")
            if "#789198" not in popup.styleSheet().lower():
                raise SystemExit("Disabled combo choices do not have a readable muted color")
            if popup.minimumWidth() < 220:
                raise SystemExit("A combo popup is too narrow to identify its choices")
        region = RegionEditor(("A:1", "A:2", "B:3"), "demo", parent=dialog)
        if region._current_settings().dependency_type != "independent":
            raise SystemExit("The region editor exposed an unsupported dependency type")
        if region._current_settings().centers != ("A:1",):
            raise SystemExit("A new region does not start with a valid rotation center")
        if not region.save_button.isEnabled() or "\\ncenter{1}" not in region.preview.toPlainText():
            raise SystemExit("The default graphical region does not generate a valid preview")
        region._set_all_residues(False)
        if region.save_button.isEnabled():
            raise SystemExit("An invalid empty region can be accepted")
        for combo in region.findChildren(QtWidgets.QComboBox):
            if combo.view().palette().base().color().name().lower() != "#102c36":
                raise SystemExit("A region-editor combo popup inherited a light host palette")
        region.close()
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
