#!/usr/bin/env python3
"""Smoke-test the Qt window and process path inside a PyMOL Python runtime."""

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pymol.Qt import QtCore, QtWidgets

from pymosaics.core.analysis import latest_log, read_text_file
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
        alternate_input = root / "alternate.inp"
        alternate_input.write_text(parameter_input.read_text(encoding="utf-8"), encoding="utf-8")

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

        dialog.project_edit.setText(str(root))
        dialog._scan_inputs()
        discovered_inputs = {
            Path(dialog.input_combo.itemText(index)).name
            for index in range(dialog.input_combo.count())
        }
        if discovered_inputs != {"parameters.input", "alternate.inp"}:
            raise SystemExit("The project directory did not automatically expose every input file")
        empty_project = root / "empty-project"
        empty_project.mkdir()
        dialog.project_edit.setText(str(empty_project))
        dialog._scan_inputs()
        if Path(dialog.input_combo.currentText()) != (empty_project / "mcmc.input").resolve():
            raise SystemExit("An empty project did not default to mcmc.input")
        dialog.project_edit.setText(str(root))
        dialog._scan_inputs(preferred=parameter_input)

        dialog.tabs.setCurrentWidget(dialog.run_tab)
        application.processEvents()
        run_widgets = (
            dialog.validation_output,
            dialog.command_preview,
            dialog.run_button,
            dialog.auto_load,
            dialog.log_output,
        )
        run_rectangles = []
        for widget in run_widgets:
            top_left = widget.mapTo(dialog.run_tab, QtCore.QPoint(0, 0))
            run_rectangles.append(QtCore.QRect(top_left, widget.size()))
        for index, first in enumerate(run_rectangles):
            for second in run_rectangles[index + 1 :]:
                if first.intersects(second):
                    raise SystemExit("Run-tab controls overlap at the minimum window size")
        if dialog.validation_output.height() < 72 or dialog.command_preview.height() < 72:
            raise SystemExit("Run-tab diagnostic panels are too short to inspect")
        if dialog.log_output.height() < 130:
            raise SystemExit("The live MOSAICS output panel is too short to use")
        dialog.tabs.setCurrentWidget(dialog.analysis_tab)
        dialog.analysis_pages.setCurrentIndex(0)
        application.processEvents()
        if dialog.energy_plot.width() < 300 or dialog.energy_plot.height() < 210:
            raise SystemExit("The Analysis energy plot is clipped at the minimum window size")
        if dialog.acceptance_table.width() < 250 or dialog.acceptance_table.height() < 210:
            raise SystemExit("The Analysis acceptance table is clipped at the minimum window size")
        if dialog.acceptance_table.horizontalScrollBar().maximum() != 0:
            raise SystemExit("The Analysis acceptance columns require horizontal scrolling")
        dialog.analysis_pages.setCurrentIndex(1)
        application.processEvents()
        if dialog.trajectory_combo.width() < 200:
            raise SystemExit("The structural-landscape trajectory selector is clipped")
        if dialog.landscape_plot.width() < 300 or dialog.landscape_plot.height() < 170:
            raise SystemExit("The structural-landscape plot is too small to use")
        if dialog.representative_list.width() < 250 or dialog.representative_list.height() < 100:
            raise SystemExit("The structural representative list is too small to use")
        dialog.analysis_pages.setCurrentIndex(2)
        application.processEvents()
        if dialog.output_list.height() < 260:
            raise SystemExit("The Analysis files list is too small to use")
        dialog.tabs.setCurrentWidget(dialog.build_tab)
        application.processEvents()

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
        if log.parent != (root / "logs").resolve():
            raise SystemExit("The run log is hidden instead of being stored in project/logs")
        if read_text_file(log, maximum_bytes=None).rstrip("\n") != dialog.log_output.toPlainText().rstrip("\n"):
            raise SystemExit("The Run tab does not display the complete persisted log")
        dialog.close()
        print("PASS: Qt dialog and shell-free QProcess execution")


if __name__ == "__main__":
    main()
