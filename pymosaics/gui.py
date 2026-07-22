"""Qt user interface loaded through PyMOL's supported compatibility layer."""

import os
from pathlib import Path
from typing import Optional

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from .core.config import ConfigError, ConfigStore
from .core.models import PreparedRun, RuntimeConfig
from .core.project import (
    PreparationError,
    discover_outputs,
    format_diagnostics,
    planned_parameter_input,
    prepare_run,
    validate_project,
)
from .core.runtime import has_errors, validate_runtime


class PymoSAICSDialog(QtWidgets.QDialog):
    """Configure, validate, run, stop, and inspect a MOSAICS job."""

    def __init__(self, parent=None, config_store: Optional[ConfigStore] = None):
        super().__init__(parent)
        self.setWindowTitle("PymoSAICS")
        self.resize(820, 650)
        self.config_store = config_store or ConfigStore()
        self.process = QtCore.QProcess(self)
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read_process_output)
        self.process.finished.connect(self._process_finished)
        self.process.errorOccurred.connect(self._process_error)
        self._active_run: Optional[PreparedRun] = None
        self._log_handle = None

        self.tabs = QtWidgets.QTabWidget()
        self.run_tab = QtWidgets.QWidget()
        self.setup_tab = QtWidgets.QWidget()
        self.about_tab = QtWidgets.QWidget()
        self.tabs.addTab(self.run_tab, "Run")
        self.tabs.addTab(self.setup_tab, "Setup")
        self.tabs.addTab(self.about_tab, "About")

        outer = QtWidgets.QVBoxLayout(self)
        outer.addWidget(self.tabs)
        self._build_run_tab()
        self._build_setup_tab()
        self._build_about_tab()
        self._load_configuration()

    def _path_row(self, line_edit, button_text, callback):
        container = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(line_edit, 1)
        button = QtWidgets.QPushButton(button_text)
        button.clicked.connect(callback)
        layout.addWidget(button)
        return container

    def _build_setup_tab(self):
        layout = QtWidgets.QVBoxLayout(self.setup_tab)
        explanation = QtWidgets.QLabel(
            "PymoSAICS does not download or redistribute MOSAICS. Select your "
            "own platform-compatible executable and force-field directory."
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)

        form = QtWidgets.QFormLayout()
        self.executable_edit = QtWidgets.QLineEdit()
        self.forcefield_edit = QtWidgets.QLineEdit()
        self.workspace_edit = QtWidgets.QLineEdit()
        form.addRow(
            "MOSAICS executable:",
            self._path_row(self.executable_edit, "Browse…", self._browse_executable),
        )
        form.addRow(
            "Force-field directory:",
            self._path_row(self.forcefield_edit, "Browse…", self._browse_forcefields),
        )
        form.addRow(
            "Default workspace (optional):",
            self._path_row(self.workspace_edit, "Browse…", self._browse_workspace),
        )
        layout.addLayout(form)

        self.setup_status = QtWidgets.QPlainTextEdit()
        self.setup_status.setReadOnly(True)
        self.setup_status.setMaximumHeight(150)
        layout.addWidget(QtWidgets.QLabel("Validation results"))
        layout.addWidget(self.setup_status)
        save_button = QtWidgets.QPushButton("Validate and save")
        save_button.clicked.connect(self._save_configuration)
        layout.addWidget(save_button)
        layout.addStretch(1)

    def _build_run_tab(self):
        layout = QtWidgets.QVBoxLayout(self.run_tab)
        explanation = QtWidgets.QLabel(
            "Choose a manually reviewed MOSAICS parameter input. PymoSAICS "
            "validates known input-file references before starting the process."
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)

        form = QtWidgets.QFormLayout()
        self.input_edit = QtWidgets.QLineEdit()
        self.output_edit = QtWidgets.QLineEdit()
        form.addRow(
            "Parameter input:",
            self._path_row(self.input_edit, "Browse…", self._browse_input),
        )
        form.addRow(
            "Output PDB (optional):",
            self._path_row(self.output_edit, "Browse…", self._browse_output),
        )
        layout.addLayout(form)

        placeholder_help = QtWidgets.QLabel(
            "Portable inputs may use ${PYMOSAICS_FORCEFIELD_DIR} and ${PROJECT_DIR}. "
            "The source input is never modified."
        )
        placeholder_help.setWordWrap(True)
        layout.addWidget(placeholder_help)

        self.validation_output = QtWidgets.QPlainTextEdit()
        self.validation_output.setReadOnly(True)
        self.validation_output.setMaximumHeight(125)
        layout.addWidget(QtWidgets.QLabel("Validation results"))
        layout.addWidget(self.validation_output)

        self.command_preview = QtWidgets.QPlainTextEdit()
        self.command_preview.setReadOnly(True)
        self.command_preview.setMaximumHeight(95)
        self.command_preview.setPlaceholderText("The exact program, argument, and working directory appear here before execution.")
        layout.addWidget(QtWidgets.QLabel("Execution plan"))
        layout.addWidget(self.command_preview)

        controls = QtWidgets.QHBoxLayout()
        validate_button = QtWidgets.QPushButton("Validate")
        validate_button.clicked.connect(self._validate_current_project)
        self.run_button = QtWidgets.QPushButton("Run MOSAICS")
        self.run_button.clicked.connect(self._start_process)
        self.stop_button = QtWidgets.QPushButton("Stop")
        self.stop_button.clicked.connect(self._stop_process)
        self.stop_button.setEnabled(False)
        load_button = QtWidgets.QPushButton("Load output in PyMOL")
        load_button.clicked.connect(self._load_output)
        controls.addWidget(validate_button)
        controls.addWidget(self.run_button)
        controls.addWidget(self.stop_button)
        controls.addWidget(load_button)
        controls.addStretch(1)
        layout.addLayout(controls)

        self.auto_load = QtWidgets.QCheckBox("Load the newest output PDB after a successful run")
        self.auto_load.setChecked(True)
        layout.addWidget(self.auto_load)

        self.log_output = QtWidgets.QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        layout.addWidget(QtWidgets.QLabel("MOSAICS output"))
        layout.addWidget(self.log_output, 1)

    def _build_about_tab(self):
        layout = QtWidgets.QVBoxLayout(self.about_tab)
        label = QtWidgets.QLabel(
            "<h2>PymoSAICS 0.1.1</h2>"
            "<p>A transparent, cross-platform PyMOL interface for an external "
            "MOSAICS runtime.</p>"
            "<p>This plugin contains no MOSAICS executable, force-field data, "
            "or automatic scientific parameter generator.</p>"
            "<p>Supported interface: current Qt-based PyMOL on Windows, macOS, "
            "and Linux. Simulation compatibility depends on the separately "
            "installed MOSAICS runtime.</p>"
        )
        label.setWordWrap(True)
        layout.addWidget(label)
        layout.addStretch(1)

    def _configuration_from_fields(self) -> RuntimeConfig:
        workspace = self.workspace_edit.text().strip()
        return RuntimeConfig(
            executable=Path(self.executable_edit.text().strip()).expanduser(),
            forcefield_directory=Path(self.forcefield_edit.text().strip()).expanduser(),
            default_workspace=Path(workspace).expanduser() if workspace else None,
        )

    def _load_configuration(self):
        try:
            config = self.config_store.load()
        except ConfigError as exc:
            self.setup_status.setPlainText("ERROR: {}".format(exc))
            self.tabs.setCurrentWidget(self.setup_tab)
            return
        if config is None:
            self.setup_status.setPlainText("No runtime has been configured yet.")
            self.tabs.setCurrentWidget(self.setup_tab)
            return
        self.executable_edit.setText(str(config.executable))
        self.forcefield_edit.setText(str(config.forcefield_directory))
        self.workspace_edit.setText(str(config.default_workspace or ""))
        self.setup_status.setPlainText(format_diagnostics(validate_runtime(config)))

    def _save_configuration(self):
        config = self._configuration_from_fields()
        diagnostics = validate_runtime(config)
        self.setup_status.setPlainText(format_diagnostics(diagnostics))
        if has_errors(diagnostics):
            return
        try:
            self.config_store.save(config)
        except ConfigError as exc:
            self.setup_status.appendPlainText("ERROR: {}".format(exc))
            return
        self.setup_status.appendPlainText("Saved to {}".format(self.config_store.path))

    def _browse_executable(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select MOSAICS executable")
        if path:
            self.executable_edit.setText(path)

    def _browse_forcefields(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select force-field directory")
        if path:
            self.forcefield_edit.setText(path)

    def _browse_workspace(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select default workspace")
        if path:
            self.workspace_edit.setText(path)

    def _initial_browse_directory(self) -> str:
        workspace = self.workspace_edit.text().strip()
        return workspace if workspace and Path(workspace).is_dir() else str(Path.home())

    def _browse_input(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select MOSAICS parameter input",
            self._initial_browse_directory(),
            "MOSAICS input (*.input *.txt);;All files (*)",
        )
        if path:
            self.input_edit.setText(path)
            self._validate_current_project()

    def _browse_output(self):
        initial = self._project_directory()
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select output PDB", str(initial), "PDB files (*.pdb);;All files (*)"
        )
        if path:
            self.output_edit.setText(path)

    def _project_directory(self) -> Path:
        input_text = self.input_edit.text().strip()
        return Path(input_text).expanduser().parent if input_text else Path(self._initial_browse_directory())

    def _validate_current_project(self):
        parameter_input = Path(self.input_edit.text().strip()).expanduser()
        config = self._configuration_from_fields()
        diagnostics = validate_project(parameter_input, config)
        self.validation_output.setPlainText(format_diagnostics(diagnostics))
        self._update_command_preview(parameter_input, config)
        return diagnostics

    def _update_command_preview(self, parameter_input: Path, config: RuntimeConfig):
        try:
            actual_input = planned_parameter_input(parameter_input, config)
        except (OSError, UnicodeError):
            actual_input = parameter_input
        self.command_preview.setPlainText(
            "Program: {}\nArgument 1: {}\nWorking directory: {}".format(
                config.executable.expanduser().resolve(),
                actual_input,
                parameter_input.expanduser().resolve().parent,
            )
        )

    def _start_process(self):
        if self.process.state() != QtCore.QProcess.NotRunning:
            return
        diagnostics = self._validate_current_project()
        if has_errors(diagnostics):
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Correct the validation errors before running.")
            return

        try:
            active_run = prepare_run(
                Path(self.input_edit.text().strip()), self._configuration_from_fields()
            )
        except (PreparationError, OSError) as exc:
            QtWidgets.QMessageBox.critical(self, "PymoSAICS", str(exc))
            return

        self._active_run = active_run
        try:
            self._log_handle = active_run.log_file.open("ab")
        except OSError as exc:
            QtWidgets.QMessageBox.critical(self, "PymoSAICS", "Cannot create run log: {}".format(exc))
            self._active_run = None
            return

        self.log_output.clear()
        preamble = (
            "Starting MOSAICS\n"
            "Program: {}\n"
            "Input: {}\n"
            "Working directory: {}\n"
            "Log: {}\n\n"
        ).format(
            active_run.command[0],
            active_run.command[1],
            active_run.working_directory,
            active_run.log_file,
        )
        self.log_output.setPlainText(preamble)
        self._log_handle.write(preamble.encode("utf-8"))
        self._log_handle.flush()
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.process.setWorkingDirectory(str(active_run.working_directory))
        self.process.setProgram(active_run.command[0])
        self.process.setArguments(list(active_run.command[1:]))
        self.process.start()

    def _read_process_output(self):
        data = bytes(self.process.readAllStandardOutput())
        if not data:
            return
        if self._log_handle is not None:
            self._log_handle.write(data)
            self._log_handle.flush()
        text = data.decode("utf-8", errors="replace")
        self.log_output.moveCursor(QtGui.QTextCursor.End)
        self.log_output.insertPlainText(text)
        self.log_output.moveCursor(QtGui.QTextCursor.End)

    def _process_error(self, error):
        message = "\nProcess error: {}".format(self.process.errorString())
        self.log_output.appendPlainText(message)
        if self._log_handle is not None:
            self._log_handle.write((message + "\n").encode("utf-8"))
            self._log_handle.flush()
        if self.process.state() == QtCore.QProcess.NotRunning:
            if self._log_handle is not None:
                self._log_handle.close()
                self._log_handle = None
            self.run_button.setEnabled(True)
            self.stop_button.setEnabled(False)

    def _process_finished(self, exit_code, exit_status):
        self._read_process_output()
        if self._log_handle is not None:
            self._log_handle.close()
            self._log_handle = None
        normal_exit = exit_status == QtCore.QProcess.NormalExit
        succeeded = normal_exit and exit_code == 0
        self.log_output.appendPlainText(
            "\nMOSAICS finished with exit code {} ({}).".format(
                exit_code, "success" if succeeded else "failure"
            )
        )
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        if succeeded and self.auto_load.isChecked():
            self._load_output(show_missing=False)

    def _stop_process(self):
        if self.process.state() == QtCore.QProcess.NotRunning:
            return
        self.log_output.appendPlainText("\nStopping MOSAICS…")
        self.process.terminate()
        QtCore.QTimer.singleShot(3000, self._kill_if_running)

    def _kill_if_running(self):
        if self.process.state() != QtCore.QProcess.NotRunning:
            self.log_output.appendPlainText("MOSAICS did not stop; terminating it forcefully.")
            self.process.kill()

    def _load_output(self, checked=False, show_missing=True):
        selected = self.output_edit.text().strip()
        output = Path(selected).expanduser() if selected else None
        if output is None or not output.is_file():
            candidates = discover_outputs(self._project_directory())
            output = candidates[0] if candidates else None
        if output is None or not output.is_file():
            if show_missing:
                QtWidgets.QMessageBox.information(self, "PymoSAICS", "No output PDB was found.")
            return

        object_name = "pymosaics_{}".format(output.stem)
        object_name = "".join(character if character.isalnum() or character == "_" else "_" for character in object_name)
        try:
            cmd.load(str(output), object_name)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "PymoSAICS", "PyMOL could not load the output: {}".format(exc))
            return
        self.output_edit.setText(str(output))
        self.log_output.appendPlainText("Loaded {} as {}.".format(output, object_name))

    def closeEvent(self, event):
        if self.process.state() != QtCore.QProcess.NotRunning:
            answer = QtWidgets.QMessageBox.question(
                self,
                "PymoSAICS",
                "MOSAICS is still running. Stop it and close PymoSAICS?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No,
            )
            if answer != QtWidgets.QMessageBox.Yes:
                event.ignore()
                return
            self.process.kill()
            self.process.waitForFinished(3000)
        if self._log_handle is not None:
            self._log_handle.close()
            self._log_handle = None
        event.accept()
