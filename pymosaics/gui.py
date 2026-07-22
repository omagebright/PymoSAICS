"""Qt interface for building, running, visualizing, and analyzing MOSAICS jobs."""

import hashlib
import os
import re
from pathlib import Path
import tempfile
from typing import List, Optional, Sequence, Tuple

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from .core.analysis import (
    EnergySeries,
    discover_energy_series,
    discover_pdb_outputs,
    discover_project_files,
    latest_log,
    parse_acceptance_log,
    read_text_file,
)
from .core.builder import (
    InputSettings,
    default_settings,
    generate_mcmc_input,
    stage_force_field,
    write_mcmc_input,
)
from .core.catalog import (
    ANALYSIS_PRESETS,
    FORCEFIELD_ROOT,
    FORCE_FIELD_PROFILES,
    RUNTIME_PROFILES,
    analysis_preset,
    force_field_profile,
    make_bundled_executable_runnable,
    runtime_profile,
    runtime_supports_force_field,
)
from .core.config import ConfigError, ConfigStore
from .core.landscape import LandscapeResult, build_landscape, write_landscape_table
from .core.models import PreparedRun, RuntimeConfig
from .core.project import (
    PreparationError,
    format_diagnostics,
    planned_parameter_input,
    prepare_run,
    validate_project,
)
from .core.protein_prep import discover_pdb2pqr, prepare_protein_with_pdb2pqr
from .core.regions import RegionSettings, generate_region_file, write_region_file
from .core.runtime import has_errors, validate_runtime
from .core.structures import (
    DisulfideCandidate,
    PDBMetadata,
    detect_disulfides_text,
    fetch_rcsb_pdb,
    inspect_pdb_text,
    prepare_structure,
    unambiguous_disulfide_keys,
)
from .core.topology import validate_pdb_against_rtf


CHECKED = QtCore.Qt.Checked if hasattr(QtCore.Qt, "Checked") else QtCore.Qt.CheckState.Checked
UNCHECKED = QtCore.Qt.Unchecked if hasattr(QtCore.Qt, "Unchecked") else QtCore.Qt.CheckState.Unchecked
ITEM_USER_CHECKABLE = (
    QtCore.Qt.ItemIsUserCheckable
    if hasattr(QtCore.Qt, "ItemIsUserCheckable")
    else QtCore.Qt.ItemFlag.ItemIsUserCheckable
)
USER_ROLE = int(QtCore.Qt.UserRole if hasattr(QtCore.Qt, "UserRole") else QtCore.Qt.ItemDataRole.UserRole)
TOOLTIP_ROLE = int(
    QtCore.Qt.ToolTipRole
    if hasattr(QtCore.Qt, "ToolTipRole")
    else QtCore.Qt.ItemDataRole.ToolTipRole
)
SCROLLBAR_ALWAYS_OFF = (
    QtCore.Qt.ScrollBarAlwaysOff
    if hasattr(QtCore.Qt, "ScrollBarAlwaysOff")
    else QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff
)

THEME = {
    "shell": "#071820",
    "panel": "#0d252f",
    "card": "#102c36",
    "field": "#071b24",
    "line": "#2c5662",
    "text": "#f0f8f9",
    "muted": "#a9c0c6",
    "accent": "#36c6b4",
    "accent_dark": "#197e75",
}


def _enum_value(owner, direct_name: str, scoped_name: str, member_name: str):
    direct = getattr(owner, direct_name, None)
    if direct is not None:
        return direct
    return getattr(getattr(owner, scoped_name), member_name)


ALIGN_LEFT = _enum_value(QtCore.Qt, "AlignLeft", "AlignmentFlag", "AlignLeft")
ALIGN_TOP = _enum_value(QtCore.Qt, "AlignTop", "AlignmentFlag", "AlignTop")
ALIGN_VCENTER = _enum_value(QtCore.Qt, "AlignVCenter", "AlignmentFlag", "AlignVCenter")
FRAME_NO_FRAME = _enum_value(QtWidgets.QFrame, "NoFrame", "Shape", "NoFrame")
FORM_ALL_NON_FIXED_GROW = _enum_value(
    QtWidgets.QFormLayout,
    "AllNonFixedFieldsGrow",
    "FieldGrowthPolicy",
    "AllNonFixedFieldsGrow",
)
FORM_DONT_WRAP_ROWS = _enum_value(
    QtWidgets.QFormLayout,
    "DontWrapRows",
    "RowWrapPolicy",
    "DontWrapRows",
)
SIZE_POLICY_IGNORED = _enum_value(
    QtWidgets.QSizePolicy, "Ignored", "Policy", "Ignored"
)


PAINTER_ANTIALIASING = _enum_value(QtGui.QPainter, "Antialiasing", "RenderHint", "Antialiasing")
PLAIN_TEXT_NO_WRAP = _enum_value(QtWidgets.QPlainTextEdit, "NoWrap", "LineWrapMode", "NoWrap")
DIALOG_BUTTON_CLOSE = _enum_value(
    QtWidgets.QDialogButtonBox, "Close", "StandardButton", "Close"
)
DIALOG_BUTTON_ACTION = _enum_value(
    QtWidgets.QDialogButtonBox, "ActionRole", "ButtonRole", "ActionRole"
)
MESSAGE_YES = _enum_value(QtWidgets.QMessageBox, "Yes", "StandardButton", "Yes")
MESSAGE_NO = _enum_value(QtWidgets.QMessageBox, "No", "StandardButton", "No")
EXTENDED_SELECTION = _enum_value(
    QtWidgets.QAbstractItemView, "ExtendedSelection", "SelectionMode", "ExtendedSelection"
)
COMBO_ADJUST_MINIMUM = _enum_value(
    QtWidgets.QComboBox,
    "AdjustToMinimumContentsLengthWithIcon",
    "SizeAdjustPolicy",
    "AdjustToMinimumContentsLengthWithIcon",
)
DIALOG_ACCEPTED = _enum_value(QtWidgets.QDialog, "Accepted", "DialogCode", "Accepted")
PROCESS_NOT_RUNNING = _enum_value(QtCore.QProcess, "NotRunning", "ProcessState", "NotRunning")
PROCESS_NORMAL_EXIT = _enum_value(QtCore.QProcess, "NormalExit", "ExitStatus", "NormalExit")
TEXT_CURSOR_END = _enum_value(QtGui.QTextCursor, "End", "MoveOperation", "End")
PALETTE_ACTIVE = _enum_value(QtGui.QPalette, "Active", "ColorGroup", "Active")
PALETTE_INACTIVE = _enum_value(QtGui.QPalette, "Inactive", "ColorGroup", "Inactive")
PALETTE_DISABLED = _enum_value(QtGui.QPalette, "Disabled", "ColorGroup", "Disabled")


def _safe_name(value: str) -> str:
    return "".join(character if character.isalnum() or character == "_" else "_" for character in value)


def _palette_role(direct_name: str, member_name: str):
    return _enum_value(QtGui.QPalette, direct_name, "ColorRole", member_name)


def _install_combo_popup_theme(combo):
    """Give combo popups a host-independent palette and readable geometry."""

    view = QtWidgets.QListView(combo)
    view.setObjectName("pymosaicsComboPopup")
    view.setUniformItemSizes(True)
    combo.setView(view)

    palette = QtGui.QPalette(combo.palette())
    standard_colors = {
        "Window": THEME["card"],
        "WindowText": THEME["text"],
        "Base": THEME["card"],
        "AlternateBase": THEME["card"],
        "Text": THEME["text"],
        "Button": THEME["card"],
        "ButtonText": THEME["text"],
        "Highlight": THEME["accent_dark"],
        "HighlightedText": "#ffffff",
    }
    for group in (PALETTE_ACTIVE, PALETTE_INACTIVE):
        for role_name, color in standard_colors.items():
            palette.setColor(group, _palette_role(role_name, role_name), QtGui.QColor(color))
    disabled_colors = dict(standard_colors)
    disabled_colors.update(
        {
            "WindowText": "#789198",
            "Text": "#789198",
            "ButtonText": "#789198",
            "HighlightedText": "#c8d8db",
        }
    )
    for role_name, color in disabled_colors.items():
        palette.setColor(
            PALETTE_DISABLED,
            _palette_role(role_name, role_name),
            QtGui.QColor(color),
        )
    combo.setPalette(palette)
    for surface in (view, view.viewport(), view.window()):
        surface.setPalette(palette)
        surface.setAutoFillBackground(True)
    view.setStyleSheet(
        """
        QListView#pymosaicsComboPopup {
            background-color: #102c36;
            color: #edf7f8;
            border: 1px solid #3a7180;
            outline: 0;
            padding: 3px 0;
            selection-background-color: #197e75;
            selection-color: #ffffff;
        }
        QListView#pymosaicsComboPopup::item {
            background-color: #102c36;
            color: #edf7f8;
            min-height: 24px;
            padding: 3px 10px;
        }
        QListView#pymosaicsComboPopup::item:hover {
            background-color: #17414b;
            color: #ffffff;
        }
        QListView#pymosaicsComboPopup::item:selected {
            background-color: #197e75;
            color: #ffffff;
        }
        QListView#pymosaicsComboPopup::item:disabled {
            background-color: #0d252f;
            color: #789198;
        }
        """
    )
    metrics = combo.fontMetrics()
    measure = getattr(metrics, "horizontalAdvance", None)
    if measure is None:
        measure = metrics.width
    longest = max((measure(combo.itemText(index)) for index in range(combo.count())), default=0)
    view.setMinimumWidth(min(680, max(220, longest + 42)))


class EnergyPlot(QtWidgets.QWidget):
    """Small dependency-free line plot for MOSAICS energy output."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._series: Optional[EnergySeries] = None
        self.setMinimumHeight(210)

    def set_series(self, series: Optional[EnergySeries]):
        self._series = series
        self.update()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.fillRect(self.rect(), QtGui.QColor(THEME["field"]))
        if self._series is None or not self._series.values:
            painter.setPen(QtGui.QColor(THEME["muted"]))
            painter.drawText(20, 35, "No energy series selected")
            return

        values = self._series.values
        left, top, right, bottom = 64, 24, 18, 44
        width = max(1, self.width() - left - right)
        height = max(1, self.height() - top - bottom)
        low, high = min(values), max(values)
        span = high - low if high != low else 1.0
        axis_pen = QtGui.QPen(QtGui.QColor(THEME["line"]))
        painter.setPen(axis_pen)
        painter.drawLine(left, top, left, top + height)
        painter.drawLine(left, top + height, left + width, top + height)

        line_pen = QtGui.QPen(QtGui.QColor(THEME["accent"]))
        line_pen.setWidth(2)
        painter.setPen(line_pen)
        path = QtGui.QPainterPath()
        maximum_points = max(2, width * 2)
        if len(values) > maximum_points:
            indices = tuple(
                round(index * (len(values) - 1) / (maximum_points - 1))
                for index in range(maximum_points)
            )
        else:
            indices = tuple(range(len(values)))
        for position, index in enumerate(indices):
            value = values[index]
            x = left + (width * index / max(1, len(values) - 1))
            y = top + height - (height * (value - low) / span)
            if position == 0:
                path.moveTo(x, y)
            else:
                path.lineTo(x, y)
        painter.drawPath(path)

        painter.setPen(QtGui.QColor(THEME["muted"]))
        painter.drawText(4, top + 5, "{:.3g}".format(high))
        painter.drawText(4, top + height, "{:.3g}".format(low))
        painter.drawText(left, self.height() - 12, "sample 1")
        painter.drawText(left + width - 90, self.height() - 12, "sample {}".format(len(values)))


class LandscapePlot(QtWidgets.QWidget):
    """Clickable two-dimensional map of aligned structural RMSD distances."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._result: Optional[LandscapeResult] = None
        self._energies: Optional[Tuple[float, ...]] = None
        self._screen_points = []
        self._selected_frame = None
        self.on_frame_selected = None
        self.setMinimumHeight(230)

    def set_result(self, result: Optional[LandscapeResult], energies: Optional[Tuple[float, ...]] = None):
        self._result = result
        self._energies = energies
        self._selected_frame = None
        self.update()

    def _project(self):
        if self._result is None or not self._result.coordinates:
            return []
        left, top, right, bottom = 62.0, 28.0, 26.0, 48.0
        width = max(1.0, self.width() - left - right)
        height = max(1.0, self.height() - top - bottom)
        xs = [point[0] for point in self._result.coordinates]
        ys = [point[1] for point in self._result.coordinates]
        xlow, xhigh = min(xs), max(xs)
        ylow, yhigh = min(ys), max(ys)
        xspan = xhigh - xlow or 1.0
        yspan = yhigh - ylow or 1.0
        return [
            (
                left + width * (x - xlow) / xspan,
                top + height - height * (y - ylow) / yspan,
                frame,
            )
            for (x, y), frame in zip(self._result.coordinates, self._result.frame_numbers)
        ]

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(PAINTER_ANTIALIASING)
        painter.fillRect(self.rect(), QtGui.QColor(THEME["field"]))
        if self._result is None:
            painter.setPen(QtGui.QColor(THEME["muted"]))
            painter.drawText(20, 35, "Select a multi-model PDB trajectory and build the structural map")
            return
        left, top, right, bottom = 62, 28, 26, 48
        width = max(1, self.width() - left - right)
        height = max(1, self.height() - top - bottom)
        painter.setPen(QtGui.QPen(QtGui.QColor("#8aa0ad"), 1))
        painter.drawLine(left, top + height, left + width, top + height)
        painter.drawLine(left, top, left, top + height)
        painter.drawText(left + width // 2 - 45, self.height() - 13, "RMSD coordinate 1")
        painter.save()
        painter.translate(16, top + height // 2 + 45)
        painter.rotate(-90)
        painter.drawText(0, 0, "RMSD coordinate 2")
        painter.restore()

        representatives = set(self._result.representative_frames)
        self._screen_points = self._project()
        count = max(1, len(self._screen_points) - 1)
        energy_low = min(self._energies) if self._energies else 0.0
        energy_high = max(self._energies) if self._energies else 1.0
        energy_span = energy_high - energy_low or 1.0
        for index, (x, y, frame) in enumerate(self._screen_points):
            fraction = (
                (self._energies[index] - energy_low) / energy_span
                if self._energies
                else index / count
            )
            color = QtGui.QColor.fromRgbF(0.10 + 0.72 * fraction, 0.62 - 0.28 * fraction, 0.72 - 0.42 * fraction)
            radius = 7.5 if frame in representatives else 4.0
            if frame == self._selected_frame:
                painter.setPen(QtGui.QPen(QtGui.QColor("#f5b700"), 3))
                radius += 3
            else:
                painter.setPen(QtGui.QPen(QtGui.QColor(THEME["field"]), 1))
            painter.setBrush(color)
            painter.drawEllipse(QtCore.QPointF(x, y), radius, radius)
        painter.setPen(QtGui.QColor(THEME["muted"]))
        painter.drawText(left + 8, top + 16, "Low energy" if self._energies else "Early frames")
        painter.drawText(left + width - 88, top + 16, "High energy" if self._energies else "Late frames")

    def mousePressEvent(self, event):
        if not self._screen_points:
            return
        position = event.position() if hasattr(event, "position") else event.pos()
        closest = min(
            self._screen_points,
            key=lambda point: (point[0] - position.x()) ** 2 + (point[1] - position.y()) ** 2,
        )
        distance_squared = (closest[0] - position.x()) ** 2 + (closest[1] - position.y()) ** 2
        if distance_squared <= 16 ** 2:
            self._selected_frame = closest[2]
            self.update()
            if self.on_frame_selected:
                self.on_frame_selected(closest[2])


class TextFileDialog(QtWidgets.QDialog):
    """View a project text file and optionally save explicit user edits."""

    def __init__(self, path: Path, editable: bool, parent=None):
        super().__init__(parent)
        self.path = path.expanduser().resolve()
        self.editable = editable
        self.setWindowTitle(("Edit " if editable else "View ") + self.path.name)
        self.resize(900, 680)
        layout = QtWidgets.QVBoxLayout(self)
        location = QtWidgets.QLineEdit(str(self.path))
        location.setReadOnly(True)
        layout.addWidget(location)
        self.editor = QtWidgets.QPlainTextEdit()
        self.editor.setLineWrapMode(PLAIN_TEXT_NO_WRAP)
        self.editor.setReadOnly(not editable)
        self.editor.document().setModified(False)
        layout.addWidget(self.editor, 1)
        self.status = QtWidgets.QLabel()
        layout.addWidget(self.status)
        buttons = QtWidgets.QDialogButtonBox(DIALOG_BUTTON_CLOSE)
        reload_button = buttons.addButton("Reload", DIALOG_BUTTON_ACTION)
        reload_button.clicked.connect(self.reload)
        if editable:
            save_button = buttons.addButton("Save changes", DIALOG_BUTTON_ACTION)
            save_button.clicked.connect(self.save)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.reload()

    def reload(self, checked=False):
        try:
            content = read_text_file(self.path)
        except (OSError, ValueError) as exc:
            self.editor.clear()
            self.status.setText(str(exc))
            return
        self.editor.setPlainText(content)
        self.editor.document().setModified(False)
        self.status.setText("{} bytes · UTF-8 text view".format(self.path.stat().st_size))

    def save(self, checked=False):
        if not self.editable:
            return
        temporary = None
        try:
            original_mode = self.path.stat().st_mode
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=str(self.path.parent),
                prefix=".pymosaics-edit-",
                suffix=".tmp",
                delete=False,
            ) as handle:
                handle.write(self.editor.toPlainText())
                temporary = Path(handle.name)
            os.replace(str(temporary), str(self.path))
            self.path.chmod(original_mode & 0o7777)
        except OSError as exc:
            QtWidgets.QMessageBox.critical(self, "PymoSAICS", "Cannot save {}: {}".format(self.path.name, exc))
            return
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()
        self.editor.document().setModified(False)
        self.status.setText("Saved {}".format(self.path))

    def reject(self):
        if self.editable and self.editor.document().isModified():
            answer = QtWidgets.QMessageBox.question(
                self,
                "Unsaved changes",
                "Close without saving changes to {}?".format(self.path.name),
                MESSAGE_YES | MESSAGE_NO,
                MESSAGE_NO,
            )
            if answer != MESSAGE_YES:
                return
        super().reject()


class RegionEditor(QtWidgets.QDialog):
    """Scientific workbench for one explicit residue-level MOSAICS region."""

    PRESETS = {
        "WP2 balanced pilot": {
            "translation_sigma": "0.0",
            "rotation_sigma": "0.0",
            "free_translation_sigma": ".5e-5",
            "free_rotation_sigma": ".5e-6",
            "pair_translation_sigma": ".5e-5",
            "pair_rotation_sigma": ".5e-6",
        },
        "WP2 paired-residue motion": {
            "translation_sigma": "0.0",
            "rotation_sigma": "0.0",
            "free_translation_sigma": "0.0",
            "free_rotation_sigma": "0.0",
            "pair_translation_sigma": ".5e-5",
            "pair_rotation_sigma": ".5e-6",
        },
    }

    def __init__(self, residues: Sequence[str], pymol_object: str, initial=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MOSAICS region workbench")
        self.resize(1080, 790)
        self.setMinimumSize(920, 690)
        self.pymol_object = pymol_object
        self.settings: Optional[RegionSettings] = None
        self._changing_table = False
        self._changing_sigma = False

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(18, 16, 18, 16)
        outer.setSpacing(12)
        explanation = QtWidgets.QLabel(
            "Define one residue-level natural-move region. Choose its members, at least one "
            "rotation center, optional non-overlapping residue pairs, and explicit proposal widths."
        )
        explanation.setWordWrap(True)
        outer.addWidget(explanation)
        self.selection_summary = QtWidgets.QLabel()
        self.selection_summary.setObjectName("contextBadge")
        outer.addWidget(self.selection_summary)

        workspace = QtWidgets.QSplitter()
        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 4, 0)
        left_layout.setSpacing(10)

        selection_group = QtWidgets.QGroupBox("1 · Region membership and rotation centers")
        selection_layout = QtWidgets.QVBoxLayout(selection_group)
        selection_note = QtWidgets.QLabel(
            "Move includes a residue. Center marks a residue MOSAICS may use as the pivot "
            "for whole-region rotation; at least one center is required."
        )
        selection_note.setWordWrap(True)
        selection_layout.addWidget(selection_note)
        self.residue_table = QtWidgets.QTableWidget(len(residues), 3)
        self.residue_table.setHorizontalHeaderLabels(("Move", "Center", "Residue (chain:number)"))
        self.residue_table.horizontalHeader().setStretchLastSection(True)
        self.residue_table.setColumnWidth(0, 72)
        self.residue_table.setColumnWidth(1, 76)
        initial_residues = set(initial.residues if initial else residues)
        initial_centers = set(initial.centers if initial else residues[:1])
        for row, selector in enumerate(residues):
            include = QtWidgets.QTableWidgetItem()
            include.setFlags(include.flags() | ITEM_USER_CHECKABLE)
            include.setCheckState(CHECKED if selector in initial_residues else UNCHECKED)
            center = QtWidgets.QTableWidgetItem()
            center.setFlags(center.flags() | ITEM_USER_CHECKABLE)
            center.setCheckState(CHECKED if selector in initial_centers else UNCHECKED)
            self.residue_table.setItem(row, 0, include)
            self.residue_table.setItem(row, 1, center)
            self.residue_table.setItem(row, 2, QtWidgets.QTableWidgetItem(selector))
        self.residue_table.itemChanged.connect(self._table_changed)
        selection_layout.addWidget(self.residue_table, 1)
        selection_actions = QtWidgets.QHBoxLayout()
        select_all = QtWidgets.QPushButton("Select all")
        select_all.clicked.connect(lambda _checked=False: self._set_all_residues(True))
        clear_all = QtWidgets.QPushButton("Clear")
        clear_all.clicked.connect(lambda _checked=False: self._set_all_residues(False))
        from_pymol = QtWidgets.QPushButton("Use current PyMOL selection")
        from_pymol.setToolTip("Replace membership with residues in PyMOL's current 'sele' selection.")
        from_pymol.clicked.connect(self._use_pymol_selection)
        selection_actions.addWidget(select_all)
        selection_actions.addWidget(clear_all)
        selection_actions.addStretch(1)
        selection_actions.addWidget(from_pymol)
        selection_layout.addLayout(selection_actions)
        left_layout.addWidget(selection_group, 3)

        pair_group = QtWidgets.QGroupBox("2 · Residue pairs")
        pair_layout = QtWidgets.QVBoxLayout(pair_group)
        pair_note = QtWidgets.QLabel(
            "Pair two selected residues that should move as one coupled unit, such as a "
            "nucleic-acid base pair. A residue can belong to only one pair."
        )
        pair_note.setWordWrap(True)
        pair_layout.addWidget(pair_note)
        pair_row = QtWidgets.QGridLayout()
        self.pair_first = QtWidgets.QComboBox()
        self.pair_second = QtWidgets.QComboBox()
        add_pair = QtWidgets.QPushButton("Add pair")
        add_pair.clicked.connect(self._add_pair)
        remove_pair = QtWidgets.QPushButton("Remove pair")
        remove_pair.clicked.connect(self._remove_pair)
        pair_row.addWidget(QtWidgets.QLabel("First"), 0, 0)
        pair_row.addWidget(QtWidgets.QLabel("Second"), 0, 1)
        pair_row.addWidget(self.pair_first, 1, 0)
        pair_row.addWidget(self.pair_second, 1, 1)
        pair_row.addWidget(add_pair, 1, 2)
        pair_row.addWidget(remove_pair, 1, 3)
        pair_layout.addLayout(pair_row)
        self.pair_list = QtWidgets.QListWidget()
        self.pair_list.setFixedHeight(88)
        if initial:
            for first, second in initial.residue_pairs:
                self.pair_list.addItem("{} — {}".format(first, second))
        pair_layout.addWidget(self.pair_list)
        left_layout.addWidget(pair_group)

        right = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(right)
        right_layout.setContentsMargins(4, 0, 0, 0)
        right_layout.setSpacing(10)
        model_group = QtWidgets.QGroupBox("3 · Scientific model")
        model_form = QtWidgets.QFormLayout(model_group)
        model_form.setHorizontalSpacing(12)
        self.element_type = QtWidgets.QLineEdit("Residue (supported)")
        self.element_type.setReadOnly(True)
        self.element_type.setToolTip("Writes \\element_top_type{residue}.")
        model_form.addRow("Element type:", self.element_type)
        self.dependency = QtWidgets.QLineEdit("Independent (supported)")
        self.dependency.setReadOnly(True)
        self.dependency.setToolTip(
            "Writes \\dependency_type{independent}; dependent hierarchy is not valid in this editor."
        )
        model_form.addRow("Dependency:", self.dependency)
        self.propagation = QtWidgets.QLineEdit("Superimpose (single region)")
        self.propagation.setReadOnly(True)
        self.propagation.setToolTip(
            "Writes \\prop_regions_type{superimpose} in mcmc.input. With one region, "
            "superimpose and onebyone select the same sole region."
        )
        model_form.addRow("Propagation:", self.propagation)
        self.preset = QtWidgets.QComboBox()
        self.preset.addItems(tuple(self.PRESETS) + ("Custom",))
        self.preset.setToolTip(
            "WP2 presets reproduce the documented MOSAICS tutorial widths. Run a short pilot "
            "and tune proposal widths from measured acceptance."
        )
        self.preset.currentTextChanged.connect(self._apply_preset)
        model_form.addRow("Proposal preset:", self.preset)
        right_layout.addWidget(model_group)

        motion_group = QtWidgets.QGroupBox("4 · Natural-move proposal widths")
        motion_layout = QtWidgets.QGridLayout(motion_group)
        motion_layout.setHorizontalSpacing(10)
        motion_layout.setVerticalSpacing(7)
        motion_layout.addWidget(QtWidgets.QLabel("Motion level"), 0, 0)
        motion_layout.addWidget(QtWidgets.QLabel("Translation σ (Å)"), 0, 1)
        motion_layout.addWidget(QtWidgets.QLabel("Rotation σ (rad)"), 0, 2)
        self.sigma_edits = {}
        motion_rows = (
            ("Whole region", "translation_sigma", "rotation_sigma",
             "Moves all included residues together about a selected center."),
            ("Free residues", "free_translation_sigma", "free_rotation_sigma",
             "Moves individual included residues not assigned to a pair."),
            ("Residue pairs", "pair_translation_sigma", "pair_rotation_sigma",
             "Moves each declared pair as a coupled unit."),
        )
        defaults = self.PRESETS["WP2 balanced pilot"]
        validator = QtGui.QDoubleValidator(0.0, 1.0e9, 12, self)
        compact_definitions = {
            "Whole region": "all members together",
            "Free residues": "unpaired members individually",
            "Residue pairs": "each pair as one unit",
        }
        for row, (label, trans_name, rot_name, definition) in enumerate(motion_rows, start=1):
            motion_layout.setRowMinimumHeight(row, 40)
            level = QtWidgets.QWidget()
            level_layout = QtWidgets.QVBoxLayout(level)
            level_layout.setContentsMargins(0, 0, 0, 0)
            level_layout.setSpacing(0)
            level_title = QtWidgets.QLabel(label)
            level_font = level_title.font()
            level_font.setBold(True)
            level_title.setFont(level_font)
            level_definition = QtWidgets.QLabel(compact_definitions[label])
            level_definition.setStyleSheet("color: #a9c0c6;")
            level_layout.addWidget(level_title)
            level_layout.addWidget(level_definition)
            level.setToolTip(definition)
            motion_layout.addWidget(level, row, 0)
            for column, attribute in ((1, trans_name), (2, rot_name)):
                edit = QtWidgets.QLineEdit(getattr(initial, attribute) if initial else defaults[attribute])
                edit.setValidator(validator)
                edit.setToolTip(definition)
                edit.textEdited.connect(self._sigma_edited)
                edit.textChanged.connect(self._update_preview)
                self.sigma_edits[attribute] = edit
                motion_layout.addWidget(edit, row, column)
        motion_layout.setColumnStretch(0, 3)
        motion_layout.setColumnStretch(1, 1)
        motion_layout.setColumnStretch(2, 1)
        right_layout.addWidget(motion_group)

        preview_group = QtWidgets.QGroupBox("5 · Generated region.data")
        preview_layout = QtWidgets.QVBoxLayout(preview_group)
        self.preview = QtWidgets.QPlainTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setMinimumHeight(85)
        preview_layout.addWidget(self.preview, 1)
        self.validation_status = QtWidgets.QLabel()
        self.validation_status.setWordWrap(True)
        self.validation_status.setMaximumHeight(38)
        preview_layout.addWidget(self.validation_status)
        right_layout.addWidget(preview_group, 1)

        workspace.addWidget(left)
        workspace.addWidget(right)
        workspace.setStretchFactor(0, 1)
        workspace.setStretchFactor(1, 1)
        outer.addWidget(workspace, 1)

        controls = QtWidgets.QHBoxLayout()
        show_button = QtWidgets.QPushButton("Show region in PyMOL")
        show_button.clicked.connect(self._show_in_pymol)
        controls.addWidget(show_button)
        documentation = QtWidgets.QLabel(
            '<a style="color:#75dfd1" href="https://www.cs.ox.ac.uk/mosaics/Documentation.php">'
            "MOSAICS region documentation</a>"
        )
        documentation.setOpenExternalLinks(True)
        controls.addWidget(documentation)
        controls.addStretch(1)
        cancel = QtWidgets.QPushButton("Cancel")
        self.save_button = QtWidgets.QPushButton("Use this region")
        self.save_button.setObjectName("primaryAction")
        cancel.clicked.connect(self.reject)
        self.save_button.clicked.connect(self._accept_settings)
        controls.addWidget(cancel)
        controls.addWidget(self.save_button)
        outer.addLayout(controls)
        for combo in self.findChildren(QtWidgets.QComboBox):
            _install_combo_popup_theme(combo)
        self._sync_pair_choices()
        self._select_matching_preset()
        self._update_preview()

    def _selected(self, column: int) -> Tuple[str, ...]:
        values = []
        for row in range(self.residue_table.rowCount()):
            if self.residue_table.item(row, column).checkState() == CHECKED:
                values.append(self.residue_table.item(row, 2).text())
        return tuple(values)

    def _pairs(self) -> Tuple[Tuple[str, str], ...]:
        pairs = []
        for row in range(self.pair_list.count()):
            first, second = self.pair_list.item(row).text().split(" — ", 1)
            pairs.append((first, second))
        return tuple(pairs)

    def _current_settings(self) -> RegionSettings:
        return RegionSettings(
            residues=self._selected(0),
            centers=self._selected(1),
            residue_pairs=self._pairs(),
            dependency_type="independent",
            **{name: edit.text().strip() for name, edit in self.sigma_edits.items()}
        )

    def _update_preview(self, *_args):
        selected, centers, pairs = self._selected(0), self._selected(1), self._pairs()
        self.selection_summary.setText(
            "{} residues  ·  {} rotation center{}  ·  {} residue pair{}".format(
                len(selected), len(centers), "" if len(centers) == 1 else "s",
                len(pairs), "" if len(pairs) == 1 else "s",
            )
        )
        try:
            generated = generate_region_file(self._current_settings())
        except ValueError as exc:
            self.preview.setPlainText("Cannot generate region: {}".format(exc))
            self.validation_status.setText("⚠ {}".format(exc))
            self.validation_status.setStyleSheet("color: #f2ad72;")
            self.save_button.setEnabled(False)
            return
        self.preview.setPlainText(generated)
        note = "Valid. Review region.data before running."
        if not pairs and any(
            float(self.sigma_edits[name].text() or 0) > 0
            for name in ("pair_translation_sigma", "pair_rotation_sigma")
        ):
            note += " Pair widths are inactive: no pair is declared."
        self.validation_status.setText("✓ " + note)
        self.validation_status.setStyleSheet("color: #75dfd1;")
        self.save_button.setEnabled(True)

    def _table_changed(self, item):
        if self._changing_table:
            return
        self._changing_table = True
        try:
            row = item.row()
            if item.column() == 1 and item.checkState() == CHECKED:
                self.residue_table.item(row, 0).setCheckState(CHECKED)
            elif item.column() == 0 and item.checkState() != CHECKED:
                self.residue_table.item(row, 1).setCheckState(UNCHECKED)
                removed = self.residue_table.item(row, 2).text()
                for pair_row in reversed(range(self.pair_list.count())):
                    if removed in self.pair_list.item(pair_row).text().split(" — "):
                        self.pair_list.takeItem(pair_row)
        finally:
            self._changing_table = False
        self._sync_pair_choices()
        self._update_preview()

    def _set_all_residues(self, selected: bool):
        self._changing_table = True
        try:
            for row in range(self.residue_table.rowCount()):
                self.residue_table.item(row, 0).setCheckState(CHECKED if selected else UNCHECKED)
                if not selected:
                    self.residue_table.item(row, 1).setCheckState(UNCHECKED)
            if not selected:
                self.pair_list.clear()
        finally:
            self._changing_table = False
        self._sync_pair_choices()
        self._update_preview()

    def _sync_pair_choices(self):
        selected = self._selected(0)
        previous = (self.pair_first.currentText(), self.pair_second.currentText())
        for combo, old_value in zip((self.pair_first, self.pair_second), previous):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(selected)
            if old_value in selected:
                combo.setCurrentText(old_value)
            combo.blockSignals(False)

    def _sigma_edited(self, _text):
        if self._changing_sigma:
            return
        self.preset.blockSignals(True)
        self.preset.setCurrentText("Custom")
        self.preset.blockSignals(False)

    def _select_matching_preset(self):
        values = {name: edit.text().strip() for name, edit in self.sigma_edits.items()}
        name = next((name for name, preset in self.PRESETS.items() if values == preset), "Custom")
        self.preset.setCurrentText(name)

    def _apply_preset(self, name):
        if name not in self.PRESETS:
            return
        self._changing_sigma = True
        try:
            for attribute, value in self.PRESETS[name].items():
                self.sigma_edits[attribute].setText(value)
        finally:
            self._changing_sigma = False
        self._update_preview()

    def _use_pymol_selection(self, checked=False):
        try:
            model = cmd.get_model("({}) and sele".format(self.pymol_object))
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Cannot read the PyMOL selection: {}".format(exc))
            return
        selected = {"{}:{}".format(atom.chain or "_", atom.resi) for atom in model.atom}
        available = {
            self.residue_table.item(row, 2).text()
            for row in range(self.residue_table.rowCount())
        }
        selected &= available
        if not selected:
            QtWidgets.QMessageBox.information(
                self, "PymoSAICS",
                "PyMOL's current selection contains no residues from the active object and chains.",
            )
            return
        self._changing_table = True
        try:
            first_selected_row = None
            for row in range(self.residue_table.rowCount()):
                selector = self.residue_table.item(row, 2).text()
                include = selector in selected
                self.residue_table.item(row, 0).setCheckState(CHECKED if include else UNCHECKED)
                if include and first_selected_row is None:
                    first_selected_row = row
                if not include:
                    self.residue_table.item(row, 1).setCheckState(UNCHECKED)
            self.pair_list.clear()
            if not self._selected(1) and first_selected_row is not None:
                self.residue_table.item(first_selected_row, 1).setCheckState(CHECKED)
        finally:
            self._changing_table = False
        self._sync_pair_choices()
        self._update_preview()

    def _add_pair(self, checked=False):
        first, second = self.pair_first.currentText(), self.pair_second.currentText()
        if not first or not second:
            return
        if first == second:
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "Choose two different residues.")
            return
        used = {item for pair in self._pairs() for item in pair}
        if first in used or second in used:
            QtWidgets.QMessageBox.information(
                self, "PymoSAICS",
                "Each residue can belong to only one pair. Remove the existing pair first.",
            )
            return
        self.pair_list.addItem("{} — {}".format(first, second))
        self._update_preview()

    def _remove_pair(self, checked=False):
        for item in self.pair_list.selectedItems():
            self.pair_list.takeItem(self.pair_list.row(item))
        self._update_preview()

    def _selector_expression(self, selectors: Sequence[str]) -> str:
        parts = []
        for selector in selectors:
            chain, residue = selector.split(":", 1)
            chain_expression = "chain {}".format(chain) if chain != "_" else "chain ''"
            parts.append("({} and resi {})".format(chain_expression, residue))
        return " or ".join(parts) or "none"

    def _show_in_pymol(self, checked=False):
        try:
            selected = self._selected(0)
            centers = self._selected(1)
            paired = tuple(item for pair in self._pairs() for item in pair)
            cmd.select(
                "pymosaics_region",
                "({}) and ({})".format(self.pymol_object, self._selector_expression(selected)),
            )
            cmd.show("sticks", "pymosaics_region")
            cmd.color("cyan", "pymosaics_region")
            cmd.select(
                "pymosaics_region_pairs",
                "({}) and ({})".format(self.pymol_object, self._selector_expression(paired)),
            )
            cmd.show("sticks", "pymosaics_region_pairs")
            cmd.color("magenta", "pymosaics_region_pairs")
            cmd.select(
                "pymosaics_region_centers",
                "({}) and ({})".format(self.pymol_object, self._selector_expression(centers)),
            )
            cmd.show("spheres", "pymosaics_region_centers and name CA+C1*")
            cmd.color("yellow", "pymosaics_region_centers")
            cmd.zoom("pymosaics_region")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Cannot show region: {}".format(exc))

    def _accept_settings(self, checked=False):
        try:
            self.settings = self._current_settings()
            generate_region_file(self.settings)
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", str(exc))
            return
        self.accept()


class PymoSAICSDialog(QtWidgets.QDialog):
    """Build, run, visualize, and analyze a reproducible MOSAICS project."""

    def __init__(self, parent=None, config_store: Optional[ConfigStore] = None):
        super().__init__(parent)
        self.setWindowTitle("PymoSAICS")
        self.resize(1120, 860)
        self.setMinimumSize(900, 710)
        self.config_store = config_store or ConfigStore()
        self.process = QtCore.QProcess(self)
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read_process_output)
        self.process.finished.connect(self._process_finished)
        self.process.errorOccurred.connect(self._process_error)
        self._active_run: Optional[PreparedRun] = None
        self._log_handle = None
        self._source_path: Optional[Path] = None
        self._source_text = ""
        self._metadata: Optional[PDBMetadata] = None
        self._region_settings: Optional[RegionSettings] = None
        self._live_digest = ""
        self._live_object_name = ""
        self._loaded_input_structure: Optional[Path] = None
        self._energy_series: Tuple[EnergySeries, ...] = ()
        self._landscape_result: Optional[LandscapeResult] = None
        self._landscape_path: Optional[Path] = None
        self._landscape_object = ""
        self._preview_signature = None

        self.tabs = QtWidgets.QTabWidget()
        self.build_tab = QtWidgets.QWidget()
        self.build_tab.setObjectName("buildTab")
        self.run_tab = QtWidgets.QWidget()
        self.run_tab.setObjectName("runTab")
        self.analysis_tab = QtWidgets.QWidget()
        self.analysis_tab.setObjectName("analysisTab")
        self.setup_tab = QtWidgets.QWidget()
        self.setup_tab.setObjectName("setupTab")
        self.about_tab = QtWidgets.QWidget()
        self.about_tab.setObjectName("aboutTab")
        self.tabs.addTab(self.build_tab, "Build")
        self.tabs.addTab(self.run_tab, "Run")
        self.tabs.addTab(self.analysis_tab, "Analysis")
        self.tabs.addTab(self.setup_tab, "Setup")
        self.tabs.addTab(self.about_tab, "About")
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(18, 16, 18, 18)
        outer.setSpacing(12)
        outer.addWidget(self._build_product_header())
        outer.addWidget(self.tabs)

        self._build_build_tab()
        self._build_run_tab()
        self._build_analysis_tab()
        self._build_setup_tab()
        self._build_about_tab()
        self._apply_style()
        self._load_configuration()
        self._refresh_pymol_objects()

        self.live_timer = QtCore.QTimer(self)
        self.live_timer.timeout.connect(self._poll_pymol)
        self.live_timer.start(1000)

    def _build_product_header(self):
        header = QtWidgets.QFrame()
        header.setObjectName("productHeader")
        layout = QtWidgets.QHBoxLayout(header)
        layout.setContentsMargins(0, 0, 18, 0)
        layout.setSpacing(16)
        identity_rail = QtWidgets.QFrame()
        identity_rail.setObjectName("identityRail")
        identity_rail.setFixedWidth(5)
        layout.addWidget(identity_rail)
        title_column = QtWidgets.QVBoxLayout()
        title_column.setContentsMargins(0, 13, 0, 13)
        title_column.setSpacing(3)
        title = QtWidgets.QLabel("PymoSAICS")
        title.setObjectName("productTitle")
        subtitle = QtWidgets.QLabel("Prepare · run · inspect MOSAICS without hidden inputs")
        subtitle.setObjectName("productSubtitle")
        title_column.addWidget(title)
        title_column.addWidget(subtitle)
        layout.addLayout(title_column)
        layout.addStretch(1)
        context = QtWidgets.QVBoxLayout()
        context.setContentsMargins(0, 11, 0, 11)
        context.setSpacing(6)
        self.header_runtime = QtWidgets.QLabel("Runtime · not selected")
        self.header_runtime.setObjectName("contextBadge")
        self.header_science = QtWidgets.QLabel("Force field · not selected")
        self.header_science.setObjectName("contextBadge")
        context.addWidget(self.header_runtime)
        context.addWidget(self.header_science)
        layout.addLayout(context)
        return header

    @staticmethod
    def _palette_role(direct_name: str, member_name: str):
        return _palette_role(direct_name, member_name)

    def _apply_theme_palette(self):
        """Keep native Qt themes from replacing plugin surfaces with light colors."""

        palette = QtGui.QPalette(self.palette())
        colors = {
            "Window": THEME["shell"],
            "WindowText": THEME["text"],
            "Base": THEME["field"],
            "AlternateBase": THEME["card"],
            "Text": THEME["text"],
            "Button": THEME["card"],
            "ButtonText": THEME["text"],
            "Highlight": THEME["accent_dark"],
            "HighlightedText": "#ffffff",
            "ToolTipBase": THEME["card"],
            "ToolTipText": THEME["text"],
        }
        for role_name, color in colors.items():
            palette.setColor(self._palette_role(role_name, role_name), QtGui.QColor(color))
        self.setPalette(palette)

    def _apply_style(self):
        self._apply_theme_palette()
        self.setStyleSheet(
            """
            QDialog { background: #071820; color: #f0f8f9; font-family: "Aptos", "Segoe UI", "Helvetica Neue", sans-serif; font-size: 13px; }
            QWidget#buildTab, QWidget#runTab, QWidget#analysisTab, QWidget#setupTab, QWidget#aboutTab { background: #0d252f; }
            QLabel, QCheckBox { color: #dcebed; background: transparent; }
            QLabel:disabled, QCheckBox:disabled { color: #6f8991; }
            QFrame#productHeader { background: #0d2934; border: 1px solid #2c5662; border-radius: 10px; }
            QFrame#identityRail { background: #36c6b4; border: 0; border-top-left-radius: 9px; border-bottom-left-radius: 9px; }
            QLabel#productTitle { color: #f5fbfc; font-size: 24px; font-weight: 700; letter-spacing: 1px; }
            QLabel#productSubtitle { color: #a9c0c6; font-size: 12px; }
            QLabel#contextBadge { color: #d8fbf5; background: #123944; border: 1px solid #34707a; border-radius: 6px; padding: 4px 9px; }
            QLabel#sectionLabel { color: #a9c0c6; font-size: 11px; font-weight: 700; letter-spacing: 1px; }
            QTabWidget::pane { border: 1px solid #2c5662; border-radius: 8px; background: #0d252f; top: -1px; }
            QTabBar::tab { background: #0b2029; color: #9eb6bc; border: 1px solid #294d58; padding: 9px 18px; margin-right: 3px; border-top-left-radius: 6px; border-top-right-radius: 6px; }
            QTabBar::tab:selected { background: #16414b; color: #f4fffd; border-color: #39717b; border-bottom-color: #16414b; }
            QTabBar::tab:hover { color: #ffffff; background: #123640; }
            QTabWidget#analysisPages QTabBar::tab { min-width: 142px; padding-left: 14px; padding-right: 14px; }
            QGroupBox { background: #102c36; color: #e7f1f3; font-weight: 600; border: 1px solid #315963; border-radius: 8px; margin-top: 12px; padding: 14px 11px 11px 11px; }
            QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 7px; color: #75dfd1; background: #102c36; }
            QLineEdit, QPlainTextEdit, QComboBox, QListWidget, QTableWidget, QSpinBox, QDoubleSpinBox { background: #071b24; color: #edf7f8; border: 1px solid #37616d; border-radius: 5px; selection-background-color: #197e75; selection-color: #ffffff; padding: 5px 7px; }
            QLineEdit:hover, QPlainTextEdit:hover, QComboBox:hover, QListWidget:hover, QTableWidget:hover, QSpinBox:hover, QDoubleSpinBox:hover { border-color: #4a7f8b; }
            QLineEdit:focus, QPlainTextEdit:focus, QComboBox:focus, QListWidget:focus, QTableWidget:focus, QSpinBox:focus, QDoubleSpinBox:focus { border: 1px solid #36c6b4; }
            QLineEdit:read-only, QPlainTextEdit:read-only { background: #091f28; color: #c9dadd; }
            QPlainTextEdit { font-family: Menlo, "Cascadia Code", Consolas, monospace; font-size: 12px; }
            QComboBox::drop-down { border: 0; width: 25px; }
            QComboBox QAbstractItemView { background: #102c36; color: #edf7f8; border: 1px solid #3a7180; outline: 0; selection-background-color: #197e75; }
            QHeaderView::section { background: #153943; color: #d6e7ea; border: 0; border-right: 1px solid #315963; border-bottom: 1px solid #315963; padding: 7px; }
            QPushButton { background: #173b47; color: #edf7f8; border: 1px solid #3b6975; border-radius: 6px; min-height: 18px; padding: 6px 12px; font-weight: 600; }
            QPushButton:hover { background: #205160; border-color: #548894; }
            QPushButton:focus { border: 1px solid #36c6b4; }
            QPushButton:pressed { background: #0d2a34; }
            QPushButton:disabled { background: #14272e; color: #688087; border-color: #28414a; }
            QPushButton#primaryAction { background: #197e75; border-color: #43cdbc; color: #ffffff; }
            QPushButton#primaryAction:hover { background: #22978b; }
            QPushButton#stopAction { background: #67323a; border-color: #a95862; }
            QPushButton#stopAction:disabled { background: #2a2226; color: #765d62; border-color: #49383d; }
            QCheckBox { spacing: 7px; }
            QSplitter::handle { background: #315963; width: 8px; }
            QScrollArea#buildScroll, QScrollArea#setupScroll { background: #0d252f; border: 0; }
            QWidget#buildScrollContent, QWidget#setupScrollContent { background: #0d252f; }
            QScrollBar:vertical { background: #091c24; width: 11px; margin: 0; }
            QScrollBar::handle:vertical { background: #356471; border-radius: 5px; min-height: 34px; }
            QScrollBar::handle:vertical:hover { background: #477d89; }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
            QScrollBar:horizontal { background: #091c24; height: 11px; margin: 0; }
            QScrollBar::handle:horizontal { background: #356471; border-radius: 5px; min-width: 34px; }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
            QToolTip { background: #153943; color: #f0f8f9; border: 1px solid #4b7e89; padding: 4px; }
            """
        )
        for combo in self.findChildren(QtWidgets.QComboBox):
            _install_combo_popup_theme(combo)

    @staticmethod
    def _configure_form(form):
        form.setFieldGrowthPolicy(FORM_ALL_NON_FIXED_GROW)
        form.setRowWrapPolicy(FORM_DONT_WRAP_ROWS)
        form.setLabelAlignment(ALIGN_LEFT | ALIGN_VCENTER)
        form.setFormAlignment(ALIGN_LEFT | ALIGN_TOP)
        form.setHorizontalSpacing(14)
        form.setVerticalSpacing(9)

    @staticmethod
    def _allow_label_to_wrap(label):
        policy = label.sizePolicy()
        policy.setHorizontalPolicy(SIZE_POLICY_IGNORED)
        label.setSizePolicy(policy)
        label.setMinimumWidth(0)

    @staticmethod
    def _allow_text_area_to_shrink(editor, minimum_height=56):
        policy = editor.sizePolicy()
        policy.setVerticalPolicy(SIZE_POLICY_IGNORED)
        editor.setSizePolicy(policy)
        editor.setMinimumHeight(minimum_height)

    def _path_row(self, line_edit, button_text, callback):
        container = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(line_edit, 1)
        button = QtWidgets.QPushButton(button_text)
        button.clicked.connect(callback)
        layout.addWidget(button)
        return container

    @staticmethod
    def _populate_force_field_combo(combo):
        combo.setSizeAdjustPolicy(COMBO_ADJUST_MINIMUM)
        combo.setMinimumContentsLength(12)
        combo.setMaxVisibleItems(12)
        for profile in FORCE_FIELD_PROFILES:
            combo.addItem(profile.label, profile.identifier)
            combo.setItemData(combo.count() - 1, profile.description, TOOLTIP_ROLE)

    @staticmethod
    def _force_field_bundle_summary(profile):
        entries = (
            ("mol_parm_file", profile.rtf),
            ("bond_database_file", profile.bond),
            ("bend_database_file", profile.bend),
            ("tors_database_file", profile.torsion),
            ("onfo_database_file", profile.one_four),
            ("inter_database_file", profile.nonbonded),
        )
        return "\n".join(
            "\\{}{{forcefield/{}}}".format(directive, Path(value).name)
            for directive, value in entries
        )

    def _build_build_tab(self):
        outer = QtWidgets.QVBoxLayout(self.build_tab)
        outer.setContentsMargins(12, 12, 12, 12)
        splitter = QtWidgets.QSplitter()
        outer.addWidget(splitter)
        left = QtWidgets.QScrollArea()
        left.setObjectName("buildScroll")
        left.setFrameShape(FRAME_NO_FRAME)
        left.setWidgetResizable(True)
        left.setHorizontalScrollBarPolicy(SCROLLBAR_ALWAYS_OFF)
        left_content = QtWidgets.QWidget()
        left_content.setObjectName("buildScrollContent")
        left_layout = QtWidgets.QVBoxLayout(left_content)
        left_layout.setContentsMargins(0, 0, 10, 0)
        left_layout.setSpacing(12)
        left.setWidget(left_content)
        surface_palette = QtGui.QPalette(left.palette())
        surface_palette.setColor(
            self._palette_role("Window", "Window"), QtGui.QColor(THEME["panel"])
        )
        surface_palette.setColor(
            self._palette_role("Base", "Base"), QtGui.QColor(THEME["panel"])
        )
        for surface in (left, left.viewport(), left_content):
            surface.setPalette(surface_palette)
            surface.setAutoFillBackground(True)
        right = QtWidgets.QWidget()
        right.setObjectName("buildEvidencePane")
        right_layout = QtWidgets.QVBoxLayout(right)
        right_layout.setContentsMargins(4, 0, 0, 0)
        right_layout.setSpacing(10)
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes((510, 510))

        project_group = QtWidgets.QGroupBox("Project")
        project_form = QtWidgets.QFormLayout(project_group)
        self._configure_form(project_form)
        self.project_edit = QtWidgets.QLineEdit()
        project_form.addRow("Directory:", self._path_row(self.project_edit, "Browse…", self._browse_project))
        self.input_name_edit = QtWidgets.QLineEdit("mcmc.input")
        self.output_stem_edit = QtWidgets.QLineEdit("simulation")
        project_form.addRow("Input filename:", self.input_name_edit)
        project_form.addRow("Output prefix:", self.output_stem_edit)
        left_layout.addWidget(project_group)

        structure_group = QtWidgets.QGroupBox("Structure")
        structure_form = QtWidgets.QFormLayout(structure_group)
        self._configure_form(structure_form)
        self.source_mode = QtWidgets.QComboBox()
        self.source_mode.addItems(("Local PDB", "RCSB PDB", "Live PyMOL object"))
        self.source_mode.currentIndexChanged.connect(self._source_mode_changed)
        structure_form.addRow("Source:", self.source_mode)
        self.local_pdb_edit = QtWidgets.QLineEdit()
        structure_form.addRow("Local file:", self._path_row(self.local_pdb_edit, "Browse…", self._browse_local_pdb))
        rcsb_widget = QtWidgets.QWidget()
        rcsb_layout = QtWidgets.QHBoxLayout(rcsb_widget)
        rcsb_layout.setContentsMargins(0, 0, 0, 0)
        self.rcsb_edit = QtWidgets.QLineEdit()
        self.rcsb_edit.setMaxLength(4)
        fetch_button = QtWidgets.QPushButton("Fetch")
        fetch_button.clicked.connect(self._fetch_rcsb)
        rcsb_layout.addWidget(self.rcsb_edit)
        rcsb_layout.addWidget(fetch_button)
        structure_form.addRow("RCSB ID:", rcsb_widget)
        pymol_widget = QtWidgets.QWidget()
        pymol_layout = QtWidgets.QHBoxLayout(pymol_widget)
        pymol_layout.setContentsMargins(0, 0, 0, 0)
        self.pymol_object_combo = QtWidgets.QComboBox()
        self.pymol_object_combo.currentIndexChanged.connect(self._read_live_pymol_object)
        refresh_objects = QtWidgets.QPushButton("Refresh")
        refresh_objects.clicked.connect(self._refresh_pymol_objects)
        pymol_layout.addWidget(self.pymol_object_combo, 1)
        pymol_layout.addWidget(refresh_objects)
        structure_form.addRow("PyMOL object:", pymol_widget)
        self.follow_pymol_edits = QtWidgets.QCheckBox("Use current PyMOL coordinates")
        self.follow_pymol_edits.setToolTip(
            "Prepare the next run from the coordinates currently displayed in PyMOL."
        )
        self.follow_pymol_edits.setChecked(True)
        structure_form.addRow("Synchronization:", self.follow_pymol_edits)
        self.model_combo = QtWidgets.QComboBox()
        structure_form.addRow("PDB model:", self.model_combo)
        self.chain_list = QtWidgets.QListWidget()
        self.chain_list.setSelectionMode(EXTENDED_SELECTION)
        self.chain_list.setMaximumHeight(80)
        self.chain_list.itemSelectionChanged.connect(self._update_header_preview)
        structure_form.addRow("Chains:", self.chain_list)
        self.header_combo = QtWidgets.QComboBox()
        self.header_combo.addItem("No CBLC header", "none")
        self.header_combo.addItem("Regular CBLC", "regular")
        self.header_combo.addItem("Successive CBLC (SCBLC)", "successive")
        self.header_combo.currentIndexChanged.connect(self._update_header_preview)
        structure_form.addRow("Move header:", self.header_combo)
        self.header_preview = QtWidgets.QLineEdit()
        self.header_preview.setReadOnly(True)
        structure_form.addRow("Header written:", self.header_preview)
        self.structure_status = QtWidgets.QLabel("No structure loaded")
        self.structure_status.setWordWrap(True)
        structure_form.addRow("Live status:", self.structure_status)
        left_layout.addWidget(structure_group)

        self.disulfide_group = QtWidgets.QGroupBox("Protein disulfides")
        disulfide_layout = QtWidgets.QVBoxLayout(self.disulfide_group)
        self.disulfide_table = QtWidgets.QTableWidget(0, 5)
        self.disulfide_table.setHorizontalHeaderLabels(("Use", "Residue 1", "Residue 2", "SG–SG Å", "Evidence"))
        self.disulfide_table.horizontalHeader().setStretchLastSection(True)
        self.disulfide_table.setMaximumHeight(120)
        disulfide_layout.addWidget(self.disulfide_table)
        show_disulfides = QtWidgets.QPushButton("Show selected disulfides in PyMOL")
        show_disulfides.clicked.connect(self._show_disulfides)
        disulfide_layout.addWidget(show_disulfides)
        left_layout.addWidget(self.disulfide_group)

        settings_group = QtWidgets.QGroupBox("Runtime, force field, and analysis")
        settings_form = QtWidgets.QFormLayout(settings_group)
        self._configure_form(settings_form)
        self.runtime_combo = QtWidgets.QComboBox()
        self.runtime_combo.setSizeAdjustPolicy(COMBO_ADJUST_MINIMUM)
        self.runtime_combo.setMinimumContentsLength(12)
        for profile in RUNTIME_PROFILES:
            label = profile.label if profile.available() else profile.label + " — unavailable on this computer"
            self.runtime_combo.addItem(label, profile.identifier)
        self.runtime_combo.currentIndexChanged.connect(self._runtime_changed)
        settings_form.addRow("Runtime:", self.runtime_combo)
        self.runtime_description = QtWidgets.QLabel()
        self.runtime_description.setWordWrap(True)
        self._allow_label_to_wrap(self.runtime_description)
        settings_form.addRow("Runtime scope:", self.runtime_description)
        self.forcefield_combo = QtWidgets.QComboBox()
        self._populate_force_field_combo(self.forcefield_combo)
        self.forcefield_combo.currentIndexChanged.connect(self._forcefield_changed)
        settings_form.addRow("Force field / topology:", self.forcefield_combo)
        self.forcefield_description = QtWidgets.QLabel()
        self.forcefield_description.setWordWrap(True)
        self._allow_label_to_wrap(self.forcefield_description)
        settings_form.addRow("Definition:", self.forcefield_description)
        self.protein_prep_mode = QtWidgets.QComboBox()
        self.protein_prep_mode.addItem(
            "Automatic AMBER preparation (PDB2PQR + PROPKA)", "pdb2pqr"
        )
        self.protein_prep_mode.addItem(
            "Require an existing ff14SB-ready all-atom PDB", "strict"
        )
        self.protein_prep_mode.setToolTip(
            "Automatic mode repairs missing heavy atoms, assigns AMBER residue names, "
            "and adds hydrogens before strict RTF validation."
        )
        self.protein_prep_mode.currentIndexChanged.connect(
            self._update_protein_prep_controls
        )
        settings_form.addRow("Protein preparation:", self.protein_prep_mode)
        self.protein_ph = QtWidgets.QDoubleSpinBox()
        self.protein_ph.setRange(0.0, 14.0)
        self.protein_ph.setDecimals(2)
        self.protein_ph.setSingleStep(0.1)
        self.protein_ph.setValue(7.0)
        self.protein_ph.setToolTip(
            "PROPKA uses this pH when assigning titration states for PDB2PQR."
        )
        settings_form.addRow("Protonation pH:", self.protein_ph)
        self.protein_prep_status = QtWidgets.QLabel()
        self.protein_prep_status.setWordWrap(True)
        self._allow_label_to_wrap(self.protein_prep_status)
        settings_form.addRow("Preparation tool:", self.protein_prep_status)
        self.preset_combo = QtWidgets.QComboBox()
        self.preset_combo.setSizeAdjustPolicy(COMBO_ADJUST_MINIMUM)
        self.preset_combo.setMinimumContentsLength(12)
        for preset in ANALYSIS_PRESETS:
            self.preset_combo.addItem(preset.label, preset.identifier)
        self.preset_combo.currentIndexChanged.connect(self._preset_changed)
        settings_form.addRow("Analysis preset:", self.preset_combo)
        self.preset_description = QtWidgets.QLabel()
        self.preset_description.setWordWrap(True)
        self._allow_label_to_wrap(self.preset_description)
        settings_form.addRow("Purpose:", self.preset_description)
        self.temperature = QtWidgets.QDoubleSpinBox()
        self.temperature.setRange(0.0, 1000000.0)
        self.temperature.setDecimals(3)
        self.steps = QtWidgets.QSpinBox()
        self.steps.setRange(0, 2147483647)
        self.statistics = QtWidgets.QSpinBox()
        self.statistics.setRange(1, 2147483647)
        self.seed = QtWidgets.QLineEdit("-7143580450")
        self.closure_sigma = QtWidgets.QLineEdit("0.001")
        self.replica_count = QtWidgets.QSpinBox()
        self.replica_count.setRange(1, 64)
        self.top_temperature = QtWidgets.QDoubleSpinBox()
        self.top_temperature.setRange(0.0, 1000000.0)
        self.top_temperature.setDecimals(3)
        self.ladder_preview = QtWidgets.QLabel()
        self.ladder_preview.setWordWrap(True)
        self._allow_label_to_wrap(self.ladder_preview)
        self.temperature.valueChanged.connect(self._update_ladder)
        self.replica_count.valueChanged.connect(self._update_ladder)
        self.top_temperature.valueChanged.connect(self._update_ladder)
        settings_form.addRow("Temperature (K):", self.temperature)
        settings_form.addRow("Monte Carlo steps:", self.steps)
        settings_form.addRow("Statistics frequency:", self.statistics)
        settings_form.addRow("Random seed:", self.seed)
        settings_form.addRow("Closure σ:", self.closure_sigma)
        settings_form.addRow("Number of PT replicas:", self.replica_count)
        settings_form.addRow("Highest temperature (K):", self.top_temperature)
        settings_form.addRow("Temperature ladder:", self.ladder_preview)
        self.use_region = QtWidgets.QCheckBox("Use file")
        self.use_region.setToolTip("Write and reference region/region.data in this project.")
        region_button = QtWidgets.QPushButton("Edit…")
        region_button.clicked.connect(self._edit_region)
        region_row = QtWidgets.QWidget()
        region_layout = QtWidgets.QHBoxLayout(region_row)
        region_layout.setContentsMargins(0, 0, 0, 0)
        region_layout.addWidget(self.use_region)
        region_layout.addWidget(region_button)
        settings_form.addRow("Region:", region_row)
        left_layout.addWidget(settings_group)
        left_layout.addStretch(1)

        preview_title = QtWidgets.QLabel("MCMC INPUT · EDITABLE SOURCE OF TRUTH")
        preview_title.setObjectName("sectionLabel")
        right_layout.addWidget(preview_title)
        self.input_preview = QtWidgets.QPlainTextEdit()
        self.input_preview.setLineWrapMode(PLAIN_TEXT_NO_WRAP)
        right_layout.addWidget(self.input_preview, 1)
        buttons = QtWidgets.QHBoxLayout()
        generate = QtWidgets.QPushButton("Generate from visible settings")
        generate.clicked.connect(self._generate_preview)
        prepare = QtWidgets.QPushButton("Prepare project")
        prepare.setObjectName("primaryAction")
        prepare.clicked.connect(self._prepare_project)
        buttons.addWidget(generate)
        buttons.addWidget(prepare)
        right_layout.addLayout(buttons)
        self.build_status = QtWidgets.QPlainTextEdit()
        self.build_status.setReadOnly(True)
        self.build_status.setMaximumHeight(120)
        right_layout.addWidget(self.build_status)
        self.output_stem_edit.textChanged.connect(self._mark_preview_stale)
        self.seed.textChanged.connect(self._mark_preview_stale)
        self.closure_sigma.textChanged.connect(self._mark_preview_stale)
        self.temperature.valueChanged.connect(self._regenerate_unedited_preview)
        self.steps.valueChanged.connect(self._regenerate_unedited_preview)
        self.statistics.valueChanged.connect(self._regenerate_unedited_preview)
        self.replica_count.valueChanged.connect(self._regenerate_unedited_preview)
        self.top_temperature.valueChanged.connect(self._regenerate_unedited_preview)
        self.use_region.toggled.connect(self._regenerate_unedited_preview)
        self._runtime_changed()
        self._preset_changed()
        self._forcefield_changed()
        self._source_mode_changed()
        self._generate_preview()

    def _build_run_tab(self):
        layout = QtWidgets.QVBoxLayout(self.run_tab)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)
        form = QtWidgets.QFormLayout()
        self._configure_form(form)
        self.input_combo = QtWidgets.QComboBox()
        self.input_combo.setEditable(True)
        self.output_edit = QtWidgets.QLineEdit()
        input_row = QtWidgets.QWidget()
        input_layout = QtWidgets.QHBoxLayout(input_row)
        input_layout.setContentsMargins(0, 0, 0, 0)
        input_layout.addWidget(self.input_combo, 1)
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self._browse_input)
        scan = QtWidgets.QPushButton("Scan project")
        scan.clicked.connect(self._scan_inputs)
        view_input = QtWidgets.QPushButton("View / edit")
        view_input.clicked.connect(self._view_run_input)
        input_layout.addWidget(browse)
        input_layout.addWidget(scan)
        input_layout.addWidget(view_input)
        form.addRow("Input file:", input_row)
        form.addRow("Specific output PDB (optional):", self._path_row(self.output_edit, "Browse…", self._browse_output))
        view_structure = QtWidgets.QPushButton("View / edit input PDB")
        view_structure.clicked.connect(self._view_prepared_structure)
        form.addRow("Prepared structure:", view_structure)
        layout.addLayout(form)

        diagnostics = QtWidgets.QHBoxLayout()
        diagnostics.setSpacing(10)
        validation_group = QtWidgets.QGroupBox("Validation")
        validation_layout = QtWidgets.QVBoxLayout(validation_group)
        self.validation_output = QtWidgets.QPlainTextEdit()
        self.validation_output.setReadOnly(True)
        self.validation_output.setMinimumHeight(72)
        self.validation_output.setMaximumHeight(112)
        validation_layout.addWidget(self.validation_output)
        diagnostics.addWidget(validation_group, 1)

        command_group = QtWidgets.QGroupBox("Exact execution plan")
        command_layout = QtWidgets.QVBoxLayout(command_group)
        self.command_preview = QtWidgets.QPlainTextEdit()
        self.command_preview.setReadOnly(True)
        self.command_preview.setMinimumHeight(72)
        self.command_preview.setMaximumHeight(112)
        command_layout.addWidget(self.command_preview)
        diagnostics.addWidget(command_group, 1)
        layout.addLayout(diagnostics)

        run_controls = QtWidgets.QGroupBox("Run controls")
        run_controls_layout = QtWidgets.QVBoxLayout(run_controls)
        run_controls_layout.setSpacing(7)
        controls = QtWidgets.QHBoxLayout()
        controls.setSpacing(8)
        validate = QtWidgets.QPushButton("Validate")
        validate.clicked.connect(self._validate_current_project)
        self.run_button = QtWidgets.QPushButton("Run MOSAICS")
        self.run_button.setObjectName("primaryAction")
        self.run_button.clicked.connect(self._start_process)
        self.stop_button = QtWidgets.QPushButton("Stop")
        self.stop_button.setObjectName("stopAction")
        self.stop_button.clicked.connect(self._stop_process)
        self.stop_button.setEnabled(False)
        load = QtWidgets.QPushButton("Load all PDB outputs")
        load.clicked.connect(self._load_all_outputs)
        open_log = QtWidgets.QPushButton("View latest log")
        open_log.clicked.connect(self._view_latest_log)
        controls.addWidget(validate)
        controls.addWidget(self.run_button)
        controls.addWidget(self.stop_button)
        controls.addWidget(load)
        controls.addWidget(open_log)
        controls.addStretch(1)
        run_controls_layout.addLayout(controls)
        self.auto_load = QtWidgets.QCheckBox("Load all final structures and trajectories after a successful run")
        self.auto_load.setChecked(True)
        run_controls_layout.addWidget(self.auto_load)
        layout.addWidget(run_controls)

        output_group = QtWidgets.QGroupBox("Live MOSAICS output")
        output_layout = QtWidgets.QVBoxLayout(output_group)
        self.log_output = QtWidgets.QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setLineWrapMode(PLAIN_TEXT_NO_WRAP)
        self.log_output.setMinimumHeight(130)
        output_layout.addWidget(self.log_output)
        layout.addWidget(output_group, 1)

    def _build_analysis_tab(self):
        layout = QtWidgets.QVBoxLayout(self.analysis_tab)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)
        top = QtWidgets.QHBoxLayout()
        refresh = QtWidgets.QPushButton("Refresh project analysis")
        refresh.clicked.connect(self._refresh_analysis)
        top.addWidget(refresh)
        top.addWidget(QtWidgets.QLabel("Inspect the exact files, energies, acceptance, and sampled conformations."))
        top.addStretch(1)
        layout.addLayout(top)

        analysis_pages = QtWidgets.QTabWidget()
        analysis_pages.setObjectName("analysisPages")
        layout.addWidget(analysis_pages, 1)

        energy_page = QtWidgets.QWidget()
        energy_layout = QtWidgets.QVBoxLayout(energy_page)
        energy_row = QtWidgets.QHBoxLayout()
        self.energy_combo = QtWidgets.QComboBox()
        self.energy_combo.currentIndexChanged.connect(self._energy_changed)
        energy_row.addWidget(QtWidgets.QLabel("Energy series:"))
        energy_row.addWidget(self.energy_combo, 1)
        energy_layout.addLayout(energy_row)
        self.energy_summary = QtWidgets.QLabel("No analysis loaded")
        self.energy_summary.setWordWrap(True)
        energy_layout.addWidget(self.energy_summary)
        self.energy_plot = EnergyPlot()
        energy_layout.addWidget(self.energy_plot, 1)
        energy_layout.addWidget(QtWidgets.QLabel("Natural-move acceptance reported by MOSAICS"))
        self.acceptance_table = QtWidgets.QTableWidget(0, 4)
        self.acceptance_table.setHorizontalHeaderLabels(("Chain", "Accepted", "Attempted", "Ratio"))
        self.acceptance_table.horizontalHeader().setStretchLastSection(True)
        self.acceptance_table.setMaximumHeight(180)
        energy_layout.addWidget(self.acceptance_table)
        analysis_pages.addTab(energy_page, "Energy & acceptance")

        landscape_page = QtWidgets.QWidget()
        landscape_layout = QtWidgets.QVBoxLayout(landscape_page)
        landscape_intro = QtWidgets.QLabel(
            "The map preserves pairwise, rigid-body-aligned RMSD as closely as two dimensions allow. "
            "It is a structural projection—not proof of thermodynamic convergence or a free-energy surface."
        )
        landscape_intro.setWordWrap(True)
        landscape_layout.addWidget(landscape_intro)
        landscape_controls = QtWidgets.QHBoxLayout()
        self.trajectory_combo = QtWidgets.QComboBox()
        self.landscape_representatives = QtWidgets.QSpinBox()
        self.landscape_representatives.setRange(1, 24)
        self.landscape_representatives.setValue(6)
        self.landscape_maximum_frames = QtWidgets.QSpinBox()
        self.landscape_maximum_frames.setRange(2, 300)
        self.landscape_maximum_frames.setValue(200)
        self.landscape_energy_combo = QtWidgets.QComboBox()
        self.landscape_energy_combo.addItem("No energy coloring", "")
        build_map = QtWidgets.QPushButton("Build structural map")
        build_map.clicked.connect(self._build_landscape)
        landscape_controls.addWidget(QtWidgets.QLabel("Trajectory:"))
        landscape_controls.addWidget(self.trajectory_combo, 1)
        landscape_controls.addWidget(QtWidgets.QLabel("Representatives:"))
        landscape_controls.addWidget(self.landscape_representatives)
        landscape_controls.addWidget(QtWidgets.QLabel("Max frames:"))
        landscape_controls.addWidget(self.landscape_maximum_frames)
        landscape_controls.addWidget(build_map)
        landscape_layout.addLayout(landscape_controls)
        energy_color_row = QtWidgets.QHBoxLayout()
        energy_color_row.addWidget(QtWidgets.QLabel("Point color:"))
        energy_color_row.addWidget(self.landscape_energy_combo, 1)
        landscape_layout.addLayout(energy_color_row)
        self.landscape_status = QtWidgets.QLabel("No structural map built")
        self.landscape_status.setWordWrap(True)
        landscape_layout.addWidget(self.landscape_status)
        landscape_splitter = QtWidgets.QSplitter()
        self.landscape_plot = LandscapePlot()
        self.landscape_plot.on_frame_selected = self._landscape_frame_selected
        self.representative_list = QtWidgets.QListWidget()
        self.representative_list.itemDoubleClicked.connect(self._representative_selected)
        representative_widget = QtWidgets.QWidget()
        representative_layout = QtWidgets.QVBoxLayout(representative_widget)
        representative_layout.addWidget(QtWidgets.QLabel("Representative frames"))
        representative_layout.addWidget(self.representative_list)
        representative_layout.addWidget(QtWidgets.QLabel("Click a point or double-click a frame to show it in PyMOL."))
        landscape_splitter.addWidget(self.landscape_plot)
        landscape_splitter.addWidget(representative_widget)
        landscape_splitter.setSizes((760, 240))
        landscape_layout.addWidget(landscape_splitter, 1)
        analysis_pages.addTab(landscape_page, "Structural landscape")

        files_page = QtWidgets.QWidget()
        output_layout = QtWidgets.QVBoxLayout(files_page)
        self.output_list = QtWidgets.QListWidget()
        self.output_list.itemDoubleClicked.connect(self._open_analysis_item)
        output_layout.addWidget(QtWidgets.QLabel("Project files and simulation outputs"))
        output_layout.addWidget(self.output_list)
        output_buttons = QtWidgets.QHBoxLayout()
        view_selected = QtWidgets.QPushButton("View selected as text")
        view_selected.clicked.connect(self._view_selected_analysis_output)
        load_selected = QtWidgets.QPushButton("Load selected PDB in PyMOL")
        load_selected.clicked.connect(self._load_selected_analysis_outputs)
        output_buttons.addWidget(view_selected)
        output_buttons.addWidget(load_selected)
        output_layout.addLayout(output_buttons)
        analysis_pages.addTab(files_page, "Files & logs")

    def _build_setup_tab(self):
        outer = QtWidgets.QVBoxLayout(self.setup_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setObjectName("setupScroll")
        scroll.setFrameShape(FRAME_NO_FRAME)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(SCROLLBAR_ALWAYS_OFF)
        content = QtWidgets.QWidget()
        content.setObjectName("setupScrollContent")
        content_policy = content.sizePolicy()
        content_policy.setVerticalPolicy(SIZE_POLICY_IGNORED)
        content.setSizePolicy(content_policy)
        layout = QtWidgets.QVBoxLayout(content)
        layout.setContentsMargins(18, 14, 18, 14)
        layout.setSpacing(12)
        scroll.setWidget(content)
        surface_palette = QtGui.QPalette(scroll.palette())
        surface_palette.setColor(
            self._palette_role("Window", "Window"), QtGui.QColor(THEME["panel"])
        )
        surface_palette.setColor(
            self._palette_role("Base", "Base"), QtGui.QColor(THEME["panel"])
        )
        for surface in (scroll, scroll.viewport(), content):
            surface.setPalette(surface_palette)
            surface.setAutoFillBackground(True)
        outer.addWidget(scroll)
        explanation = QtWidgets.QLabel(
            "Stable and experimental executables are bundled for Apple-Silicon macOS. "
            "Windows, Linux, and Intel-macOS users can select a compatible custom executable."
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)
        runtime_group = QtWidgets.QGroupBox("MOSAICS runtime")
        form = QtWidgets.QFormLayout(runtime_group)
        self._configure_form(form)
        self.custom_executable_edit = QtWidgets.QLineEdit()
        self.workspace_edit = QtWidgets.QLineEdit()
        form.addRow("Custom MOSAICS executable:", self._path_row(self.custom_executable_edit, "Browse…", self._browse_executable))
        form.addRow("Default workspace:", self._path_row(self.workspace_edit, "Browse…", self._browse_workspace))
        layout.addWidget(runtime_group)

        protein_tools = QtWidgets.QGroupBox("Protein preparation")
        protein_tools_form = QtWidgets.QFormLayout(protein_tools)
        self._configure_form(protein_tools_form)
        self.pdb2pqr_edit = QtWidgets.QLineEdit()
        detected_pdb2pqr = discover_pdb2pqr()
        if detected_pdb2pqr is not None:
            self.pdb2pqr_edit.setText(str(detected_pdb2pqr))
        protein_tools_form.addRow(
            "PDB2PQR executable:",
            self._path_row(self.pdb2pqr_edit, "Browse…", self._browse_pdb2pqr),
        )
        protein_note = QtWidgets.QLabel(
            "Used only when automatic protein preparation is selected. PymoSAICS runs "
            "PDB2PQR with AMBER naming and PROPKA at the visible pH, retains the PQR/log, "
            "and still requires an exact ff14SB topology match before MOSAICS can run."
        )
        protein_note.setWordWrap(True)
        protein_tools_form.addRow("Method:", protein_note)
        layout.addWidget(protein_tools)

        force_field_group = QtWidgets.QGroupBox("Force-field profile")
        force_field_layout = QtWidgets.QVBoxLayout(force_field_group)
        force_field_form = QtWidgets.QFormLayout()
        self._configure_form(force_field_form)
        self.setup_forcefield_combo = QtWidgets.QComboBox()
        self._populate_force_field_combo(self.setup_forcefield_combo)
        current_index = self.setup_forcefield_combo.findData(self.forcefield_combo.currentData())
        if current_index >= 0:
            self.setup_forcefield_combo.setCurrentIndex(current_index)
        self.setup_forcefield_combo.currentIndexChanged.connect(self._setup_forcefield_changed)
        force_field_form.addRow("Profile:", self.setup_forcefield_combo)
        bundled = QtWidgets.QLineEdit(str(FORCEFIELD_ROOT))
        bundled.setReadOnly(True)
        force_field_form.addRow("Bundled library:", bundled)
        force_field_layout.addLayout(force_field_form)
        profile_count = QtWidgets.QLabel(
            "{} complete profiles are bundled: {} DNA/RNA and {} protein. "
            "This includes AMBER99-based legacy profiles.".format(
                len(FORCE_FIELD_PROFILES),
                sum(profile.chemistry == "nucleic_acid" for profile in FORCE_FIELD_PROFILES),
                sum(profile.chemistry == "protein" for profile in FORCE_FIELD_PROFILES),
            )
        )
        profile_count.setWordWrap(True)
        force_field_layout.addWidget(profile_count)
        self.setup_forcefield_description = QtWidgets.QLabel()
        self.setup_forcefield_description.setWordWrap(True)
        force_field_layout.addWidget(self.setup_forcefield_description)
        force_field_layout.addWidget(QtWidgets.QLabel("Files written automatically into mcmc.input"))
        self.setup_forcefield_files = QtWidgets.QPlainTextEdit()
        self.setup_forcefield_files.setReadOnly(True)
        self.setup_forcefield_files.setMaximumHeight(118)
        force_field_layout.addWidget(self.setup_forcefield_files)
        layout.addWidget(force_field_group)
        self._update_setup_force_field_details(self._current_forcefield())

        self.setup_status = QtWidgets.QPlainTextEdit()
        self.setup_status.setReadOnly(True)
        self.setup_status.setMaximumHeight(78)
        layout.addWidget(self.setup_status)
        save = QtWidgets.QPushButton("Validate and save")
        save.setObjectName("primaryAction")
        save.clicked.connect(self._save_configuration)
        save_row = QtWidgets.QHBoxLayout()
        save_row.addStretch(1)
        save_row.addWidget(save)
        layout.addLayout(save_row)
        layout.addStretch(1)
        self._runtime_changed()
        self._update_protein_prep_controls()

    def _build_about_tab(self):
        layout = QtWidgets.QVBoxLayout(self.about_tab)
        layout.setContentsMargins(24, 22, 24, 22)
        label = QtWidgets.QLabel(
            "<h2>PymoSAICS 0.2.3</h2>"
            "<p>A transparent PyMOL workbench for MOSAICS structure preparation, "
            "input generation, execution, visualization, and analysis.</p>"
            "<p><b>MOSAICS</b> was created by Peter Minary. Obtain official executables, "
            "documentation, and licensing information from the "
            "<a href='https://www.cs.ox.ac.uk/mosaics/Downloads.php'>MOSAICS Downloads page</a>.</p>"
            "<p>The original PyMOSAICS graphical interface is credited to Konrad Krawczyk. "
            "This Python 3 / Qt redesign and current PymoSAICS release are by "
            "Folorunsho Bright Omage.</p>"
            "<p>The plugin source is MIT licensed. Bundled MOSAICS executables and "
            "force-field data retain their own provenance and redistribution terms.</p>"
            "<p>Every generated PDB header, region file, and mcmc.input remains visible and editable.</p>"
        )
        label.setWordWrap(True)
        label.setOpenExternalLinks(True)
        layout.addWidget(label)
        layout.addStretch(1)

    def _current_runtime(self):
        return runtime_profile(self.runtime_combo.currentData())

    def _current_forcefield(self):
        return force_field_profile(self.forcefield_combo.currentData())

    def _current_preset(self):
        return analysis_preset(self.preset_combo.currentData())

    def _configuration_from_fields(self) -> RuntimeConfig:
        runtime = self._current_runtime()
        executable = runtime.executable() if runtime.identifier != "custom" else Path(self.custom_executable_edit.text().strip()).expanduser()
        if executable is None:
            executable = Path("")
        workspace = self.workspace_edit.text().strip()
        pdb2pqr_text = (
            self.pdb2pqr_edit.text().strip() if hasattr(self, "pdb2pqr_edit") else ""
        )
        return RuntimeConfig(
            executable=executable,
            forcefield_directory=FORCEFIELD_ROOT,
            default_workspace=Path(workspace).expanduser() if workspace else None,
            runtime_id=runtime.identifier,
            force_field_id=self._current_forcefield().identifier,
            pdb2pqr_executable=(Path(pdb2pqr_text).expanduser() if pdb2pqr_text else None),
        )

    def _load_configuration(self):
        try:
            config = self.config_store.load()
        except ConfigError as exc:
            self.setup_status.setPlainText("ERROR: {}".format(exc))
            return
        if config is None:
            preferred = "mosaics-experimental-2026-07-21"
            index = self.runtime_combo.findData(preferred)
            if index >= 0 and runtime_profile(preferred).available():
                self.runtime_combo.setCurrentIndex(index)
            self.setup_status.setPlainText("Bundled defaults are ready. Choose a workspace and save if desired.")
            return
        runtime_index = self.runtime_combo.findData(config.runtime_id)
        if runtime_index >= 0:
            self.runtime_combo.setCurrentIndex(runtime_index)
        forcefield_index = self.forcefield_combo.findData(config.force_field_id)
        if forcefield_index >= 0:
            self.forcefield_combo.setCurrentIndex(forcefield_index)
        if config.runtime_id == "custom":
            self.custom_executable_edit.setText(str(config.executable))
        if config.pdb2pqr_executable is not None:
            self.pdb2pqr_edit.setText(str(config.pdb2pqr_executable))
        self.workspace_edit.setText(str(config.default_workspace or ""))
        if not self.project_edit.text() and config.default_workspace:
            self.project_edit.setText(str(config.default_workspace))
        self.setup_status.setPlainText(format_diagnostics(validate_runtime(self._configuration_from_fields())))
        self._update_protein_prep_controls()

    def _save_configuration(self, checked=False):
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

    def _runtime_changed(self, *_args):
        profile = self._current_runtime()
        self.runtime_combo.setToolTip(profile.label)
        self.runtime_description.setText(profile.description)
        self.header_runtime.setText("Runtime · {}".format(profile.label))
        if profile.relative_executable and profile.available():
            make_bundled_executable_runnable(profile)
        if profile.identifier == "custom":
            self.tabs.setTabToolTip(self.tabs.indexOf(self.setup_tab), "Choose the custom executable here")
        force_field_combos = [self.forcefield_combo]
        if hasattr(self, "setup_forcefield_combo"):
            force_field_combos.append(self.setup_forcefield_combo)
        for combo in force_field_combos:
            for index, force_field in enumerate(FORCE_FIELD_PROFILES):
                compatible = runtime_supports_force_field(profile.identifier, force_field.identifier)
                item = combo.model().item(index)
                if item is not None:
                    item.setEnabled(compatible)
        current = self._current_forcefield()
        if not runtime_supports_force_field(profile.identifier, current.identifier):
            compatible_index = next(
                index
                for index, force_field in enumerate(FORCE_FIELD_PROFILES)
                if runtime_supports_force_field(profile.identifier, force_field.identifier)
            )
            self.forcefield_combo.setCurrentIndex(compatible_index)

    def _forcefield_changed(self, *_args):
        profile = self._current_forcefield()
        self.forcefield_combo.setToolTip(profile.label)
        if hasattr(self, "setup_forcefield_combo"):
            setup_index = self.setup_forcefield_combo.findData(profile.identifier)
            if setup_index >= 0 and setup_index != self.setup_forcefield_combo.currentIndex():
                self.setup_forcefield_combo.blockSignals(True)
                self.setup_forcefield_combo.setCurrentIndex(setup_index)
                self.setup_forcefield_combo.blockSignals(False)
            self._update_setup_force_field_details(profile)
        self.header_science.setText("Force field · {}".format(profile.label))
        self.forcefield_description.setText("{} {}".format(profile.description, profile.validation))
        self.disulfide_group.setVisible(profile.chemistry == "protein")
        self._update_protein_prep_controls()
        compatible = []
        for index, preset in enumerate(ANALYSIS_PRESETS):
            compatible.append(index if preset.chemistry in ("any", profile.chemistry) else -1)
            item = self.preset_combo.model().item(index)
            if item is not None:
                item.setEnabled(preset.chemistry in ("any", profile.chemistry))
        current = self._current_preset()
        if current.chemistry not in ("any", profile.chemistry):
            first = next(index for index in compatible if index >= 0)
            self.preset_combo.setCurrentIndex(first)
        self._refresh_disulfides()
        self._regenerate_unedited_preview()

    def _update_protein_prep_controls(self, *_args):
        if not hasattr(self, "protein_prep_mode"):
            return
        is_protein = self._current_forcefield().chemistry == "protein"
        fields = (self.protein_prep_mode, self.protein_ph, self.protein_prep_status)
        form = self.protein_prep_mode.parentWidget().layout()
        for field in fields:
            field.setVisible(is_protein)
            if isinstance(form, QtWidgets.QFormLayout):
                label = form.labelForField(field)
                if label is not None:
                    label.setVisible(is_protein)
        if not is_protein:
            return
        automatic = self.protein_prep_mode.currentData() == "pdb2pqr"
        self.protein_ph.setEnabled(automatic)
        explicit = None
        if hasattr(self, "pdb2pqr_edit") and self.pdb2pqr_edit.text().strip():
            explicit = Path(self.pdb2pqr_edit.text().strip()).expanduser()
        tool = discover_pdb2pqr(explicit)
        if not automatic:
            self.protein_prep_status.setText(
                "Strict mode: the selected PDB must already contain exact ff14SB atoms and names."
            )
        elif tool is None:
            self.protein_prep_status.setText(
                "PDB2PQR was not found. Choose its executable in Setup before preparing a heavy-atom protein."
            )
        else:
            self.protein_prep_status.setText("Detected: {}".format(tool))

    def _setup_forcefield_changed(self, *_args):
        identifier = self.setup_forcefield_combo.currentData()
        build_index = self.forcefield_combo.findData(identifier)
        if build_index >= 0 and build_index != self.forcefield_combo.currentIndex():
            self.forcefield_combo.setCurrentIndex(build_index)
        else:
            self._update_setup_force_field_details(self._current_forcefield())

    def _update_setup_force_field_details(self, profile):
        if not hasattr(self, "setup_forcefield_description"):
            return
        self.setup_forcefield_combo.setToolTip(profile.label)
        self.setup_forcefield_description.setText(
            "{} {} Selecting this profile populates all six force-field directives; "
            "individual database files do not need to be chosen manually.".format(
                profile.description, profile.validation
            )
        )
        self.setup_forcefield_files.setPlainText(self._force_field_bundle_summary(profile))

    def _preset_changed(self, *_args):
        preset = self._current_preset()
        self.preset_combo.setToolTip(preset.label)
        self.preset_description.setText(preset.description)
        settings = default_settings(preset)
        self.temperature.setValue(settings.temperature)
        self.steps.setValue(settings.total_steps)
        self.statistics.setValue(settings.statistics_frequency)
        self.seed.setText(str(settings.random_seed))
        self.closure_sigma.setText(settings.closure_sigma)
        self.replica_count.setValue(preset.replica_count)
        self.top_temperature.setValue(preset.top_temperature)
        if preset.identifier in ("landscape-pilot", "protein-landscape-pilot"):
            self._apply_landscape_suggestion()
        else:
            self._update_ladder()
        index = self.header_combo.findData(preset.recommended_header)
        if index >= 0:
            self.header_combo.setCurrentIndex(index)
        self._regenerate_unedited_preview()

    def _source_mode_changed(self, *_args):
        mode = self.source_mode.currentText()
        if mode == "Live PyMOL object":
            self._read_live_pymol_object()

    def _browse_project(self, checked=False):
        start = self.project_edit.text().strip() or self.workspace_edit.text().strip() or str(Path.home())
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select MOSAICS project directory", start)
        if path:
            self.project_edit.setText(path)
            self._scan_inputs()

    def _browse_local_pdb(self, checked=False):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select PDB", str(Path.home()), "PDB files (*.pdb);;All files (*)")
        if path:
            self.local_pdb_edit.setText(path)
            self.source_mode.setCurrentText("Local PDB")
            self._inspect_source(Path(path), load_into_pymol=True)

    def _fetch_rcsb(self, checked=False):
        project = self._project_directory()
        if project is None:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Choose a project directory before fetching from RCSB.")
            return
        identifier = self.rcsb_edit.text().strip().upper()
        destination = project / ".pymosaics" / "sources" / (identifier + ".pdb")
        try:
            path = fetch_rcsb_pdb(identifier, destination)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "PymoSAICS", "RCSB fetch failed: {}".format(exc))
            return
        self.source_mode.setCurrentText("RCSB PDB")
        self._inspect_source(path, load_into_pymol=True)

    def _refresh_pymol_objects(self, checked=False):
        try:
            objects = list(cmd.get_object_list("all"))
        except Exception:
            objects = []
        displayed = [
            self.pymol_object_combo.itemText(index)
            for index in range(self.pymol_object_combo.count())
        ]
        if objects == displayed:
            return
        current = self.pymol_object_combo.currentText()
        self.pymol_object_combo.blockSignals(True)
        self.pymol_object_combo.clear()
        self.pymol_object_combo.addItems(objects)
        if current in objects:
            self.pymol_object_combo.setCurrentText(current)
        self.pymol_object_combo.blockSignals(False)

    def _read_live_pymol_object(self, *_args):
        if self.source_mode.currentText() != "Live PyMOL object" and not self.follow_pymol_edits.isChecked():
            return
        name = self.pymol_object_combo.currentText()
        if not name:
            return
        try:
            text = cmd.get_pdbstr(name, state=0)
        except Exception as exc:
            self.structure_status.setText("Cannot read PyMOL object: {}".format(exc))
            return
        digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        object_changed = name != self._live_object_name
        if not object_changed and self._live_digest and digest == self._live_digest:
            return
        try:
            metadata = inspect_pdb_text(text)
        except Exception as exc:
            self.structure_status.setText("Cannot read PyMOL object: {}".format(exc))
            return
        changed = bool(self._live_digest) and digest != self._live_digest
        self._source_path = None
        self._source_text = text
        self._metadata = metadata
        self._live_digest = digest
        self._live_object_name = name
        self._populate_structure_choices(preserve=changed and not object_changed)
        if self._current_preset().identifier in ("landscape-pilot", "protein-landscape-pilot") and not changed:
            self._apply_landscape_suggestion()
        self.structure_status.setText(
            "Live object {} synchronized{} — {} model(s)".format(name, " (updated)" if changed else "", len(metadata.models))
        )

    def _poll_pymol(self):
        self._refresh_pymol_objects()
        if self.source_mode.currentText() == "Live PyMOL object" or self.follow_pymol_edits.isChecked():
            self._read_live_pymol_object()

    def _inspect_source(self, path: Path, load_into_pymol: bool):
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
            metadata = inspect_pdb_text(text)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "PymoSAICS", "Cannot read structure: {}".format(exc))
            return
        self._source_path = path.resolve()
        self._source_text = text
        self._metadata = metadata
        self._live_digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        self._live_object_name = ""
        self._populate_structure_choices()
        if self._current_preset().identifier in ("landscape-pilot", "protein-landscape-pilot"):
            self._apply_landscape_suggestion()
        self.structure_status.setText("Loaded {} — {} model(s)".format(path.name, len(metadata.models)))
        if load_into_pymol:
            name = _safe_name("pymosaics_source_" + path.stem)
            try:
                cmd.load(str(path), name)
                self._refresh_pymol_objects()
                self.pymol_object_combo.setCurrentText(name)
            except Exception as exc:
                self.structure_status.setText(self.structure_status.text() + "; PyMOL load failed: " + str(exc))

    def _populate_structure_choices(self, preserve=False):
        if self._metadata is None:
            return
        previous_model = self.model_combo.currentText()
        previous_chains = set(self._selected_chains())
        self.model_combo.blockSignals(True)
        self.model_combo.clear()
        self.model_combo.addItems(self._metadata.models)
        if preserve and previous_model in self._metadata.models:
            self.model_combo.setCurrentText(previous_model)
        self.model_combo.blockSignals(False)
        try:
            self.model_combo.currentTextChanged.disconnect(self._model_changed)
        except (TypeError, RuntimeError):
            pass
        self.model_combo.currentTextChanged.connect(self._model_changed)
        self._populate_chains(previous_chains if preserve else None)

    def _model_changed(self, *_args):
        self._populate_chains()

    def _populate_chains(self, preserve=None):
        if self._metadata is None:
            return
        chains = self._metadata.chains_by_model.get(self.model_combo.currentText(), ())
        self.chain_list.clear()
        for chain in chains:
            item = QtWidgets.QListWidgetItem(chain or "(blank)")
            item.setData(USER_ROLE, chain)
            self.chain_list.addItem(item)
            if preserve is None or chain in preserve:
                item.setSelected(True)
        self._update_header_preview()
        self._refresh_disulfides()

    def _selected_chains(self) -> Tuple[str, ...]:
        return tuple(item.data(USER_ROLE) for item in self.chain_list.selectedItems())

    def _update_header_preview(self, *_args):
        mode = self.header_combo.currentData()
        chains = "".join(chain or "A" for chain in self._selected_chains())
        self.header_preview.setText("" if mode == "none" else "CBLC {}{}".format(">" if mode == "regular" else "~", chains))

    def _refresh_disulfides(self):
        self.disulfide_table.setRowCount(0)
        if self._current_forcefield().chemistry != "protein" or not self._source_text:
            return
        try:
            candidates = detect_disulfides_text(self._source_text, self.model_combo.currentText() or "1")
        except ValueError:
            candidates = ()
        recommended = set(unambiguous_disulfide_keys(candidates))
        for row, candidate in enumerate(candidates):
            self.disulfide_table.insertRow(row)
            use = QtWidgets.QTableWidgetItem()
            use.setFlags(use.flags() | ITEM_USER_CHECKABLE)
            use.setCheckState(CHECKED if candidate.key in recommended else UNCHECKED)
            use.setData(USER_ROLE, candidate)
            self.disulfide_table.setItem(row, 0, use)
            self.disulfide_table.setItem(row, 1, QtWidgets.QTableWidgetItem(candidate.first.selector))
            self.disulfide_table.setItem(row, 2, QtWidgets.QTableWidgetItem(candidate.second.selector))
            self.disulfide_table.setItem(row, 3, QtWidgets.QTableWidgetItem("{:.3f}".format(candidate.distance_angstrom)))
            evidence = candidate.evidence
            if candidate.key not in recommended:
                evidence += "; ambiguous sulfur neighborhood—choose at most one pairing"
            self.disulfide_table.setItem(row, 4, QtWidgets.QTableWidgetItem(evidence))

    def _selected_disulfides(self) -> Tuple[Tuple[str, str], ...]:
        pairs = []
        for row in range(self.disulfide_table.rowCount()):
            item = self.disulfide_table.item(row, 0)
            if item.checkState() == CHECKED:
                candidate = item.data(USER_ROLE)
                pairs.append(candidate.key)
        return tuple(pairs)

    def _pymol_source_object(self) -> str:
        return self.pymol_object_combo.currentText()

    def _show_disulfides(self, checked=False):
        object_name = self._pymol_source_object()
        if not object_name:
            return
        try:
            cmd.delete("pymosaics_disulfides")
            cmd.delete("pymosaics_disulfide_bonds")
            residue_expressions = []
            for index, (first, second) in enumerate(self._selected_disulfides(), start=1):
                first_chain, first_resi = first.split(":", 1)
                second_chain, second_resi = second.split(":", 1)
                a = "({} and chain {} and resi {} and name SG)".format(object_name, first_chain, first_resi)
                b = "({} and chain {} and resi {} and name SG)".format(object_name, second_chain, second_resi)
                residue_expressions.extend((a, b))
                cmd.distance("pymosaics_disulfide_bonds", a, b)
            expression = " or ".join(residue_expressions) or "none"
            cmd.select("pymosaics_disulfides", expression)
            cmd.show("sticks", "byres pymosaics_disulfides")
            cmd.color("yellow", "pymosaics_disulfides")
            cmd.zoom("byres pymosaics_disulfides")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Cannot show disulfides: {}".format(exc))

    def _available_residues(self) -> Tuple[str, ...]:
        if self._metadata is None:
            return ()
        model = self.model_combo.currentText()
        chains = set(self._selected_chains())
        return tuple(item.selector for item in self._metadata.residues_by_model.get(model, ()) if item.chain in chains)

    def _edit_region(self, checked=False):
        residues = self._available_residues()
        if not residues:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Load a structure and select chains first.")
            return
        dialog = RegionEditor(residues, self._pymol_source_object(), self._region_settings, self)
        if dialog.exec() == DIALOG_ACCEPTED and dialog.settings:
            self._region_settings = dialog.settings
            self.use_region.setChecked(True)

    def _current_input_settings(self) -> InputSettings:
        count = self.replica_count.value()
        low = self.temperature.value()
        high = self.top_temperature.value()
        if count > 1 and (low <= 0 or high < low):
            raise ValueError("highest temperature must be at least the positive base temperature")
        gap = (high / low) ** (1.0 / (count - 1)) if count > 1 else 1.1
        return InputSettings(
            temperature=self.temperature.value(),
            total_steps=self.steps.value(),
            statistics_frequency=self.statistics.value(),
            random_seed=int(self.seed.text().strip()),
            closure_sigma=self.closure_sigma.text().strip(),
            replica_number=count - 1,
            energy_gap=gap,
        )

    def _apply_landscape_suggestion(self):
        atom_count = sum(
            1 for line in self._source_text.splitlines() if line[:6].strip().upper() in ("ATOM", "HETATM")
        )
        if atom_count <= 1000:
            replicas = 6
        elif atom_count <= 5000:
            replicas = 8
        else:
            replicas = 10
        top = 529.0 if self._current_forcefield().chemistry == "protein" else 500.0
        self.replica_count.setValue(replicas)
        self.top_temperature.setValue(top)
        self._update_ladder()
        if atom_count:
            self.preset_description.setText(
                self._current_preset().description
                + " Initial suggestion: {} replicas for {} atoms; validate exchange rates before production.".format(
                    replicas, atom_count
                )
            )

    def _update_ladder(self, *_args):
        count = self.replica_count.value()
        low = self.temperature.value()
        high = self.top_temperature.value()
        if count == 1:
            self.ladder_preview.setText("{:g} K (single replica)".format(low))
            return
        if low <= 0 or high < low:
            self.ladder_preview.setText("Highest temperature must be at least the base temperature.")
            return
        gap = (high / low) ** (1.0 / (count - 1))
        values = [low * (gap ** index) for index in range(count)]
        self.ladder_preview.setText(
            " · ".join("{:g}".format(value) for value in values) + " K  (energy_gap={:.6g})".format(gap)
        )

    def _preview_settings_signature(self):
        settings = self._current_input_settings()
        return (
            self._current_preset().identifier,
            self._current_forcefield().identifier,
            self.output_stem_edit.text().strip(),
            settings.temperature,
            settings.total_steps,
            settings.statistics_frequency,
            settings.random_seed,
            settings.closure_sigma,
            settings.replica_number,
            settings.energy_gap,
            self.use_region.isChecked(),
        )

    def _mark_preview_stale(self, *_args):
        if hasattr(self, "build_status"):
            self.build_status.setPlainText(
                "Visible settings changed. Choose Generate from visible settings before preparing the project."
            )

    def _regenerate_unedited_preview(self, *_args):
        if not hasattr(self, "input_preview"):
            return
        if self.input_preview.document().isModified():
            self._mark_preview_stale()
            return
        self._generate_preview()

    def _generate_preview(self, checked=False):
        try:
            signature = self._preview_settings_signature()
            content = generate_mcmc_input(
                self._current_preset(),
                self._current_forcefield(),
                "structure.pdb",
                self.output_stem_edit.text().strip(),
                self._current_input_settings(),
                "region/region.data" if self.use_region.isChecked() else "",
            )
        except (ValueError, KeyError) as exc:
            self.build_status.setPlainText("Cannot generate input: {}".format(exc))
            return
        self.input_preview.setPlainText(content)
        self.input_preview.document().setModified(False)
        self._preview_signature = signature
        self.build_status.setPlainText("Generated from the displayed settings. Review and edit before preparing the project.")

    def _project_directory(self) -> Optional[Path]:
        text = self.project_edit.text().strip()
        return Path(text).expanduser().resolve() if text else None

    def _source_for_preparation(self, project: Path) -> Path:
        use_pymol = self.source_mode.currentText() == "Live PyMOL object" or (
            self.follow_pymol_edits.isChecked() and bool(self._pymol_source_object())
        )
        if use_pymol:
            self._read_live_pymol_object()
            if not self._source_text:
                raise ValueError("no live PyMOL object is selected")
            source = project / ".pymosaics" / "sources" / (_safe_name(self._pymol_source_object()) + ".pdb")
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_text(self._source_text, encoding="utf-8")
            return source
        if self._source_path and self._source_path.is_file():
            return self._source_path
        local = Path(self.local_pdb_edit.text().strip()).expanduser()
        if local.is_file():
            return local.resolve()
        raise ValueError("load or fetch a PDB structure first")

    @staticmethod
    def _concise_topology_issues(issues, maximum=4):
        missing = sum(len(issue.missing_atoms) for issue in issues)
        unexpected = sum(len(issue.extra_atoms) for issue in issues)
        duplicates = sum(len(issue.duplicate_atoms) for issue in issues)
        lines = [
            "{} residue mismatch(es): {} missing, {} unexpected, and {} duplicate atom name(s).".format(
                len(issues), missing, unexpected, duplicates
            )
        ]
        for issue in issues[:maximum]:
            details = []
            for label, atoms in (
                ("missing", issue.missing_atoms),
                ("unexpected", issue.extra_atoms),
                ("duplicate", issue.duplicate_atoms),
            ):
                if atoms:
                    preview = ", ".join(atoms[:5])
                    suffix = " …" if len(atoms) > 5 else ""
                    details.append("{} {}{}".format(label, preview, suffix))
            lines.append(
                "{} {} → {}: {}".format(
                    issue.selector,
                    issue.pdb_residue,
                    issue.topology_residue,
                    "; ".join(details),
                )
            )
        if len(issues) > maximum:
            lines.append("… {} additional residue mismatch(es)".format(len(issues) - maximum))
        return "\n".join(lines)

    def _prepare_project(self, checked=False):
        project = self._project_directory()
        if project is None:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Choose a project directory.")
            return None
        project.mkdir(parents=True, exist_ok=True)
        if not self._selected_chains():
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Select at least one chain.")
            return None
        try:
            input_name = self.input_name_edit.text().strip()
            if not input_name or Path(input_name).name != input_name or "\\" in input_name:
                raise ValueError("input filename must be a simple filename inside the project")
            if self._preview_signature != self._preview_settings_signature():
                raise ValueError(
                    "visible settings changed after this mcmc.input was generated; "
                    "choose Generate from visible settings, then reapply any deliberate text edits"
                )
            source = self._source_for_preparation(project)
            force_field = self._current_forcefield()
            staged_force_field = stage_force_field(project, force_field)
            selected_chains = self._selected_chains()
            selected_disulfides = self._selected_disulfides()
            prepared = prepare_structure(
                source,
                project / "structure.pdb",
                self.model_combo.currentText(),
                selected_chains,
                force_field.chemistry,
                force_field.topology_profile,
                self.header_combo.currentData(),
                selected_disulfides,
            )
            topology_issues = validate_pdb_against_rtf(
                prepared.pdb_path,
                staged_force_field / Path(force_field.rtf).name,
                force_field.chemistry,
            )
            protein_preparation = None
            if (
                topology_issues
                and force_field.chemistry == "protein"
                and self.protein_prep_mode.currentData() == "pdb2pqr"
            ):
                explicit = (
                    Path(self.pdb2pqr_edit.text().strip()).expanduser()
                    if self.pdb2pqr_edit.text().strip()
                    else None
                )
                pdb2pqr = discover_pdb2pqr(explicit)
                if pdb2pqr is None:
                    raise ValueError(
                        "This protein is not ff14SB-ready and PDB2PQR was not found. "
                        "Choose a PDB2PQR executable in Setup, or select strict mode and load "
                        "an AMBER-prepared all-atom PDB. PyMOL's generic Add Hydrogens command "
                        "does not produce ff14SB atom names."
                    )
                protein_preparation = prepare_protein_with_pdb2pqr(
                    source,
                    project / ".pymosaics" / "protein-preparation",
                    pdb2pqr,
                    self.model_combo.currentText(),
                    selected_chains,
                    ph=self.protein_ph.value(),
                )
                prepared = prepare_structure(
                    protein_preparation.pdb_path,
                    project / "structure.pdb",
                    "1",
                    selected_chains,
                    force_field.chemistry,
                    force_field.topology_profile,
                    self.header_combo.currentData(),
                    selected_disulfides,
                )
                topology_issues = validate_pdb_against_rtf(
                    prepared.pdb_path,
                    staged_force_field / Path(force_field.rtf).name,
                    force_field.chemistry,
                )
            if topology_issues:
                summary = self._concise_topology_issues(topology_issues)
                if force_field.chemistry == "protein" and protein_preparation is None:
                    raise ValueError(
                        "The selected protein is not an ff14SB-ready all-atom structure. "
                        "Choose automatic PDB2PQR/PROPKA preparation, or provide an AMBER-prepared PDB.\n{}".format(
                            summary
                        )
                    )
                if protein_preparation is not None:
                    raise ValueError(
                        "Automatic AMBER protein preparation completed, but exact ff14SB validation still failed. "
                        "Inspect {}.\n{}".format(protein_preparation.log_path, summary)
                    )
                raise ValueError(
                    "The prepared PDB does not match the selected all-atom topology.\n{}".format(
                        summary
                    )
                )
            if self.use_region.isChecked():
                if self._region_settings is None:
                    residues = self._available_residues()
                    self._region_settings = RegionSettings(residues, residues[:1], ())
                write_region_file(project / "region" / "region.data", self._region_settings)
            if not self.input_preview.toPlainText().strip():
                self._generate_preview()
            if not self.input_preview.toPlainText().strip():
                raise ValueError("mcmc.input is empty; generate or enter the input before preparing")
            input_path = write_mcmc_input(
                project / input_name, self.input_preview.toPlainText()
            )
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "PymoSAICS", "Project preparation failed: {}".format(exc))
            return None
        self._loaded_input_structure = prepared.pdb_path
        self.input_combo.setEditText(str(input_path))
        self._scan_inputs(preferred=input_path)
        try:
            cmd.load(str(prepared.pdb_path), "pymosaics_input")
            cmd.orient("pymosaics_input")
        except Exception as exc:
            self.build_status.setPlainText("Prepared files, but PyMOL could not load the input: {}".format(exc))
        else:
            details = "Prepared {} atoms in {} residues. Header: {}".format(
                prepared.atom_count, prepared.residue_count, prepared.header_line or "none"
            )
            details += "\nForce field: {}".format(staged_force_field)
            if protein_preparation is not None:
                details += (
                    "\nProtein preparation: PDB2PQR AMBER + PROPKA at pH {:.2f}"
                    "\nPrepared all-atom source: {}"
                    "\nPreparation log: {}"
                ).format(
                    protein_preparation.ph,
                    protein_preparation.pdb_path,
                    protein_preparation.log_path,
                )
            if prepared.mapping_path:
                details += "\nNaming map: {}".format(prepared.mapping_path)
            self.build_status.setPlainText(details)
        self._save_configuration()
        self._validate_current_project()
        return input_path

    def _scan_inputs(self, checked=False, preferred=None):
        project = self._project_directory()
        if project is None or not project.is_dir():
            return
        paths = sorted(list(project.glob("*.input")) + list(project.glob("*.inp")))
        current = str(preferred or self.input_combo.currentText())
        self.input_combo.clear()
        for path in paths:
            self.input_combo.addItem(str(path.resolve()))
        if current:
            self.input_combo.setEditText(current)

    def _browse_input(self, checked=False):
        start = str(self._project_directory() or Path.home())
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select MOSAICS input", start, "MOSAICS input (*.input *.inp *.txt);;All files (*)")
        if path:
            self.input_combo.setEditText(path)
            self.project_edit.setText(str(Path(path).parent))
            self._validate_current_project()

    def _browse_output(self, checked=False):
        start = str(self._project_directory() or Path.home())
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select output PDB", start, "PDB files (*.pdb);;All files (*)")
        if path:
            self.output_edit.setText(path)

    def _show_text_file(self, path: Path, editable: bool):
        try:
            if not path.expanduser().is_file():
                raise ValueError("file does not exist: {}".format(path))
            dialog = TextFileDialog(path, editable, self)
            dialog.exec()
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", str(exc))

    def _view_run_input(self, checked=False):
        text = self.input_combo.currentText().strip()
        if not text:
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "Select or prepare a MOSAICS input first.")
            return
        self._show_text_file(Path(text), editable=True)
        self._validate_current_project()

    def _view_prepared_structure(self, checked=False):
        project = self._project_directory()
        candidates = [self._loaded_input_structure]
        if project:
            candidates.append(project / "structure.pdb")
        path = next((item for item in candidates if item and item.is_file()), None)
        if path is None:
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "Prepare the project before editing its input PDB.")
            return
        self._show_text_file(path, editable=True)
        try:
            cmd.delete("pymosaics_input")
            cmd.load(str(path), "pymosaics_input")
            cmd.orient("pymosaics_input")
        except Exception as exc:
            self.log_output.appendPlainText("Edited PDB was saved, but PyMOL could not reload it: {}".format(exc))

    def _view_latest_log(self, checked=False):
        project = self._project_directory()
        path = latest_log(project) if project else None
        if path is None:
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "No MOSAICS run log was found in this project.")
            return
        self._show_text_file(path, editable=False)

    def _browse_executable(self, checked=False):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select MOSAICS executable")
        if path:
            self.custom_executable_edit.setText(path)
            self.runtime_combo.setCurrentIndex(self.runtime_combo.findData("custom"))

    def _browse_pdb2pqr(self, checked=False):
        start = self.pdb2pqr_edit.text().strip() or str(Path.home())
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select PDB2PQR executable", start
        )
        if path:
            self.pdb2pqr_edit.setText(path)
            self._update_protein_prep_controls()

    def _browse_workspace(self, checked=False):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select default workspace", str(Path.home()))
        if path:
            self.workspace_edit.setText(path)

    def _validate_current_project(self, checked=False):
        parameter_input = Path(self.input_combo.currentText().strip()).expanduser()
        config = self._configuration_from_fields()
        diagnostics = validate_project(parameter_input, config)
        self.validation_output.setPlainText(format_diagnostics(diagnostics))
        try:
            actual_input = planned_parameter_input(parameter_input, config)
        except (OSError, UnicodeError):
            actual_input = parameter_input
        self.command_preview.setPlainText(
            "Program: {}\nArgument 1: {}\nWorking directory: {}".format(
                config.executable.expanduser().resolve(), actual_input, parameter_input.expanduser().resolve().parent
            )
        )
        return diagnostics

    def _start_process(self, checked=False):
        if self.process.state() != PROCESS_NOT_RUNNING:
            return
        use_live_coordinates = self.source_mode.currentText() == "Live PyMOL object" or (
            self.follow_pymol_edits.isChecked() and bool(self._pymol_source_object())
        )
        if use_live_coordinates and self._project_directory():
            if self._prepare_project() is None:
                return
        diagnostics = self._validate_current_project()
        if has_errors(diagnostics):
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Correct the validation errors before running.")
            return
        try:
            active_run = prepare_run(Path(self.input_combo.currentText().strip()), self._configuration_from_fields())
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
        preamble = "Starting MOSAICS\nProgram: {}\nInput: {}\nWorking directory: {}\nLog: {}\n\n".format(
            active_run.command[0], active_run.command[1], active_run.working_directory, active_run.log_file
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
        self.log_output.moveCursor(TEXT_CURSOR_END)
        self.log_output.insertPlainText(data.decode("utf-8", errors="replace"))
        self.log_output.moveCursor(TEXT_CURSOR_END)

    def _process_error(self, error):
        message = "\nProcess error: {}".format(self.process.errorString())
        self.log_output.appendPlainText(message)
        if self._log_handle is not None:
            self._log_handle.write((message + "\n").encode("utf-8"))
            self._log_handle.flush()
        if self.process.state() == PROCESS_NOT_RUNNING:
            self._close_run_log()
            self._active_run = None

    def _close_run_log(self):
        if self._log_handle is not None:
            self._log_handle.close()
            self._log_handle = None
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def _process_finished(self, exit_code, exit_status):
        self._read_process_output()
        succeeded = exit_status == PROCESS_NORMAL_EXIT and exit_code == 0
        completion = "\nMOSAICS finished with exit code {} ({}).".format(
            exit_code, "success" if succeeded else "failure"
        )
        self.log_output.appendPlainText(completion)
        if self._log_handle is not None:
            self._log_handle.write((completion + "\n").encode("utf-8"))
            self._log_handle.flush()
        self._close_run_log()
        self._active_run = None
        if succeeded and self.auto_load.isChecked():
            self._load_all_outputs()
        self._refresh_analysis()

    def _stop_process(self, checked=False):
        if self.process.state() == PROCESS_NOT_RUNNING:
            return
        self.log_output.appendPlainText("\nStopping MOSAICS…")
        self.process.terminate()
        QtCore.QTimer.singleShot(3000, self._kill_if_running)

    def _kill_if_running(self):
        if self.process.state() != PROCESS_NOT_RUNNING:
            self.log_output.appendPlainText("MOSAICS did not stop; terminating it forcefully.")
            self.process.kill()

    def _load_pdb_paths(self, paths: Sequence[Path]):
        loaded = []
        for path in paths:
            try:
                name = _safe_name("pymosaics_" + path.stem)
                cmd.load(str(path), name)
                loaded.append("{} → {}".format(path.name, name))
            except Exception as exc:
                self.log_output.appendPlainText("Could not load {}: {}".format(path, exc))
        if loaded:
            self.log_output.appendPlainText("Loaded into PyMOL:\n" + "\n".join(loaded))

    def _load_all_outputs(self, checked=False):
        selected = self.output_edit.text().strip()
        if selected and Path(selected).is_file():
            paths = (Path(selected).resolve(),)
        else:
            project = self._project_directory()
            paths = discover_pdb_outputs(project, self._loaded_input_structure) if project else ()
        if not paths:
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "No output PDB or trajectory was found.")
            return
        self._load_pdb_paths(paths)

    def _refresh_analysis(self, checked=False):
        project = self._project_directory()
        if project is None or not project.is_dir():
            return
        self._energy_series = discover_energy_series(project)
        self.energy_combo.clear()
        current_landscape_energy = self.landscape_energy_combo.currentData()
        self.landscape_energy_combo.clear()
        self.landscape_energy_combo.addItem("Sampling order (no matched energy)", "")
        for series in self._energy_series:
            self.energy_combo.addItem(series.label, str(series.path))
            self.landscape_energy_combo.addItem(series.label, str(series.path))
        if current_landscape_energy:
            index = self.landscape_energy_combo.findData(current_landscape_energy)
            if index >= 0:
                self.landscape_energy_combo.setCurrentIndex(index)
        elif self.landscape_energy_combo.count() > 1:
            self.landscape_energy_combo.setCurrentIndex(1)
        self._energy_changed()
        log = latest_log(project)
        summaries = parse_acceptance_log(log) if log else ()
        self.acceptance_table.setRowCount(len(summaries))
        for row, summary in enumerate(summaries):
            values = (summary.chain, summary.accepted, summary.attempted, "{:.2%}".format(summary.ratio))
            for column, value in enumerate(values):
                self.acceptance_table.setItem(row, column, QtWidgets.QTableWidgetItem(str(value)))
        project_files = discover_project_files(project)
        self.output_list.clear()
        for project_file in project_files:
            try:
                relative = project_file.path.relative_to(project)
            except ValueError:
                relative = project_file.path
            item = QtWidgets.QListWidgetItem("{}  ·  {}".format(relative, project_file.kind))
            item.setData(USER_ROLE, str(project_file.path))
            item.setData(USER_ROLE + 1, project_file.text_readable)
            item.setData(USER_ROLE + 2, project_file.loadable_in_pymol)
            self.output_list.addItem(item)
        current_trajectory = self.trajectory_combo.currentData()
        self.trajectory_combo.clear()
        for path in discover_pdb_outputs(project, self._loaded_input_structure):
            self.trajectory_combo.addItem(path.name, str(path))
        if current_trajectory:
            index = self.trajectory_combo.findData(current_trajectory)
            if index >= 0:
                self.trajectory_combo.setCurrentIndex(index)

    def _energy_changed(self, *_args):
        index = self.energy_combo.currentIndex()
        if index < 0 or index >= len(self._energy_series):
            self.energy_plot.set_series(None)
            self.energy_summary.setText("No energy output found")
            return
        series = self._energy_series[index]
        self.energy_plot.set_series(series)
        self.energy_summary.setText(
            "{} samples · min {:.6g} · mean {:.6g} · max {:.6g} · {}".format(
                len(series.values), series.minimum, series.mean, series.maximum, series.path
            )
        )

    def _build_landscape(self, checked=False):
        value = self.trajectory_combo.currentData()
        if not value:
            QtWidgets.QMessageBox.information(
                self, "PymoSAICS", "Run a simulation or select a multi-model PDB trajectory first."
            )
            return
        path = Path(value)
        self.landscape_status.setText("Computing aligned pairwise RMSD and the two-dimensional projection…")
        QtWidgets.QApplication.processEvents()
        try:
            result = build_landscape(
                path,
                representatives=self.landscape_representatives.value(),
                maximum_frames=self.landscape_maximum_frames.value(),
            )
            energies = self._landscape_energies(result)
            project = self._project_directory()
            table = write_landscape_table(
                (project or path.parent) / "analysis" / (path.stem + ".structural_landscape.tsv"),
                result,
                energies,
            )
        except (OSError, RuntimeError, ValueError) as exc:
            self._landscape_result = None
            self.landscape_plot.set_result(None)
            self.landscape_status.setText("Structural map could not be built: {}".format(exc))
            return
        self._landscape_result = result
        self._landscape_path = path.resolve()
        self.landscape_plot.set_result(result, energies)
        self.representative_list.clear()
        for frame in result.representative_frames:
            item = QtWidgets.QListWidgetItem("Frame {}".format(frame))
            item.setData(USER_ROLE, frame)
            self.representative_list.addItem(item)
        coloring = "energy" if energies is not None else "sampling order"
        self.landscape_status.setText(
            "Mapped {} sampled frames using {} common atoms; points colored by {}. Coordinates: {}".format(
                len(result.frame_numbers), result.atom_count, coloring, table
            )
        )

    def _landscape_energies(self, result: LandscapeResult) -> Optional[Tuple[float, ...]]:
        selected = self.landscape_energy_combo.currentData()
        if not selected:
            return None
        series = next((item for item in self._energy_series if str(item.path) == selected), None)
        if series is None:
            return None
        if result.frame_numbers and max(result.frame_numbers) <= len(series.values):
            return tuple(series.values[frame - 1] for frame in result.frame_numbers)
        if len(series.values) == len(result.frame_numbers):
            return tuple(series.values)
        return None

    def _ensure_landscape_in_pymol(self) -> Optional[str]:
        if self._landscape_path is None or not self._landscape_path.is_file():
            return None
        name = _safe_name("pymosaics_landscape_" + self._landscape_path.stem)
        try:
            objects = set(cmd.get_object_list("all"))
            if name not in objects:
                cmd.load(str(self._landscape_path), name)
                cmd.show("cartoon", name)
                cmd.orient(name)
            self._landscape_object = name
            return name
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Cannot load trajectory into PyMOL: {}".format(exc))
            return None

    def _landscape_frame_selected(self, frame: int):
        name = self._ensure_landscape_in_pymol()
        if not name:
            return
        try:
            states = int(cmd.count_states(name))
            if frame < 1 or frame > states:
                raise ValueError("trajectory has {} PyMOL states, not frame {}".format(states, frame))
            cmd.frame(frame)
            cmd.enable(name)
            self.landscape_status.setText(
                "Showing frame {} of {} in PyMOL object {}. The map is an RMSD projection.".format(
                    frame, states, name
                )
            )
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "PymoSAICS", "Cannot show frame {}: {}".format(frame, exc))

    def _representative_selected(self, item):
        self._landscape_frame_selected(int(item.data(USER_ROLE)))

    def _load_selected_analysis_outputs(self, checked=False):
        paths = [
            Path(item.data(USER_ROLE))
            for item in self.output_list.selectedItems()
            if bool(item.data(USER_ROLE + 2))
        ]
        if not paths:
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "Select one or more PDB files first.")
            return
        self._load_pdb_paths(paths)

    def _view_selected_analysis_output(self, checked=False):
        items = self.output_list.selectedItems()
        if len(items) != 1:
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "Select one text-readable project file.")
            return
        item = items[0]
        if not bool(item.data(USER_ROLE + 1)):
            QtWidgets.QMessageBox.information(self, "PymoSAICS", "The selected output is binary and cannot be shown as text.")
            return
        self._show_text_file(Path(item.data(USER_ROLE)), editable=False)

    def _open_analysis_item(self, item):
        path = Path(item.data(USER_ROLE))
        if bool(item.data(USER_ROLE + 2)):
            self._load_pdb_paths((path,))
        elif bool(item.data(USER_ROLE + 1)):
            self._show_text_file(path, editable=False)

    def closeEvent(self, event):
        if self.process.state() != PROCESS_NOT_RUNNING:
            answer = QtWidgets.QMessageBox.question(
                self,
                "PymoSAICS",
                "MOSAICS is still running. Stop it and close PymoSAICS?",
                MESSAGE_YES | MESSAGE_NO,
                MESSAGE_NO,
            )
            if answer != MESSAGE_YES:
                event.ignore()
                return
            self.process.kill()
            self.process.waitForFinished(3000)
        self._close_run_log()
        event.accept()
