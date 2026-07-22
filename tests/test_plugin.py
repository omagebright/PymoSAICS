import sys
import types
import unittest
from unittest.mock import patch

import pymosaics


class PluginRegistrationTests(unittest.TestCase):
    def test_modern_plugin_entrypoint_registers_qt_menu(self):
        labels = []
        pymol_module = types.ModuleType("pymol")
        plugins_module = types.ModuleType("pymol.plugins")
        plugins_module.addmenuitemqt = lambda label, callback: labels.append((label, callback))
        pymol_module.plugins = plugins_module

        with patch.dict(sys.modules, {"pymol": pymol_module, "pymol.plugins": plugins_module}):
            pymosaics.__init_plugin__()

        self.assertEqual(labels[0][0], "PymoSAICS")
        self.assertTrue(callable(labels[0][1]))


if __name__ == "__main__":
    unittest.main()
