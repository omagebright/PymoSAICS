import unittest

from pymosaics.core.nucleic_builder import (
    hydrogen_names_for_parent,
    plan_nucleic_acid_build,
    reverse_complement,
)


class NucleicAcidBuilderTests(unittest.TestCase):
    def test_force_field_hydrogen_names_cover_dna_rna_and_terminals(self):
        self.assertEqual(hydrogen_names_for_parent("DA", "C2'"), ("H2'", "H2''"))
        self.assertEqual(hydrogen_names_for_parent("A", "C2'"), ("H2'",))
        self.assertEqual(hydrogen_names_for_parent("A", "O2'"), ("HO2'",))
        self.assertEqual(
            hydrogen_names_for_parent("DT", "C5M"), ("H51", "H52", "H53")
        )
        self.assertEqual(hydrogen_names_for_parent("DG", "N2"), ("H21", "H22"))
        self.assertEqual(hydrogen_names_for_parent("DC", "O2P"), ())
        self.assertEqual(hydrogen_names_for_parent("A", "O3'", False), ())

    def test_hybrid_plan_builds_an_antiparallel_rna_complement(self):
        plan = plan_nucleic_acid_build(
            "dna-rna-hybrid", "ACG", strand1_form="A", strand2_form="A"
        )
        self.assertEqual(plan.strands[0].polymer, "DNA")
        self.assertEqual(plan.strands[1].polymer, "RNA")
        self.assertEqual(plan.strands[1].sequence, "CGU")
        self.assertEqual(plan.strands[1].chain, "B")

    def test_duplex_complement_and_noncanonical_rna_b_warning(self):
        self.assertEqual(reverse_complement("AUTG", "DNA"), "CAAT")
        plan = plan_nucleic_acid_build(
            "rna-duplex", "ACGU", strand1_form="B", strand2_form="B"
        )
        self.assertTrue(any("noncanonical" in warning.lower() for warning in plan.warnings))

    def test_invalid_sequence_and_mismatched_custom_complement_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "unsupported"):
            plan_nucleic_acid_build("single-dna", "ACX", strand1_form="B")
        with self.assertRaisesRegex(ValueError, "reverse complement"):
            plan_nucleic_acid_build(
                "dna-duplex", "ACG", strand1_form="B", strand2_form="B", sequence2="AAA"
            )


if __name__ == "__main__":
    unittest.main()
