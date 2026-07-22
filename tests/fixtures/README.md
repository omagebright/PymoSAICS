# Runtime regression fixtures

`single_gud_cblc.pdb` is the canonical deoxyguanosine structure produced by
the PymoSAICS single-DNA builder and normalized by its MOSAICS preparation
step. The `GUD` residue and `CBLC >A` header exercise both terminal backbone
closure and the purine glycosidic chi move. It is used only for deterministic,
short runtime regression tests.
