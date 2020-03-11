
class Modification(object):
    def __init__(self, location, monoisotopicMassDelta, name):
        self.location = location  # 0 = at the beginning; 1 = first amino acid; modification to c terminal end = N+1
        self.name = name
        self.monoisotopicMassDelta = monoisotopicMassDelta


class Modification_on_Residue(Modification):
    def __init__(self, location, monoisotopicMassDelta, name, residues):
        self.residues = residues  # if multiple residues = can not assign modification with 100%
        super().__init__(location, monoisotopicMassDelta, name)
