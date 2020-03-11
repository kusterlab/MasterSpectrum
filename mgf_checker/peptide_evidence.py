import math
from mgf_checker.modification import Modification
from mgf_checker.modification import Modification_on_Residue


class PeptideEvidence(object):
    def __init__(self, isDecoy, start, end, peptide_ref, id, pre, dBSequence_ref, post):
        self.isDecoy = isDecoy  # bool
        self.start = start  # as_position
        self.end = end  # as_position for end
        self.peptide_ref = peptide_ref  # peptide_ref
        self.id = id  # protein id
        self.pre = pre  # AS before cut
        self.dBSequence_ref = dBSequence_ref  # database ref
        self.post = post  # amino acid after cut
        self.modification = {}
        self.peptide_sequence = ''

    def set_peptide_modification(self, data):
        self.peptide_sequence = data['PeptideSequence']
        if 'Modification' in data:
            for mod in data['Modification']:
                if 'residues' in mod:  # is a modification connected to a amino acid
                    self.modification[mod['location']] = Modification_on_Residue(location=mod['location'],
                                                                                 monoisotopicMassDelta=mod['monoisotopicMassDelta'],
                                                                                 name=mod['name'],
                                                                                 residues=mod['residues'])
                else:  # is a fixed modification
                    self.modification[mod['location']] = Modification(location=mod['location'],
                                                                      monoisotopicMassDelta=mod['monoisotopicMassDelta'],
                                                                      name=mod['name'])

    def get_annotated_positions(self):
        """
        my algorithm works just for TMT
        returns:
        b_tmt = hashset with amount of tag at position (position starts at 1)
        y_tmt = hashset with amount of tag at position (position starts at 1)
        """
        if len(self.modification) != 0:
            len_peptide = len(self.peptide_sequence)  # N
            tmt_position = []  # starting at N-Terminus; first amino acid is 1
            for mod in self.modification:
                if self.modification[mod].name == 'TMT6plex':
                    tmt_position.append(self.modification[mod].location)

            tmt_amount = 0

            if 0 in tmt_position:  # tmt amount at start
                tmt_amount += 1

            b_tmt = {}
            for i in range(1, len_peptide + 1):
                if i in tmt_position:
                    tmt_amount += 1
                b_tmt[i] = tmt_amount

            tmt_position_y = [min(int(math.fabs(i - len_peptide - 1)), len_peptide) for i in tmt_position]

            tmt_amount = 0
            y_tmt = {}
            for i in range(1, len_peptide + 1):
                if i in tmt_position_y:
                    tmt_amount += 1
                y_tmt[i] = tmt_amount

            return b_tmt, y_tmt
        else:
            return {}, {}
