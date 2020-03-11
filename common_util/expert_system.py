import json
from pyteomics import mass


class Expert_system(object):
    def __init__(self):
        self.data = {}
        self.masses = []
        self.load_data()
        self.parse_data()

    def load_data(self):
        with open('common_util/rules.json') as f:
            self.data = json.load(f)

    def parse_data(self):
        for key in self.data:
            for rules in self.data[key]:
                if rules["annotation"] == "standard":
                    for masses in rules["losses"]:
                        self.masses.append(mass.calculate_mass(formula=masses))
        self.masses = list(set(self.masses))
