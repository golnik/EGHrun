import re

class Reader(object):
    #def __init__(self):
    def get_gradient(self, string):
        pass

    def get_energy(self, string):
        energy_pattern = r'[^\S\r\n]*\$energy\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)$'
        energy_regex = re.compile(energy_pattern,re.IGNORECASE|re.MULTILINE)

        match = re.search(energy_regex,string)

        energy = float(match.group(1))

        return energy