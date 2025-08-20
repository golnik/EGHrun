import re

class Reader(object):
    #def __init__(self):
    def get_gradient(self, string):
        pass

    def get_energy(self, string):
        energy_pattern = r'[^\S\r\n]*\$energy(?:\_)?(\w+)?\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)$'
        energy_regex = re.compile(energy_pattern,re.IGNORECASE|re.MULTILINE)

        energies = []
        states   = []
        for match in re.finditer(energy_regex, string):
            if len(match.groups()) == 2:
                energy = float(match.group(1))
            elif len(match.groups()) == 3:
                state  = match.group(1)
                energy = float(match.group(2))
                states.append(state)

            energies.append(energy)

        return energies, states