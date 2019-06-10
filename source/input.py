import re

class Input(object):
    def __init__(self, geom, tmp_dir):
        self.geom = geom
        self.tmp_dir = tmp_dir

    def subst_template(self, inp_fname, out_fname):
        #open file and read content
        with open(inp_fname,'r') as file:
            content = file.read().strip()

        #geometry section pattern
        geometry_pattern = r'[^\S\r\n]*\$geometry(?:|\s+)\{\n([\s\S]*?)\};$'
        geometry_regex = re.compile(geometry_pattern,re.IGNORECASE|re.MULTILINE)

        #replace geometry section by template specification
        for match in geometry_regex.finditer(content):
            geometry_template_section = match.group(0)
            geometry_template_str = match.group(1)
            geometry_str = self.subst_geometry(geometry_template_str)

            #replace geometry section in content
            content = content.replace(geometry_template_section,geometry_str)

        #replace %TMP_DIR variable
        tmp_dir_pattern = "%TMP_DIR"
        content = self.subst_string_pattern(content,tmp_dir_pattern,self.tmp_dir)

        #write content to new output file
        with open(out_fname,"w") as file:
            file.write(content)

    def subst_geometry(self, string):
        res = ''
        n_atoms = self.geom.get_n_atoms()
        for i_atom in range(n_atoms):
            x = self.geom.get_i_coord(3*i_atom)
            y = self.geom.get_i_coord(3*i_atom+1)
            z = self.geom.get_i_coord(3*i_atom+2)

            subst = string
            subst = self.subst_string_pattern(subst,"%X",str(x))
            subst = self.subst_string_pattern(subst,"%Y",str(y))
            subst = self.subst_string_pattern(subst,"%Z",str(z))

            AL = self.geom.atoms[i_atom]
            subst = self.subst_string_pattern(subst,"%AL",AL)

            AM = self.geom.atomic_masses[i_atom]
            subst = self.subst_string_pattern(subst,"%AM",str(AM))

            AN = self.geom.atomic_numbers[i_atom]
            subst = self.subst_string_pattern(subst,"%AN",str(AN))

            res += subst

        return res

    def subst_string_pattern(self, string, pattern, subst):
        res = re.sub(pattern,subst,string)
        return res