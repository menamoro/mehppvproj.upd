
def gen_lj_dict(path):
    g = open(path+"OPLS-AA-lj-dihedral.txt", "r")
    oplsaa_list = [y.split() for y in g]

    oplsaa_dict ={}

    for k in oplsaa_list:
        if len(k) > 1:
            oplsaa_dict[(k[0])] = k[1:]
    g.close()
    return oplsaa_dict

def gen_harm_bond_dict(path):
    f = open(path+"OPLS-AA-harmbond.txt", "r")
    harm_bond_dict = {}
    for op12 in f:
        harm_bond_dict[op12[:5]] = op12[5:].split()
    f.close()

    return harm_bond_dict

def gen_angle_dict(path):
    f = open(path+"OPLS-AA-bondangles.txt", "r")
    angle_dict= {}
    for xe in f:
        rest_of_list = xe[8:].split()
        angle_dict[xe[:8]] = rest_of_list

    f.close()

    return angle_dict

def gen_dihedral_dict(path):
    f = open(path+"OPLS-AA-dihedral.txt", "r")
    dihedral_dict = {}
    for op19 in f:
        key = op19[47:58]
        str_of_concern = op19[3:39].split()
        dihedral_dict[key] = str_of_concern

    f.close()
    return dihedral_dict
