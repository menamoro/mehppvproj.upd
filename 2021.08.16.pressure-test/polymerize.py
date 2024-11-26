import espressomd
from espressomd import interactions
from espressomd.io.writer import vtf
from espressomd.interactions import HarmonicBond
from types import *
import ffparams.generate_lists as generate_lists
import math
import numpy as np
"""
running this module requires you to set the put any changes in atomtypes in the .ids file, along with the !Pairs (pairs of points in monomer that will polymerize), !Leave (possible leaving groups for each point), and !Change (change in atomtypes when polymerizing) lists. other than that, this will parametrize the monomers and the oligomers as they form. no parametrization occurs for planar dihedrals, may be something to do in the future, if data becomes available. no preference for conformations (ie, cis/trans)

see packages for needed pacakages. needed rin yung opls-aa files and the script "generate_lists". you may need to edit the opls-aa files (dihedrals in particular)
things to rewrite:
    the garbage code you made to create bond_dict
    the implementation of "!Change" in condense()
things to add:
changing atomids depending on the polymerization target
    implementation: .ids file should be:
        !Change
        [7,3]
        7 id1
        3 id2
        6 id3
        [3,5]
if the added species if different from the base species
CIS/TRANS RNGs

pwede rin options for addition polymerization siguro
"""
#initialize all the globals. always call this first.
def init(system_inp, no_monomers_inp,script):
    global system, no_monomers, pairs,leave,change,harm_bond_dict,lj_dict,angle_dict,dihedral_dict,total_poly_atoms,adj_list,atom_type_dict, unique_bond_type_dict,unique_bond_angle_dict, unique_dihedral_dict,atom_type_str_dict,reverse_atom_type_dict,change_new,pairs_new,leave_new
    harm_bond_dict, lj_dict, angle_dict,dihedral_dict = generate_lists.gen_harm_bond_dict("ffparams/"), generate_lists.gen_lj_dict("ffparams/"), generate_lists.gen_angle_dict("ffparams/"), generate_lists.gen_dihedral_dict("ffparams/")
    system = system_inp
    no_monomers = no_monomers_inp
    exec(open(f"{script}-params.py").read(),globals())

    #read the .ids file to read the pairs being bonded, the leaving groups (this is a condensation polymerization), and the change in atom types.
    pairs = []
    leave = {}
    change = {}
    with open(f"{script}.ids","r") as f:
        curr_str =""
        for line in f:
            if line == "!Pairs\n" or line == "!Leave\n" or line == "!Change\n":
                curr_str = line
                continue
            if curr_str == "!Pairs\n":
                temp_list = [int(x) for x in line.split()]
                pairs.append(temp_list)
            elif curr_str == "!Leave\n":
                temp_list = [int(x) for x in line.split()]
                leave[temp_list[0]] = temp_list[1:]
            elif curr_str == "!Change\n":
                temp_list = [int(x) for x in line.split()]
                change[temp_list[0]] = temp_list[1]
            if line[:3] == "!No": no_atom_per_monomer = int(line.split()[1])
            if line[:3] == "!Na": name = line.split()[1]

    total_poly_atoms = no_monomers* no_atom_per_monomer

    adj_list = {}
    create_adj_list()
    #read the .dict file to get the dictionaries for already declared bonds.atom_type_dict is for OPLS-AA type to system type.
    atom_type_dict, unique_bond_type_dict,unique_bond_angle_dict, unique_dihedral_dict = {},{},{}, {}
    with open(f"{script}.dict","r") as f:
        charlength = 2
        for line in f:
            if line[0] == "!":
                charlength += 3
                continue
            if charlength == 5:
                unique_bond_type_dict[line[:5]] = int(line[6:])
            elif charlength == 8:
                unique_bond_angle_dict[line[:8]] = int(line[9:])
            elif charlength == 11:
                unique_dihedral_dict[line[:11]] = line[12:].split()
            else:
                curr_line = line.split()
                atom_type_dict[int(curr_line[0])] = int(curr_line[1])
    #add atomtypes that are not currently present in the system
    max_type = np.amax(list(atom_type_dict.values()))
    max_type += 1
    for op1 in change.values():
        if op1 not in atom_type_dict:
            atom_type_dict[op1] = max_type
            max_type += 1
            for op2 in atom_type_dict:
                lj_1 = lj_dict[str(op2)]
                lj_2 = lj_dict[str(op1)]
                ep_in_e = (float(lj_1[4])*float(lj_2[4]))**(1/2) * 4.1858 * 100
                sigma = (float(lj_1[3]) + float(lj_2[3]) )/ 2
                exec(f"system.non_bonded_inter[{atom_type_dict[op2]},{atom_type_dict[op1]}].lennard_jones.set_params(epsilon={ep_in_e}, sigma = {sigma},cutoff = 12.000, shift=\"auto\")")


    #atom_type_str_dict converts atomtype to string type (ie, 1 -> "H ") and reverse_atom_type_dict converts system type to OPLS-AA type.
    atom_type_str_dict = {}
    reverse_atom_type_dict = {}
    for op in atom_type_dict:
        str_temp = lj_dict[str(op)][1]
        if len(str_temp) == 1: str_temp = str_temp + " "
        atom_type_str_dict[op] = str_temp
        reverse_atom_type_dict[atom_type_dict[op]] = op

    #create the list of changes for all monomers. done because i didnt want to track by molecules
    change_new = {}
    for op1 in change:
        for op2 in range(no_monomers):
            change_new[op1 + (op2*no_atom_per_monomer)] = change[op1]
    #create list of pairs to monitor
    pairs_new = []
    for op1 in pairs:
        temp_list = []
        for loop in range(no_monomers):
            temp_list.append([op1[0]+(loop*no_atom_per_monomer), op1[1] + (loop*no_atom_per_monomer)])
        pairs_new.append(temp_list)
    #create list of leaving groups to monitor
    leave_new = {}
    for op1 in leave:
        for op2 in range(no_monomers):
            temp_list = []
            for op3 in leave[op1]:
                temp_list.append(op3 + (op2*no_atom_per_monomer))
            leave_new[op1 + (op2*no_atom_per_monomer)] = temp_list


#create an adjacency list for all atoms
def create_adj_list():
    for op1 in range(total_poly_atoms): adj_list[op1] = []

    for op1 in range(total_poly_atoms):
        bonds = system.part[op1].bonds
        for op2 in bonds:
            if isinstance(op2[0],espressomd.interactions.HarmonicBond):
                adj_list[op1].append(op2[1])
                adj_list[op2[1]].append(op1)

#updates the adj list every time a bond is broken or made
def update_adj_list(c1,c2,didBreak):
    if didBreak:
        adj_list[c1].remove(c2)
        adj_list[c2].remove(c1)
    else:
        adj_list[c1].append(c2)
        adj_list[c2].append(c1)



#inp = particle id, link = list to edit, visited = a list of -1. run DFS to append to link a list of connected nodes.
def dfs_link(inp,link,visited):
    visited[inp] = 1
    link.append(inp)
    for op1 in adj_list[inp]:
        if visited[op1] == -1:
            dfs_link(op1,link,visited)

#delete all non-harmonicbonds in the unit that atom id c1 belongs to.
def del_nonharm_bonds_in_unit(c1):
    visited = [-1]*(total_poly_atoms)
    to_del = []
    dfs_link(c1,to_del,visited)
    for op1 in to_del:
        bonds = system.part[op1].bonds
        for op2 in bonds:
            if isinstance(op2[0],espressomd.interactions.HarmonicBond): continue
            rest_of_list = op2[1:]
            system.part[op1].delete_bond((op2[0],*rest_of_list))

#delete harmonic bond between a and b.
def del_harm(a,b):
    for op1 in system.part[a].bonds:
        if isinstance(op1[0],espressomd.interactions.HarmonicBond) and op1[1] == b:
            system.part[a].delete_bond((op1[0],op1[1]))
    for op1 in system.part[b].bonds:
        if isinstance(op1[0],espressomd.interactions.HarmonicBond) and op1[1] == a:
            system.part[a].delete_bond((op1[0],op1[1]))
    update_adj_list(a,b,True)


#create a harmonic bond between a and b. if bond type doesnt exist, create one.
def cre_harm(a,b):
    a_type, b_type = reverse_atom_type_dict[system.part[a].type], reverse_atom_type_dict[system.part[b].type]

    str_bond1 = atom_type_str_dict[a_type] + "-" + atom_type_str_dict[b_type]
    str_bond2 = atom_type_str_dict[b_type] + "-" + atom_type_str_dict[a_type]

    harm_line = -1
    max_hb = len(unique_bond_type_dict)
    new_bond, new_bond2 = False,False
    if str_bond1 in unique_bond_type_dict:
        harm_line = unique_bond_type_dict[str_bond1]
        exec(f"system.part[{a}].add_bond((hb{harm_line},{b}))",globals())
    elif str_bond2 in unique_bond_type_dict:
        harm_line = unique_bond_type_dict[str_bond2]
        exec(f"system.part[{a}].add_bond((hb{harm_line},{b}))",globals())
    elif str_bond1 in harm_bond_dict:
        new_bond = True
        harm_line = harm_bond_dict[str_bond1]
    elif str_bond2 in harm_bond_dict:
        new_bond2 = True
        harm_line = harm_bond_dict[str_bond2]
    if new_bond or new_bond2:
        force_k = float(harm_line[0]) * 4.1858 * 200
        eq_dist = float(harm_line[1])
        exec(f"hb{max_hb} = espressomd.interactions.HarmonicBond(k = {force_k}, r_0 = {eq_dist})", globals())
        exec(f"system.bonded_inter.add(hb{max_hb})")
        exec(f"system.part[{a}].add_bond((hb{max_hb},{b}))")
        if new_bond:
            unique_bond_type_dict[str_bond1] = max_hb
        if new_bond2:
            unique_bond_type_dict[str_bond2] = max_hb
    update_adj_list(a,b,False)

#re-parametrize the atoms in the network that atom_in_network belongs to
def cre_harm_param(atom_in_network):
    dfs_atoms = []
    visited = [-1]*total_poly_atoms
    dfs_link(atom_in_network,dfs_atoms, visited)
    for op1 in dfs_atoms:
        for op2 in system.part[op1].bonds:
            if isinstance(op2[0], espressomd.interactions.HarmonicBond):
                del_harm(op1,op2[1])
                cre_harm(op1,op2[1])


#create the angle and dihedral parameters in the network that atom_in_network belongs to
def cre_angle_and_dihedral(atom_in_network):
    temp = []
    visited = [-1]*(total_poly_atoms)
    dfs_link(atom_in_network,temp, visited)
    bonds_dict = cre_bond_dict(temp)
    for op1 in bonds_dict:
        curr_line = bonds_dict[op1]
        if len(curr_line) < 2: continue
        for op2 in range(len(curr_line)):
            for op3 in range(op2+1,len(curr_line)):
                type1, type2, type3 = reverse_atom_type_dict[system.part[curr_line[op2]].type],  reverse_atom_type_dict[system.part[op1].type], reverse_atom_type_dict[system.part[curr_line[op3]].type]
                atom1,atom2,atom3 = atom_type_str_dict[type1],atom_type_str_dict[type2],atom_type_str_dict[type3]
                str_temp = atom1 + "-" + atom2 + "-" + atom3
                str_temp2 =  atom3 + "-" + atom2 + "-" + atom1
                new_angle,new_angle2 = False, False
                angle_line = []
                max_ah = len(unique_bond_angle_dict)
                if str_temp in unique_bond_angle_dict:
                    exec(f"system.part[{op1}].add_bond((ah{unique_bond_angle_dict[str_temp]},{curr_line[op2]},{curr_line[op3]}))")
                elif str_temp2 in unique_bond_angle_dict:
                    exec(f"system.part[{op1}].add_bond((ah{unique_bond_angle_dict[str_temp2]},{curr_line[op3]},{curr_line[op2]}))")
                elif str_temp in angle_dict:
                    angle_line = angle_dict[str_temp]
                    new_angle = True
                elif str_temp2 in angle_dict:
                    angle_line = angle_dict[str_temp2]
                    new_angle2 = True
                if new_angle or new_angle2:
                    force_k = float(angle_line[0]) * 4.1858 * 200
                    eq_angle = float(angle_line[1]) * math.pi/180
                    exec(f"ah{max_ah} = espressomd.interactions.AngleHarmonic(bend = {force_k}, phi0 = {eq_angle})", globals())
                    exec(f"system.bonded_inter.add(ah{max_ah})")
                    if new_angle:
                        unique_bond_angle_dict[str_temp] = max_ah
                        exec(f"system.part[{op1}].add_bond((ah{max_ah},{curr_line[op2]},{curr_line[op3]}))")
                    elif new_angle2:
                        unique_bond_angle_dict[str_temp2] = max_ah
                        exec(f"system.part[{op1}].add_bond((ah{max_ah},{curr_line[op3]},{curr_line[op2]}))")


    max_da = -1
    for count in unique_dihedral_dict.values():
        for count2 in count:
            if int(count2) > max_da: max_da = int(count2)
    max_da+=1
    for op1 in bonds_dict:
        curr_line = bonds_dict[op1]
        if len(curr_line) < 2: continue
        for op2 in range(len(curr_line)):
            for op3 in range(len(curr_line)):
                if op2 == op3: continue
                if curr_line[op3] not in bonds_dict:continue
                for op4 in bonds_dict[curr_line[op3]]:
                    if op4 == op2: continue
                    if op4 == op1: continue
                    type1, type2, type3,type4 = reverse_atom_type_dict[system.part[curr_line[op2]].type],  reverse_atom_type_dict[system.part[op1].type], reverse_atom_type_dict[system.part[curr_line[op3]].type],reverse_atom_type_dict[system.part[op4].type]

                    atom1,atom2,atom3,atom4 = atom_type_str_dict[type1],atom_type_str_dict[type2],atom_type_str_dict[type3], atom_type_str_dict[type4]
                    str_temp = atom1 + "-" + atom2 + "-" + atom3 +  "-" + atom4
                    str_temp2 =  atom4 + "-" + atom3 + "-" + atom2 + "-" + atom1

                    new_dihedral, new_dihedral2 = False, False
                    dihedral_line =[]
                    dihedral_numbers = []

                    if str_temp in unique_dihedral_dict:
                        for op5 in unique_dihedral_dict[str_temp]:
                            exec(f"system.part[{op1}].add_bond((da{op5},{curr_line[op2]},{curr_line[op3]},{op4}))")
                    elif str_temp2 in unique_dihedral_dict:
                        for op5 in unique_dihedral_dict[str_temp2]:
                            exec(f"system.part[{curr_line[op3]}].add_bond((da{(op5)},{op4},{curr_line[op2]},{op1}))")
                    elif str_temp in dihedral_dict:
                        dihedral_line = dihedral_dict[str_temp]
                        new_dihedral = True
                    elif str_temp2 in dihedral_dict:
                        dihedral_line = dihedral_dict[str_temp2]
                        new_dihedral2 = True
                    if new_dihedral or new_dihedral2:
                        new_append = []
                        for op5 in range(len(dihedral_line)):
                            if float(dihedral_line[op5]) == 0: continue
                            phase = math.pi
                            if (op5+1)%2 == 0: phase = 0
                            force_k  = float(dihedral_line[op5]) * 4.1858 * 50
                            exec(f"da{max_da} = espressomd.interactions.Dihedral(bend = {force_k}, mult = {op5+1},phase = {phase})", globals())
                            exec(f"system.bonded_inter.add(da{max_da})")
                            if new_dihedral:
                                exec(f"system.part[{op1}].add_bond((da{max_da},{curr_line[op2]},{curr_line[op3]},{op4}))")
                            elif new_dihedral2:
                                exec(f"system.part[{curr_line[op3]}].add_bond((da{(max_da)},{op4},{curr_line[op2]},{op1}))")
                            new_append.append(max_da)
                            max_da +=1
                        if new_dihedral:
                            unique_dihedral_dict[str_temp] = new_append
                        elif new_dihedral2:
                            unique_dihedral_dict[str_temp2] = new_append


#generate an adjacency list from a list of atoms. true adjacency list this time
def cre_bond_dict(dfs_atoms):
    bonds_dict = {}
    for op1 in dfs_atoms:
        bonds_dict[op1] = []
    for op1 in range(len(dfs_atoms)):
        for op2 in system.part[dfs_atoms[op1]].bonds:
            if isinstance(op2[0],espressomd.interactions.HarmonicBond):
                if op2[1] not in bonds_dict[dfs_atoms[op1]]:
                    bonds_dict[dfs_atoms[op1]].append(op2[1])
                if dfs_atoms[op1] not in bonds_dict[op2[1]]:
                    bonds_dict[op2[1]].append(dfs_atoms[op1])

    return bonds_dict
#the ctr-edge limit dfs
def dfs_exclude(atom_in_network,ctr,visited,list):
    if ctr < 0: return
    if visited[atom_in_network] == 1: return
    visited[atom_in_network] = 1
    list.append(atom_in_network)
    for op1 in adj_list[atom_in_network]:
        dfs_exclude(op1,ctr-1,visited,list)

#create exclude using a 3-edge limit DFS
def cre_exclude(atom_in_network):
    network = []
    visited = [-1]*(total_poly_atoms)
    dfs_link(atom_in_network,network,visited)
    for op1 in network:
        exclude = []
        visited = [-1]*(total_poly_atoms)
        dfs_exclude(op1, 3, visited, exclude)
        list_exclusions = system.part[op1].exclusions
        exclude[:] = [x for x in exclude if x != op1]
        for op2 in exclude:
            list_exclusions2 = system.part[op2].exclusions
            if op1 not in list_exclusions2 and op2 not in list_exclusions:
                exec(f"system.part[{op1}].add_exclusion({op2})")

#breaks C1-H1 and C2-H2, forms C1-C2 and H1-H2, assigns parameters for all formed molecules.
def condense(c1,c2,h1,h2):
    del_harm(c1,h1)
    del_harm(c2,h2)
    # replace this with a loop
    system.part[h1].type = atom_type_dict[change_new[h1]]
    system.part[c1].type = atom_type_dict[change_new[c1]]
    system.part[c2].type = atom_type_dict[change_new[c2]]
    system.part[h2].type = atom_type_dict[change_new[h2]]
    system.part[h1].q = float(lj_dict[str(change_new[h1])][2])
    system.part[h2].q = float(lj_dict[str(change_new[h2])][2])
    system.part[c1].q = float(lj_dict[str(change_new[c1])][2])
    system.part[c2].q = float(lj_dict[str(change_new[c2])][2])
    cre_harm(h1,h2)
    cre_harm(c1,c2)
    cre_harm_param(c1)
    cre_harm_param(h1)
    del_nonharm_bonds_in_unit(c1)
    del_nonharm_bonds_in_unit(h1)
    cre_angle_and_dihedral(c1)
    cre_angle_and_dihedral(h1)
    cre_exclude(c1)
    cre_exclude(h1)

def minimize_energy(system,max_f, max_d, gam, max_std_step, time_per_interval, f_tol,e_tol):
    system.thermostat.turn_off()
    system.integrator.set_steepest_descent(f_max  = max_f, max_displacement = max_d, gamma =  gam)

    ctr = 0
    system.integrator.run(0)
    f_prev = np.max(system.part[:].f)
    curr_f = 0
    e_prev = system.analysis.energy()["total"]
    ctr_energy = 0
    for i in range(int(max_std_step/time_per_interval)):
        system.integrator.run(time_per_interval)
        curr_f = np.max(system.part[:].f)
        e_curr = system.analysis.energy()["total"]
        rel_f = abs((curr_f-f_prev)/(f_prev))
        rel_e = abs((e_curr-e_prev)/(e_prev))
        if rel_f < f_tol:
            ctr += 1
        else:
            ctr = 0
        if rel_e < e_tol:
            ctr_energy += 1
        else:
            ctr_energy = 0
        if ctr_energy > 10 and ctr > 10:
            break
        f_prev = curr_f
        e_prev = e_curr
#check, pairwise, all atoms in  pairs_new.
def check_pairs(dist_req,angle_req):
    did_bond = False
    for op1 in pairs_new:
        for op2 in range(len(op1)):
            for op3 in range(len(op1)):
                is_bonded = False
                if op3 == op2: continue
                atom1 = op1[op2][0]
                atom2 = op1[op3][1]
                if atom1 == -1 or atom2 == -1: continue
                for op4 in leave_new[atom1]:
                    for op5 in leave_new[atom2]:
                        dist = get_distance(op4,op5)
                        angle = abs(get_angle([atom1,op4],[atom2,op5]) -math.pi)
                        if dist < dist_req and angle < angle_req:
                            condense(atom1,atom2,op4,op5)
                            print(f"found appropriate polymerization at {atom1},{op4} and {atom2},{op5}")
                            op1[op2][0] = -1
                            op1[op3][1] = -1
                            print(f"minimizing", end="\r")
                            minimize_energy(system,0,0.006,0.01,10000,5,1e-3,1e-3)
                            print(f"done with minimizing", end="\n")
                            is_bonded = True
                            did_bond = True
                            break
                    if is_bonded: break
    return did_bond

#gets distance between 2 atoms
def get_distance(a,b):
    dist = np.linalg.norm(system.part[a].pos - system.part[b].pos)
    return dist

#gets angle between vector in set1 and vector in set2
def get_angle(set1,set2):
    a1,a2 = set1
    b1,b2 = set2
    vec1 = system.part[a1].pos - system.part[a2].pos
    vec2 = system.part[b1].pos - system.part[b2].pos
    angle1 = np.arccos(-abs(np.dot(vec1,vec2))/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
    return angle1
