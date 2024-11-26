from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
import time
import sys



#from LigParGen, to properly implement OPLS-AA
def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce('4*epsilon*((sigma/r)^12-(sigma/r)^6);sigma=sqrt(sigma1*sigma2);epsilon=sqrt(epsilon1*epsilon2)')

    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')

    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print p1,p2,sig,eps
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system

#creating the custom force. this defines an LJ potential that defines the interaction between the wall and the rest of the system. need to create a custom potential since I wanted the potential to be equal for all atoms. also accounts for the need of openmm to have all of its nonbondedforce objects have the same exceptions and number of particles.
def create_cust_nbdforce(sigma_inp, epsilon_inp, system,nbcutoff, interaction_group):
    sigma = sigma_inp
    epsilon = epsilon_inp
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    cust_nbdforce = CustomNonbondedForce(f'4*{epsilon}*(({sigma}/r)^12-({sigma}/r)^6)')
    cust_nbdforce.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    #cust_nbdforce.setCutoffDistance(nonbonded_force.getCutoffDistance())
    cust_nbdforce.setCutoffDistance(nbcutoff*nanometer)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        cust_nbdforce.addParticle()
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        cust_nbdforce.addExclusion(p1, p2)
        if eps._value != 0.0:
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    cust_nbdforce.addInteractionGroup(interaction_group[0], interaction_group[1])
    return cust_nbdforce


def create_cust_extforce(no_particles_in_system,box_z_size,sigma_inp,epsilon_inp):
    cust_extforce_top = CustomExternalForce(f'4*epsilon_top*((sigma_top/periodicdistance(x,y,z,x,y,{box_z_size}))^12)')
    cust_extforce_top.addPerParticleParameter('sigma_top')
    cust_extforce_top.addPerParticleParameter('epsilon_top')
    cust_extforce_bot = CustomExternalForce('4*epsilon_bot*((sigma_bot/periodicdistance(x,y,z,x,y,0))^12-(sigma_bot/periodicdistance(x,y,z,x,y,0))^6)')
    cust_extforce_bot.addPerParticleParameter('sigma_bot')
    cust_extforce_bot.addPerParticleParameter('epsilon_bot')
    #cust_extforce_top.setCutoffDistance(sigma_inp_top*(2^(1/6))*nanometer)
    for index in range(no_particles_in_system):
        paramter_vector_top = (sigma_inp*nanometer,epsilon_inp*kilojoules_per_mole)
        cust_extforce_top.addParticle(index,paramter_vector_top)
        paramter_vector_bot = (sigma_inp*nanometer,epsilon_inp*kilojoules_per_mole)
        cust_extforce_bot.addParticle(index,paramter_vector_bot)
    return cust_extforce_top, cust_extforce_bot


def create_new_system_simulation(input_files, topology, positions, wall_params, temp, pressure_params, timestep, box_vect):
    forcefield = ForceField(input_files["input_ff_file"])
    if "wall_ff_file" in input_files:
        forcefield.loadFile(input_files["wall_ff_file"])
    system = forcefield.createSystem(topology, nonbondedMethod = CutoffPeriodic, nonbondedCutoff = 1.2*nanometer, constraints = HBonds)
    system = OPLS_LJ(system)
    if "pressure" in pressure_params and "direction" in pressure_params:
        if "x" in pressure_params["direction"] and "y" in pressure_params["direction"] and "z" in pressure_params["direction"] and pressure_params["pressure"] != "":
            barostat = MonteCarloBarostat(pressure_params["pressure"]*bar,temp*kelvin,25 )
            system.addForce(barostat)
            #print("adding barostat")
        elif pressure_params["pressure"] != "":
            barostat = MonteCarloAnisotropicBarostat( Vec3(pressure_params["pressure"],pressure_params["pressure"],pressure_params["pressure"])*bar, temp*kelvin, False, False, True, 25)
            system.addForce(barostat)
    if "sigma" in wall_params:
        sigma,epsilon,nbcutoff= wall_params["sigma"], wall_params["epsilon"], wall_params["nbcutoff"]#, wall_params["interaction_group"]
        num_atoms_in_system = topology.getNumAtoms()
        wall_force_top, wall_force_bot = create_cust_extforce(num_atoms_in_system,box_vect[2].z,sigma,epsilon)
        #print(wall_force.usesPeriodicBoundaryConditions())
        system.addForce(wall_force_top)
        system.addForce(wall_force_bot)
        #cust_nbdforce = create_cust_nbdforce(sigma,epsilon,system,nbcutoff, interaction_group)
        #system.addForce(cust_nbdforce)
    integrator = LangevinIntegrator(temp*kelvin, 1/picosecond,timestep*picosecond)
    simulation = Simulation(topology,system,integrator)
    simulation.context.setPositions(positions)
    simulation.context.setPeriodicBoxVectors(*box_vect)
    return system, simulation

#create an interaction group for the custom force which defines which groups of atoms interact
def create_interaction_group(solvent_interact, topology, wall_name):
    interaction_group = [[],[]]
    for res in topology.residues():
        if (res.name != solvent_name or solvent_interact) and res.name != wall_name:
            for atom in res.atoms():
                interaction_group[0].append(atom.index)
        elif res.name == wall_name:
            for atom in res.atoms():
                interaction_group[1].append(atom.index)
    return interaction_group

#detecting which type of input
def detect_input(input_name):
    phase = -1
    if "heating" in input_name:
        phase = 0
    elif "cool" in input_name:
        phase = 1
    elif "wall" in input_name:
        phase = 2
    elif "evap" in input_name:
        phase = 3
    elif "clear" in input_name:
        phase = 4
    return phase

#getting text version of detect_input
def num_to_phase(input_num, sim_type):
    temp_str =""
    if input_num == -1:
        temp_str = "none"
    elif input_num == 0:
        temp_str = "heating"
    elif input_num == 1:
        temp_str = "cool"
    elif input_num == 2:
        temp_str = "wall"
    elif input_num == 3:
        temp_str = "evap"
    elif input_num == 4:
        temp_str = "clear"
    if input_num != 0:
        temp_str = temp_str +"_" + sim_type
    return temp_str

#a generic simulation functiion, with minimization and output.
def sim(simulation,sim_type, steps, report_steps, topology, path_output, input_name, output_type):
    curr_phase = detect_input(input_name)
    print(report_steps,steps)
    if curr_phase == -1:
        print("invalid input name, returning")
        exit(0)
    else:
        output_name = name.replace(num_to_phase(curr_phase, sim_type),output_type)
        output_cp_name = output_name.replace("pdb", "xml")
    print("minimizing energy")
    simulation.minimizeEnergy()
    writefile = open(f"{path_output}checker.pdb","w")
    PDBFile.writeFile(topology, simulation.context.getState(getPositions=True).getPositions(),writefile)
    if "final" in output_type:
        temp_output_type = output_type + "_min"
        temp_output_name = name.replace(num_to_phase(curr_phase,sim_type),temp_output_type)
        writefile = open(f"{path_output}{temp_output_name}","w")
        PDBFile.writeFile(topology, simulation.context.getState(getPositions=True).getPositions(),writefile)
    simulation.reporters.append(PDBReporter(f'{path_output}output_{output_type}.pdb', report_steps))
    simulation.reporters.append(StateDataReporter(f"{path_output}data_eq.log",reportInterval = report_steps, time = True, totalEnergy = True, temperature = True))
    print(f"done minimizing, starting simulation: {output_type}")
    simulation.step(steps)
    simulation.saveState(f'{path_output}{output_cp_name}')
    writefile = open(f"{path_output}{output_name}","w")
    PDBFile.writeFile(topology,simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(),writefile)
    print(f"created outputs in {path_output}: {output_name},{output_cp_name}")

#cooling the system until 1.5 times the ideal box length, then transitioning to an anisotropic barostat to shrink the z-vector
def cool_flat(simulation, sim_type,steps, report_steps, check_steps,target_length,pressure,temp,timestep, topology, path_output, input_name, output_type):
    curr_phase = detect_input(input_name)
    if curr_phase == -1:
        print("invalid input name, returning")
        exit(0)
    else:
        output_name = name.replace(num_to_phase(curr_phase,sim_type),output_type)
        output_cp_name = output_name.replace("pdb", "xml")
    print("minimizing energy")
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter(f'{path_output}output_{output_type}.pdb', report_steps))
    simulation.reporters.append(StateDataReporter(f"{path_output}output_{output_type}.log",reportInterval = report_steps, time = True, totalEnergy = True, temperature = True))
    print(f"done minimizing, starting simulation: {output_type}")
    loops = int(steps/check_steps)
    ctr_loop = 0
    print(loops)
    for loop in range(loops):
        simulation.step(check_steps)
        print(simulation.context.getState().getPeriodicBoxVectors()[0].x)
        ctr_loop+=1
        if simulation.context.getState().getPeriodicBoxVectors()[0].x <= 1.5*target_length:
            break
    #system, simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file":wall_ff}, topology, simulation.context.getState(getPositions =True).getPositions(), {}, temp,{"pressure":pressure,"direction":"xy"}, timestep,simulation.context.getState().getPeriodicBoxVectors())
    #remain_steps = steps - ctr_loop*check_steps
    #print(remain_steps)
    #simulation.reporters.append(PDBReporter(f'{path_output}output_{output_type}.pdb', report_steps))
    #print(f"starting: anisotropic {output_type}")
    #simulation.step(remain_steps)
    simulation.saveState(f'{path_output}{output_cp_name}')
    writefile = open(f"{path_output}{output_name}","w")
    PDBFile.writeFile(topology,simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(),writefile)
    print(f"created outputs in {path_output}: {output_name},{output_cp_name}")

#add a wall to the system input
def add_wall(system,box_vect,wall_part_dia,pdb_top, positions):
    box_vect = [Vec3(box_vect[0].x,0,0), Vec3(0,box_vect[1].y,0),Vec3(0,0,2*box_vect[2].z)]
    n_part_per_line = int(np.floor(box_vect[0].x/wall_part_dia)) #number of box particles per line. 5 angstrom diameter, might change it
    pos_line = np.linspace(wall_part_dia/2,box_vect[0].x-wall_part_dia/2,n_part_per_line) #create a list of equally spaced points in a line of length box_vect[0]
    pos_list = [] #list of positions of the atoms in the system
    total_parts = n_part_per_line ** 2 # total particles in the wall
    curr_ind = system.getNumParticles() #indexes for counting
    orig_num = curr_ind
    #adding the particles of the wall to the system object
    #for part in range(total_parts):
    #    system.addParticle(0)
    #    pos_list.append([pos_line[part%n_part_per_line], pos_line[(int(np.floor(part/n_part_per_line)))%n_part_per_line], wall_part_dia/2])
    #    positions.append(Vec3(0,0,0)*nanometer)
    #    curr_ind += 1
    #adding the particles of the wall to the Topology object and determining their positions
    #pdb_top = pdb.getTopology()
    #wall_chain = pdb_top.addChain()
    min = [np.Inf,np.Inf,np.Inf]
    for op in range(orig_num):
        if positions[op].x <= min[0]: min[0] = positions[op].x
        if positions[op].y <= min[1]: min[1] = positions[op].y
        if positions[op].z <= min[2]: min[2] = positions[op].z
    for op in range(orig_num):
     #   if op >= orig_num:
     #       wall_res = pdb_top.addResidue("WR",wall_chain)
     #       pdb_top.addAtom("W00",None,wall_res )
     #       positions[op] = Vec3(pos_list[op-orig_num][0], pos_list[op-orig_num][1],pos_list[op-orig_num][2])*nanometer
      #  else:
        positions[op] = Vec3(positions[op].x ,positions[op].y ,positions[op].z - min[2] + wall_part_dia)*nanometer
    return system, box_vect, pdb_top, positions

#generate a solvent_tuple, along with other important variables for evap algorithm
def generate_solvent_tuple(pdb_top, solvent_name):
    residue_name_list = [x for x in pdb_top.residues() if x.name == solvent_name]
    wall_name_list = [x for x in pdb_top.residues() if x.name == "WR"]
    total_wall_atoms = len(wall_name_list)
    total_poly_atoms = list(residue_name_list[0].atoms())[0].index
    num_atom = pdb_top.getNumAtoms() - total_wall_atoms
    num_poly_res = residue_name_list[0].index
    num_atom_in_solv = len(list(residue_name_list[0].atoms()))
    solvent_tuple = [total_poly_atoms,num_atom]
    no_solv_mol = (solvent_tuple[1] - solvent_tuple[0] )/ num_atom_in_solv
    return solvent_tuple, no_solv_mol, num_atom_in_solv, num_poly_res

#do an evaporation
def evap(simulation,sim_type,input_ff, wall_ff,wall_params, temp,pressure_params, timestep,target_percent,no_solvent_delete, report_evap_frames, pdb_top,solvent_name,input_name, output_type):
    print("evaporating")
    solvent_tuple, no_solv_mol, num_atom_in_solv, num_poly_res = generate_solvent_tuple(pdb_top, solvent_name)
    solvent_init = (solvent_tuple[1] - solvent_tuple[0])/num_atom_in_solv
    total_parts = 0
    for residue in pdb_top.residues():
        if residue.name == "WR":
            total_parts += 1
    ctr_reporter = 1
    if "pressure" not in pressure_params:
        pressure_params["pressure"] = ""
        pressure_params["direction"] = "xyz"
    #reporter_prod = PDBReporter(f'{path_output}output_prod_evap_{ctr_reporter}.pdb',report_evap_prod_step)
    #simulation.reporters.append(reporter_prod)\
    print(simulation.context.getState().getPeriodicBoxVectors())
    curr_phase = detect_input(input_name)
    if curr_phase == -1:
        print("invalid input name, returning")
        exit(0)
    else:
        evap_name = name.replace(num_to_phase(curr_phase,sim_type),output_type)
        evap_name = evap_name.replace(".pdb","")
    print(no_solv_mol,no_solvent_delete)
    evap_loops = int(no_solv_mol/no_solvent_delete)
    print(evap_loops,report_evap_frames)
    evap_report_freq = int(evap_loops/report_evap_frames)
    if evap_report_freq == 0:
        evap_report_freq = 1
    model = Modeller(pdb_top,simulation.context.getState(getPositions=True).getPositions())
    solvent_remain = (solvent_tuple[1] - solvent_tuple[0])/num_atom_in_solv
    no_solvent_delete = int(np.floor(solvent_remain*0.25))# NEW CODE
    with open(f"{path_output}{evap_name}.log","w") as data_file:
        data_file.write("Time (ps) - Box size (cu. nm) - Box vects (nm) \n")
        for loop in range(evap_loops):
            if loop%evap_report_freq == 0:
                simulation.saveState(f'{path_output}{evap_name}_{loop}.xml')
                writefile = open(f"{path_output}{evap_name}_{loop}.pdb","w")
                PDBFile.writeFile(model.getTopology(),simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(),writefile)
                print(f"created outputs:{evap_name}_{loop}.pdb,{evap_name}_{loop}.xml" )
            simulation.step(evap_prod_steps)
            box_vect = simulation.context.getState().getPeriodicBoxVectors()
            box_size = box_vect[0].x * box_vect[1].y * box_vect[2].z
            curr_time = simulation.context.getTime()
            data_file.write(f"{simulation.context.getTime()}: {box_size:.4f} nm3 - [{box_vect[0].x:.4f},{box_vect[1].y:.4f},{box_vect[2].z:.4f}]nm\n" )
            #curr_pos_temp = [vec.z for vec in list(simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions())][solvent_tuple[0]:solvent_tuple[1]]
            #curr_pos = []
            #temp_z = 0
            #for posi in range(len(curr_pos_temp)):
            #    temp_z += curr_pos_temp[posi]
            #    if (posi+1) % num_atom_in_solv == 0:
            #        temp_z /= num_atom_in_solv
            #        curr_pos.append(temp_z)
            #        temp_z = 0
            #max_z = [-np.Inf]*no_solvent_delete
            #ind = [-1]*no_solvent_delete
            #for pos2 in range(len(curr_pos)):
            #    pos = curr_pos[pos2]
            #    for pos_maxz in range(len(max_z)):
            #        if pos > max_z[pos_maxz]:
            #            max_z.insert(pos_maxz,pos)
            #            ind.insert(pos_maxz,pos2)
            #            del max_z[len(max_z)-1]
            #            del ind[-1]
            #            break
            #print(max_z,ind)
            model_top = model.getTopology()
            residue_list = list(model_top.residues())
            solvent_tuple, no_solv_mol, num_atom_in_solv, num_poly_res = generate_solvent_tuple(model_top, solvent_name)
            ind = np.random.randint(0, high=no_solv_mol, size=no_solvent_delete)
            to_delete = []
            for delete in ind:
                to_delete.append(residue_list[num_poly_res+delete])

            model.delete(to_delete)
            #print(to_delete)

            #impleement new system
            if "sigma" in wall_params:
                sigma,epsilon,nbcutoff,interaction_group = wall_params["sigma"], wall_params["epsilon"], wall_params["nbcutoff"], create_interaction_group(solvent_interact, model.getTopology(), "WR")
                system, simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file":wall_ff}, model.getTopology(), model.getPositions(), {"sigma":sigma, "epsilon":epsilon, "nbcutoff":nbcutoff, "interaction_group":interaction_group}, temp,{"pressure":pressure_params["pressure"], "direction":pressure_params["direction"]}, timestep,box_vect)
            else:
                system, simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file":wall_ff}, model.getTopology(), model.getPositions(), {}, temp,{"pressure":pressure_params["pressure"], "direction":pressure_params["direction"]}, timestep,box_vect)
            simulation.minimizeEnergy()
            model = Modeller(model.getTopology(),simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions())
            simulation.context.setTime(curr_time)
            ctr_reporter+=1
            #print(f"{ctr_reporter=}")
            #reporter_prod = PDBReporter(f'{path_output}output_prod_evap_{ctr_reporter}.pdb',report_evap_prod_step)
            #simulation.reporters.append(reporter_prod)
            solvent_tuple[1] = model.getTopology().getNumAtoms() - total_parts
            #print(solvent_tuple)
            solvent_remain = (solvent_tuple[1] - solvent_tuple[0])/num_atom_in_solv
            no_solvent_delete = int(np.floor(solvent_remain*0.25)) # NEW
            if solvent_remain <= solvent_init*target_percent:
                break
            if solvent_tuple[1] -solvent_tuple[0] < no_solvent_delete*num_atom_in_solv: break
    clear_name = evap_name.replace("evap","clear")
    clear_name = clear_name+".pdb"
    clear_cp_name = clear_name.replace("pdb", "xml")
    simulation.saveState(f'{path_output}{clear_cp_name}')
    writefile = open(f"{path_output}{clear_name}","w")
    PDBFile.writeFile(model.getTopology(),simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(),writefile)
    print(f"created outputs: {clear_name},{clear_cp_name}, {evap_name}.log")
    return model.getTopology(), simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()

print("running version: 05 JUL 2022")
#import the parameters from config file
print("getting system params")
with open(sys.argv[1], 'r') as f:
    exec(f.read())

print("importing pdbs")
pdb_list = [""]
pdb_list[0] = input("input start PDB file (absolute path): ")
input_state_file = input("input start state file (absolute paths): ")
path_output = input("input path to output directory: ")
pdb1 = PDBFile(pdb_list[0])
pdb = PDBFile(pdb_list[0])
pdb_top = pdb.getTopology()
init_model = Modeller(pdb1.getTopology(),pdb1.getPositions())
for pdb_entry in pdb_list[1:]:
    pdb = PDBFile(pdb_entry)
    init_model.add(pdb.getTopology(),pdb.getPositions())
print("done importing, loading force field")
forcefield = ForceField(input_ff)
forcefield.loadFile(wall_ff)
print("done loading force field, setting system params")

#things to declare for setups
name = pdb_list[0]
name = name[name.rfind("/")+1:]
final_steps = wall_eq_steps

curr_phase = detect_input(name)
print(f"sim type: {sim_type}, current phase: {num_to_phase(curr_phase,sim_type)}")
#inputs
if sim_type == "near":
    direction = "xyz"
    if curr_phase == 2:
        direction = "z"
    system,simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file": wall_ff}, pdb.getTopology(), pdb.getPositions(), {}, temp, {"pressure":pressure,"direction":direction}, timestep, pdb.getTopology().getPeriodicBoxVectors())
    simulation.loadState(input_state_file)
    
    if curr_phase == 0:
        forces = {system.getForce(index).__class__.__name__: system.getForce(
            index) for index in range(system.getNumForces())}
        barostat_f = forces["MonteCarloBarostat"]
        simulation.context.getIntegrator().setTemperature(temp*kelvin)
        simulation.context.setParameter(barostat_f.Temperature(), temp*kelvin)
        #cool_flat(simulation,sim_type, cool_steps, report_cool_step, check_steps,target_length,pressure,temp,timestep, pdb.getTopology(), path_output, name, "cool_near")
        sim(simulation, sim_type,cool_steps,report_cool_step,init_model.getTopology(),path_output,name,"cool_near")
    if curr_phase <= 1:
        box_vect = simulation.context.getState().getPeriodicBoxVectors()
        positions = simulation.context.getState(getPositions = True).getPositions()
        print(len(positions))
        system, box_vect, pdb_top, positions = add_wall(system, box_vect,wall_part_dia, pdb_top, positions)
        num_atoms_in_system = pdb_top.getNumAtoms()
        
        print(len(positions))
        #interaction_group = create_interaction_group(solvent_interact, pdb_top, "WR")
        system, simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file":wall_ff}, pdb_top, positions, {"sigma":sigma, "epsilon":epsilon, "nbcutoff":nbcutoff}, temp,{"pressure":pressure, "direction":"z"}, timestep, box_vect)
       
        forces = {system.getForce(index).__class__.__name__: system.getForce(
            index) for index in range(system.getNumForces())}
        print(forces)
        sim(simulation, sim_type, wall_eq_steps,report_wall_eq_step,pdb_top, path_output,name,"wall_near")#
    if curr_phase <=5:
        forces = {system.getForce(index).__class__.__name__: system.getForce(
            index) for index in range(system.getNumForces())}
        print(forces)
        barostat_f = forces["MonteCarloAnisotropicBarostat"]
        simulation.context.getIntegrator().setTemperature(temp*kelvin)
        simulation.context.setParameter(barostat_f.Temperature(), temp*kelvin)
        positions = simulation.context.getState(getPositions=True).getPositions()
        box_vect = simulation.context.getState().getPeriodicBoxVectors()
        box_vect = [Vec3(box_vect[0].x,0,0), Vec3(0,box_vect[1].y,0),Vec3(0,0,box_vect[2].z)]
        interaction_group = create_interaction_group(solvent_interact, pdb_top, "WR")
        print(box_vect)
        system, simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file":wall_ff}, pdb_top, positions, {"sigma":sigma, "epsilon":epsilon, "nbcutoff":nbcutoff, "interaction_group":interaction_group}, temp,{"pressure":"","direction":"xyz"}, timestep, box_vect)
        new_top, new_pos = evap(simulation,sim_type,input_ff, wall_ff,{}, temp, {"pressure":pressure,"direction":"z"},timestep,0.5,no_solvent_delete, report_evap_frames, pdb_top,solvent_name,name, "evap_near")
        box_vect = [Vec3(box_vect[0].x,0,0), Vec3(0,box_vect[1].y,0),Vec3(0,0,5*box_vect[0].x)]
        interaction_group = create_interaction_group(solvent_interact, new_top, "WR")
        system, simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file":wall_ff}, new_top, new_pos, {"sigma":sigma, "epsilon":epsilon, "nbcutoff":nbcutoff, "interaction_group":interaction_group}, temp,{}, timestep, box_vect)
        new_top, new_pos = evap(simulation,sim_type,input_ff, wall_ff,{}, temp, {},timestep,0,no_solvent_delete, report_evap_frames, new_top,solvent_name,name, "evap_near2")
elif sim_type == "far":
    system,simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file": wall_ff}, pdb.getTopology(), pdb.getPositions(), {}, temp, {"pressure":pressure,"direction":"xyz"}, timestep, pdb.getTopology().getPeriodicBoxVectors())
    simulation.loadState(input_state_file)
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    barostat_f = forces["MonteCarloBarostat"]
    simulation.context.getIntegrator().setTemperature(temp*kelvin)
    simulation.context.setParameter(barostat_f.Temperature(), temp*kelvin)
    if curr_phase == 0:
        sim(simulation, sim_type,cool_steps,report_cool_step,init_model.getTopology(),path_output,name,"cool_far")
    if curr_phase <= 2:
        pdb_top, positions =evap(simulation,sim_type,input_ff, wall_ff,{}, temp, {"pressure":pressure,"direction":"xyz"},timestep,0,no_solvent_delete, report_evap_frames, pdb_top,solvent_name,name, "evap_far") #npt for pressure differnce #save the box size #can also statistics about the evap configs
    if curr_phase <= 4:
        box_vect=  simulation.context.getState().getPeriodicBoxVectors()
        system, simulation = create_new_system_simulation({"input_ff_file":input_ff, "wall_ff_file":wall_ff}, pdb_top, positions, {}, temp,{"pressure":"","direction":"xyz"}, timestep, box_vect)
        sim(simulation, sim_type,cool_steps,report_cool_step,pdb_top,path_output,name,"final_far")
