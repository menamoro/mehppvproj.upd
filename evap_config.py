#paths
input_ff = "mehpv6_clbenz_ua.xml"
wall_ff = "evap_wall.xml"
#params

temp = 300 #in Kelvin
pressure = 1.01325 #in Bar
timestep = 0.001 #in picosecond

#name of solvent residue
solvent_name = "CLB"

#simulation type: either "near" or "far". near: cool, anisotropic barostat with wall, evaporation in nvt. far: cool, evaporation in nvt.
sim_type = "near"

#if near, these parameters are needed:
check_steps = 100000 #no. of steps per checking of periodic box length during cooling
target_length = 11 #actual ideal target, 1.5x of this length will stop the simulation and proceed to anisotropic

#simulation params
cool_steps = 20000000 # do for eq_steps timesteps
report_cool_step = 5000000 #report a pdb and log file every report_eq_step (for the equilibration phase)
wall_eq_steps = 500000 #do evap equilibration for this amount of steps
report_wall_eq_step = 25000#report a pdb every evap_eq_step (for evaporation equilibration)
#evap_loops = 2000 #number of loops to evaporate, will terminate early if solvent is all evaporated: evap_loops now automatically calculated
evap_prod_steps =  1000000 #steps per loop
no_solvent_delete = 200 #no of solvent to delete per loop
report_evap_frames = 40 # no of frames to output from the evap process
#params for the wall
sigma = 0.35 #in nm
epsilon = 10
nbcutoff = 1.2 # cutoff distance in nanometer
wall_part_dia = 0.35 # diameter of a wall particle in nanometer

#wall also attracts solvent?
solvent_interact = True
