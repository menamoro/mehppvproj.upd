
import numpy as np
import espressomd
from espressomd import electrostatics
from espressomd.io.writer import vtf
"""
np.random.seed(42)
print([-1]*50)
box_le = 150
bjerrum_length = 1 # angstrom
# C = lb * T * k_b / q^2 , lb in angstrom
syst = System(box_l =  box_le  * np.ones(3))

#tutorial says assume unit system, we gon use a real system
temp = 1 # K
k_b = 1 # in E/K
q_e = 1 # as in, 1 eV

syst.time_step =0.001 # in ps
syst.cell_system.skin = 0.4 # in angstrom
n_lines = 4
exec(open("read-lj.py").read())

rod_type = 0
ctrion_type =1
"""
def setup_rod_and_counterions(system,ion_valency, counterion_type,rod_charge_dens, N_rod_beads, rod_type): #rod_charge_dens is linear charge density
    rod_l = syst.box_l[2]
    total_rod_charge = rod_charge_dens * rod_l
    charge_per_rod_mol = total_rod_charge / N_rod_beads

    rod_beads_pos = np.column_stack(([syst.box_l[0]/2] * N_rod_beads, [syst.box_l[1]/2]* N_rod_beads, np.linspace(0,rod_l,N_rod_beads)))

    syst.part.add(type= [int(rod_type)]*int(N_rod_beads), pos = rod_beads_pos, q = np.ones(N_rod_beads) * charge_per_rod_mol, fix= [[True]*3]*N_rod_beads )

    no_ctrion = abs(int(total_rod_charge / ion_valency))

    syst.part.add(type=[int(counterion_type)]*int(no_ctrion), pos = np.random.random((no_ctrion, 3)) *[syst.box_l[0],syst.box_l[1],syst.box_l[2]], q = np.ones(no_ctrion) * ion_valency)

    return 0



#add salts


#syst.part.add(type=[2]*25, pos=np.random.random((25, 3)) *[syst.box_l[0],syst.box_l[1],syst.box_l[2]], q = [-1]*25)
#syst.part.add(type=[3]*25, pos=np.random.random((25, 3)) *[syst.box_l[0],syst.box_l[1],syst.box_l[2]], q = [1]*25)

def set_p3m(system, k_b, temp, bjerrum_length, q_e):
    p3m_params = {"prefactor": k_b * temp *bjerrum_length * (q_e ** 2), "accuracy": 1e-4 ,"check_neutrality":False}

    p3m  = electrostatics.P3M(**p3m_params)
    system.actors.add(p3m)
# steepest descent
def minimize_e(system,max_f, max_d, gam, max_std_step, time_per_interval, f_tol,e_tol):
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
        if ctr_energy > 5 and ctr > 5:
            break
        f_prev = curr_f
        e_prev = e_curr

def minimize_e_write(system,max_f, max_d, gam, max_std_step, time_per_interval, f_tol,e_tol,file):
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
        vtf.writevcf(system,file)
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
        if ctr_energy > 5 and ctr > 5:
            break
        f_prev = curr_f
        e_prev = e_curr

#production
"""
def setup_profile_calculation(system,delta_N, ion_types, r_min, n_radial_bins):
    ia = []
    for i in ion_types:
        ia.append(system.part.select(type=i).id)
    rdf_bc_list = np.empty(shape=(len(ia),len(ia)))
    rdf_acc_list = np.empty(shape=(len(ia),len(ia)))
    for x in range(len(ia)):
        for y in range(x,len(ia)):
            r_temp = observables.RDF(ids1=ia[x], ids2=ia[y], min_r  = r_min, n_r_bins = n_radial_bins)
            rdf_bc_list[x,y] = r_temp.bin_centers()
            rdf_acc_list[x,y] = accumulators.MeanVarianceCalculator(obs= r_temp, delta_N = delta_N)
            system.auto_update_accumulators.add(rdf_acc_list[x,y])
    return 0
"""

def setup_profile_calculation(system,delta_N,ion_types,r_min,n_radial_bins):
    ia = []
    ctp = math.CylindricalTransformationParameters(center= np.array(system.box_l)/2, axis =[0,0,1], orientation=[1,0,0])
    cdp_acc_list = {}

    for x_int in ion_types:
        cdp_temp = observables.CylindricalDensityProfile(ids=system.part.select(type=x_int).id, min_r = r_min, n_r_bins = n_radial_bins, max_r = syst.box_l[0]/2, max_z = syst.box_l[1]/2, min_z = 0, transform_params=ctp)
        cdp_acc_list[x_int]= accumulators.MeanVarianceCalculator(obs=cdp_temp, delta_N = delta_N)
        cdp_bc = cdp_temp.bin_centers()
        system.auto_update_accumulators.add(cdp_acc_list[x_int])

    return cdp_bc,cdp_acc_list

def capped_sim(system, cap, loops, time_per_loop,file):
    cappy = cap
    system.force_cap = cappy
    for i in range(loops):
        system.integrator.run(time_per_loop)
        vtf.writevcf(system,file)
        cappy += cap
        system.force_cap = cappy
        print(system.analysis.energy()["total"])
        print("progress: {0:0.2f}%".format ( i*100/loops ), end='\r')

def production_run(system, steps_per_t, total_t):
    system.time = 0
    last_p = 0
    total_iter = int(total_t/steps_per_t)
    for i in range(total_iter):
        system.integrator.run(steps_per_t)
        curr_p =  int( (i)*100/total_iter)
        if curr_p > last_p:
            print("progress: {0}%".format ( curr_p ), end='\r')
            last_p += 1
    print("progress: 100%")
    print("done!")

def clear_system(system):
    system.time = 0.0
    system.part.clear()
    system.thermostat.turn_off()
    system.auto_update_accumulators.clear()

def add_salt(system,anion_type,anion_charge,cation_type,cation_charge,n_pairs):
    system.part.add(type=[anion_type]*n_pairs, pos=np.random.random((n_pairs,3)) * [system.box_l[0],system.box_l[1],system.box_l[2]], q = [anion_charge]*n_pairs)
    system.part.add(type=[cation_type]*n_pairs, pos=np.random.random((n_pairs,3)) * [system.box_l[0],system.box_l[1],system.box_l[2]], q = [cation_charge]*n_pairs)
    return 0

"""
clear_system(syst)


runs = {1:{"ion_valency":-2, "counterion_type":1, "rod_charge_dens":1, "N_rod_beads":150,"rod_type":0}, 2:{"ion_valency":-1,"counterion_type":1, "rod_charge_dens":2,"N_rod_beads":150,"rod_type":0}}
fig, ax = plt.subplots(figsize=(10,6))
charge_func = np.zeros(700)
for i in runs:
    clear_system(syst)
    setup_rod_and_counterions(syst,**runs[i])
    p3m.tune()
    print("-------------minimizing--------------- system no:" +str(i)  )
    minimize_e(syst,0,0.01*1,30,10000,10,1e-3)
    syst.thermostat.set_langevin(kT = k_b * temp, gamma = 1, seed = 39)
    syst.integrator.set_vv()
    print("----------------warmup---------------- system no:" +str(i)  )
    production_run(syst,1000,10000)
    syst.time =0
    steps_per_t_v = 100
    cdp_bc_l, cdp_acc_l = setup_profile_calculation(syst,delta_N  = steps_per_t_v, ion_types=[0,1], r_min = 0, n_radial_bins= 700)
    print("--------------production-------------- system no:" +str(i) )
    production_run(syst, steps_per_t = steps_per_t_v, total_t = 10000)
    list_of_ions = [runs[i]["counterion_type"]]
    for m in list_of_ions:
       cdp_acc_ctrion = cdp_acc_l[m].mean()
       h_list = np.array(cdp_acc_ctrion[:,0,0]) * cdp_bc_l[:,0,0,0]
       cum_hist = np.cumsum(h_list)
       cum_hist /= cum_hist[-1]
       charge_func = charge_func +  (cum_hist*syst.part.select(type=m).q[0])
       ax.plot(cdp_bc_l[:,0,0,0], cum_hist, label="the " + str(m) + " type one, of " + str(i) + " run")


for i in runs:
    clear_system(syst)
    setup_rod_and_counterions(syst,**runs[i])
    add_salt(syst,2,-1,3,1,10)
    p3m.tune()
    print("-------------minimizing--------------- system no:" +str(i) + " with salt" )
    minimize_e(syst,0,0.01*1,30,10000,10,1e-3)
    syst.thermostat.set_langevin(kT = k_b * temp, gamma = 1, seed = 39)
    syst.integrator.set_vv()
    print("----------------warmup---------------- system no:" +str(i) + " with salt" )
    production_run(syst,1000,3000)
    syst.time =0
    steps_per_t_v = 100
    cdp_bc_l, cdp_acc_l = setup_profile_calculation(syst,delta_N  = steps_per_t_v, ion_types=[0,1,2,3], r_min = 1, n_radial_bins= 700)
    cdp_bc_l = np.array(cdp_bc_l[:,0,0,0])
    #cdp_bc_l /= cdp_bc_l[-1]

    print("--------------production-------------- system no:" +str(i) + " with salt")
    production_run(syst, steps_per_t = steps_per_t_v, total_t = 20000)
    list_of_ions = [runs[i]["counterion_type"],2,3]
    for m in list_of_ions:
        cdp_acc_ctrion = cdp_acc_l[m].mean()

        h_list = cdp_acc_ctrion[:,0,0]* cdp_bc_l[1:]
        cum_hist = np.cumsum(h_list)
        cum_hist /= cum_hist[-1]
        print(cum_hist)
        charge_func = charge_func +  (cum_hist*syst.part.select(type=m).q[0])

    #    ax.plot(cdp_bc_l[:,0,0,0], cum_hist, label="the " + str(m) + " type one, of " + str(i) + " run, with salt")

#ax.plot(cdp_bc_l[:,0,0,0], charge_func, label="fuck tyou")
ax.legend()
#ax.set_xscale('log')

plt.xlabel('r')
plt.ylabel('wow')
plt.show(block=True)

v = visualization.openGLLive(syst)
v.run()
"""
