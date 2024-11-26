import espressomd
import espressomd.shapes
from espressomd import interactions
from espressomd import shapes
from espressomd.io.writer import vtf
from espressomd.interactions import HarmonicBond
import numpy as np
import charged_dev
import math
import polymerize
import sys

no_monomers =100
box_z = 100
system = espressomd.System(box_l = [80,80,box_z])

system.periodicity = [True,True,True]
system.time_step = 0.002
system.cell_system.skin = 1
#system.set_random_state_PRNG()
k_b = 0.8314
temp = 298

print("executing params")
name = "pressure-test"

#logger = logging.getLogger()
#handler = logging.FileHandler(f'{name}.log')
#logger.addHandler(handler)
#exec(open("poly-test-params.py").read())

#initialize polymerize only after setting system variables
polymerize.init(system,no_monomers,f"{name}")

#floorz = shapes.Wall(normal=[0, 0, 1], dist=0)
#z1 = system.constraints.add( penetrable=False, only_positive=False, shape=floorz)
#ceilz = shapes.Wall(normal=[0, 0, -1], dist=-system.box_l[2])
#z2 = system.constraints.add(penetrable=False, only_positive=False, shape=ceilz)
#floorx = shapes.Wall(normal=[1, 0, 0], dist=0)
#ceilx = shapes.Wall(normal=[-1, 0, 0], dist=-system.box_l[0])
#x1 = system.constraints.add( penetrable=False, only_positive=False, shape=floorx)
#x2 = system.constraints.add(penetrable=False, only_positive=False, shape=ceilx)
#floory = shapes.Wall(normal=[0, 1, 0], dist=0)
#ceily = shapes.Wall(normal=[0, -1, 0], dist=-system.box_l[1])
#y1 = system.constraints.add( penetrable=False, only_positive=False, shape=floory)
#y2 = system.constraints.add(penetrable=False, only_positive=False, shape=ceily)


print("minimizing energy")
charged_dev.minimize_e(system,0,0.03,0.1,10000,5,1e-3,1e-3)

no_part = len(system.part[:].type)

f = open(f"{name}-debug.vtf","w+t")
vtf.writevsf(system,f)
vtf.writevcf(system,f)


vtf.writevcf(system,f)

print("done with minmizing energy, setting thermostat")
system.integrator.set_vv()
system.thermostat.set_langevin(kT = k_b*temp, gamma=1, seed=42)
print("activating electrostatics")
charged_dev.set_p3m(system,k_b = k_b, temp=temp,bjerrum_length=560,q_e = 1)
print("done electrostatics")

print("doing force capping")

fcap = 20
system.force_cap = fcap

for n in range(0,500):
    system.integrator.run(50)
    fcap += fcap
    #vtf.writevcf(system,f)
    system.force_cap = fcap
    #print("progress: {0:0.3f}%".format((n*100/300)), end='\r' )

print("done with force capping")
system.force_cap = 0
print("doing uncapped simulation")
pressure = 0
for n in range(0,3000):
    system.integrator.run(100)
    pressure += system.analysis.pressure()["total"]
    if n%20 == 0:
        vtf.writevcf(system,f)
    if n%1000 == 0:
        pressure /= 1000
        print(f"{pressure=}")
        pressure = 0
    #print(system.analysis.pressure()["kinetic"])
    #print(system.analysis.energy()["kinetic"])
    #print("progress: {0:0.3f}%".format((n*100/10000)), end='\r' )

print("done with uncapped simulation")
#system.thermostat.turn_off()

#system.integrator.set_isotropic_npt(ext_pressure = 6.10e-3, piston = 20)
#system.thermostat.set_npt(kT = k_b*temp, gamma0=1, gammav = 1)
vtf.writevcf(system,f)
print("doing pressure simulation")
ctr_loop = 250
while True:

 
    pressure_2 = 0

    for n in range(0,1000):
        system.integrator.run(100)
        if n%25 == 0:
            vtf.writevcf(system,f)
        pressure_2 += system.analysis.pressure()["total"]
  
    pressure_2 /= 1000
    print(f"{pressure_2=}")
    actual_pressure = pressure_2
    ideal_pressure = 6.10193e-3
    rel_diff = (ideal_pressure-actual_pressure)/ideal_pressure
    print(actual_pressure,ideal_pressure,rel_diff)
    if abs(rel_diff) <= 1:
        print("reached desired pressure")
        break

    #system.constraints.remove(x2)
    #system.constraints.remove(y2)
    #system.constraints.remove(z2)

    gap = 0
    #if abs(rel_diff) > 10:
    #    gap = 5
    if abs(rel_diff) > 200:
        gap = 1
    elif abs(rel_diff) > 100:
        gap = 0.5
    elif abs(rel_diff) > 50:
        gap = 0.1
    elif abs(rel_diff) > 10:
        gap = 0.01
    else:
        gap = 0.001
    if rel_diff < 0:
        gap *= -1
    box_z -= gap

    #ceilx = shapes.Wall(normal=[-1, 0, 0], dist=-(box_l))
    #x2 = system.constraints.add(penetrable=False, only_positive=False, shape=ceilx)
    #ceilz = shapes.Wall(normal=[0, 0, -1], dist=-(box_z))
    #z2 = system.constraints.add(penetrable=False, only_positive=False, shape=ceilz)
    #ceily = shapes.Wall(normal=[0, -1, 0], dist=-(box_l))
    #y2 = system.constraints.add( penetrable=False, only_positive=False, shape=ceily)

    system.change_volume_and_rescale_particles(box_z,"z")

    print(system.box_l)
    ctr_loop -= 1
    if ctr_loop <= 0:
        print("too many loops")
        break
f.close()
"""
fp  = open(f'{name}-1.vtf', mode='w+t')
ctr = 2

vtf.writevsf(system,fp)
vtf.writevcf(system,fp)
if_polymerized = 0

print("doing simulation with polymerization")
for n in range(0,2000):
    system.integrator.run(50)
    if_polymerized = polymerize.check_pairs(dist_req =2.5, angle_req =math.pi/6 )
    if if_polymerized:
        fp.close()
        filename = f'{name}-'+ str(ctr)+'.vtf'
        fp = open(filename, "w+t")
        vtf.writevsf(system,fp)
        vtf.writevcf(system,fp)
        ctr+=1
        #system.integrator.set_isotropic_npt(ext_pressure = 6.10e-3, piston = 20)
        system.integrator.set_vv()
        system.thermostat.set_langevin(kT = 247.76, gamma = 10, seed =42)
        #system.thermostat.set_npt(kT = k_b*temp, gamma0=1, gammav = 1)
        if_polymerized = 0
    print("progress: {0:0.3f}%".format((n*100/2000)), end='\r' )
    if n%10 == 0:
        vtf.writevcf(system,fp)
fp.close()
"""
exit(0)
