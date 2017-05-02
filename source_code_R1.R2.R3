#VPython source code for simulating the motion of 3 protons through the B-field with different resistances in
# the circuits, R1,R2, and R3

GlowScript 2.4 VPython
from visual import *

#variable definitions
q=1.6*10**-19   #charge
m=1.7*10**-27  #mass
mu0=(4*pi)*10**-7 
Vb=100     #voltage of battery source
#R=20000       #reistance of resistor
R1=1000
R2=2000
R3=4000
print("R1 (blue) is " +R1+ " ohms.")
print("R2 (purple) is " +R2+ " ohms.")
print("R3 (yellow) is " +R3+ " ohms.")
C=.003       #Capacitance
It=0     #current at time t
r0x=.1   #initial distance between proton and wire
r0y=-20
v0x=0
v0y=10

myboxL=box(pos=vector(-.2,0,0), length=.10, height=40, width=.01)
myboxR=box(pos=vector(5.2,0,0), length=.10, height=40, width=.01)
proton1=sphere(pos=vector(r0x,r0y,0),radius=0.03,color=vector(0,1,1),charge=q,mass=m,v=vector(v0x,v0y,0),trail=curve(color=vector(0,0,1)))
proton2=sphere(pos=vector(r0x,r0y,0),radius=0.03,color=vector(0,1,1),charge=q,mass=m,v=vector(v0x,v0y,0),trail=curve(color=vector(1,0,1)))
proton3=sphere(pos=vector(r0x,r0y,0),radius=0.03,color=vector(0,1,1),charge=q,mass=m,v=vector(v0x,v0y,0),trail=curve(color=vector(1,1,0)))

#Particle dynamics
t=0
dt=.001
#It=(Vb/R)*(exp(-t/(R*C)))
#print (It)
k=0
while t<=2:
    rate(800)
    t=t+dt
    It=(Vb/R1)*(exp(-t/(R1*C)))       #current at time t
    BL=vector(0,0,(mu0*It)/(4*pi*abs(proton1.pos.x)))         #magnitude of the magnetic field produced by the wire at the proton's position
    BR=vector(0,0,-(mu0*It)/(4*pi*abs(5-proton1.pos.x))) 
    B1=BL+BR
    Fm=q*cross(proton1.v,B1)     #magnetic force on the proton
    deltav=(1/proton1.mass)*Fm*dt        #change in velocity of the proton, using change in momentum
    vxavg=proton1.v+(deltav/2)                   #average velocity of the proton over the time frame of dt
    proton1.pos.x=proton1.pos.x+(vxavg.x*dt)
    proton1.pos.y=proton1.pos.y+(vxavg.y*dt)
    proton1.v=proton1.v+deltav
    proton1.trail.append(proton1.pos)
t=0
It=0    #zero out the current
while t<=2:
    rate(800)
    t=t+dt
    It=(Vb/R2)*(exp(-t/(R2*C)))       #current at time t
    BL=vector(0,0,(mu0*It)/(4*pi*abs(proton2.pos.x)))         #magnitude of the magnetic field produced by the wire at the proton's position
    BR=vector(0,0,-(mu0*It)/(4*pi*abs(5-proton2.pos.x)))
    B2=BL+BR
    Fm=q*cross(proton2.v,B2)     #magnetic force on the proton
    deltav=(1/proton2.mass)*Fm*dt        #change in velocity of the proton, using change in momentum
    vxavg=proton2.v+(deltav/2)                   #average velocity of the proton over the time frame of dt
    proton2.pos.x=proton2.pos.x+(vxavg.x*dt)
    proton2.pos.y=proton2.pos.y+(vxavg.y*dt)
    proton2.v=proton2.v+deltav
    proton2.trail.append(proton2.pos)
t=0
It=0    #zero out the current
while t<=2:
    rate(800)
    t=t+dt
    It=(Vb/R3)*(exp(-t/(R3*C)))       #current at time t 
    BL=vector(0,0,(mu0*It)/(4*pi*abs(proton3.pos.x)))         #magnitude of the magnetic field produced by the wire at the proton's position
    BR=vector(0,0,-(mu0*It)/(4*pi*abs(5-proton3.pos.x)))
    B3=BR+BL
    Fm=q*cross(proton3.v,B3)     #magnetic force on the proton
    deltav=(1/proton3.mass)*Fm*dt        #change in velocity of the proton, using change in momentum
    vxavg=proton3.v+(deltav/2)                   #average velocity of the proton over the time frame of dt
    proton3.pos.x=proton3.pos.x+(vxavg.x*dt)
    proton3.pos.y=proton3.pos.y+(vxavg.y*dt)
    proton3.v=proton3.v+deltav
    proton3.trail.append(proton3.pos)

