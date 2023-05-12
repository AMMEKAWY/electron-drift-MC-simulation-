import random
import numpy as np
from matplotlib import pyplot as plt
from tioncscalc import tion
from hest import hest
from texcs import tex
from randnumgen import rand
from scipy.optimize import curve_fit


#He, Ne, Ar

e_exc=[19.8, 1.47, 11.6]					#y
e_ion=[24.6, 15.8, 15.8]				#l


E=np.linspace(10**6,1.5*10**7, 3)

#E=np.linspace(1,9*10**6, 50)


d=5*10**-3

k=1.38*10**-23
p=760*133.3*3

t=300

rad=5*10**-3
vol=(np.pi*(rad)**2)*d

n_l=2.69*10**25

N_avo=6.022*10**23
R=8.3145
n_a=n_l*p*(273/t)*vol

#n_a=p*vol/(np.sqrt(2)*k*t)

electroncharge=1.6*10**-19

me=9*10**-31
c=3*10**8
n_a=n_l*(273/t)*p*vol


def rw():

	
	yx=np.random.uniform(0,1)
	r=-np.log(1-yx)*lamda
	
	return r

def vel(en):



	A=1/((en/(me*c**2))+1)**2
	
	return c*np.sqrt(1-A) 

	
	
def n_inelastic_per_meter(x,y):

	if (x==0):
	
		return 0
		
	else:

		return (1-np.exp(-x/y))
		#return np.exp(-y/x)
	
		
def n_p_ions(x,y,l):

	if (x==0):
	
		return 1
		
	else:

		#return ((1-np.exp(-(x)/y))*(1-np.exp(-x/l)))
		#return ((1-np.exp(-x/y))*(np.exp(-l/x)))
		return ((np.exp(-l/x)))
	
	
def mainfunc(ener, ionE, excE , ln, elas, inela, inelaion, ioni, distance):

	if (ener[-1]>=0):
		
		#ener=[energy]
		
		
		if(ener[-1]<ionE[ln]):
				
					
			if (ener[-1]<excE[ln]):
						
					
				#elastic
				
				#elas.append(elas[-1]+1)
				#elas.pop(0)
					
				ener=np.append(ener,ener[-1])
				ener=np.delete(ener, 0)	
					
				#ioni.append(ioni[-1]+e)
				
				#distance.append(distance[-1]+x)	
				#distance.pop(0)
					
				#energy=energy+e

				return inelaion, ener[-1]		
					
					
			if (ener[-1]>=excE[ln] and energy[-1]<=ionE[ln] ):
			
							
				if (n_inelastic_per_meter(ener[-1],excE[ln])<0.5):
					
					#elastic

					#elas.append(elas[-1]+1)
					#elas.pop(0)
					
					ener=np.append(ener, ener[-1])
					ener=np.delete(ener, 0)
					
					#ioni.append(ioni[-1]+e)
						
					#energy=energy+e
					
					#distance.append(distance[-1]+x)	
					
					#distance.pop(0)
						
					return inelaion, ener[-1]
					
							
					
				
				else:
					
					#excitation
					
					#inela.append(inela[-1]+1)

					#inela.pop(0)
					
					ener=np.append(ener, ener[-1]-e_exc[2])
					
					ener=np.delete(ener, 0)
							
					#ioni.append(ioni[-1]-e_exc[2]+e)		
					
					#distance.append(distance[-1]+x)
					
					#distance.pop(0)
									
					return inelaion, ener[-1]

					#print("here4")
		else:
									
			if (n_p_ions(ener[-1],excE[ln],ionE[ln])<0.5):

									
				if (n_inelastic_per_meter(ener[-1],excE[ln])<0.5):
					
					
					#elastic
					
					#elas.append(elas[-1]+1)
					#elas.pop(0)
				
					ener=np.append(ener, ener[-1])
					ener=np.delete(ener, 0)
		
					return inelaion, ener[-1]
							
					
				
				else:
					
					
					#excitation

					#inela=1
					
					ener=np.append(ener, ener[-1]-e_exc[2])
					ener=np.delete(ener, 0)

					#inela.append(inela[-1]+1)
					#inela.pop(0)
					
					#ioni.append(ioni[-1]-e_exc[2]+e)
					
					#distance.append(distance[-1]+x)	
					#distance.pop(0)
					
					return inelaion, ener[-1]
							
								
				
			else:
					
				#ionization
				
				inelaion.append(inelaion[-1]+1)
				inelaion.pop(0)
				
				ener=np.append(ener, (ener[-1]-e_ion[2]))
				ener=np.delete(ener, 0)
				
				#ioni.append(ioni[-1]-e_ion[2]+e)
				
				#distance.append(distance[-1]+x)	
				#distance.pop(0)
				
				return inelaion, ener[-1]
	
	else:
	
		return(inelaion, ener[-1] )



i=0

x2=x=np.array([0])
y2=y=np.array([0])
vx=np.array([0])
vy=np.array([0])
xd=np.array([0])
yd=np.array([0])

x=np.array([0])
y=np.array([0])
j=1

en=np.array([0])

h=10**-15*6

h=10**-8

#h=10

drift=np.array([0])

r1=0

r2=0

inelaio=[0]

#th=np.linspace(0, 2*np.pi, 10**5)

energy=[0]

l=[]

#lamda=1/(n_a*((np.pi*(180*10**-12)**2)))

#print(1/lamda)

xl=0
kr=1000



#E=[5*10**6]

#E=[0]

l=0

#hin=np.logspace(12, 2, 10, base=10)

#print(hin)

vd=[]

while l != len(E):

	energy=np.array([0.00000001])
	#vd=[np.sqrt(2*(energy[-1]*electroncharge)/me)]
	

	inelaio=[0]
	i=0
	xl=0

	#p=[0]
	
	v=[0]
	
	
	
	h=5*10**-14
	
	theta=[np.random.uniform(-np.pi/2, np.pi/2)]
	phi=np.random.uniform(0, 2*np.pi)	
		
	x2=x=np.array([0])
	y2=y=np.array([0])
	z2=z=np.array([0])
	vx=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.cos(theta)*np.cos(phi)])
	vy=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(theta)*np.cos(phi)])
	vz=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(phi)])
	xd=np.array([0])
	yd=np.array([0])
	zd=np.array([0])
	
	#x=np.array([0])
	#y=np.array([0])

	m=0
	
	#print(h)

	while x[-1] < d/1000:
	
	
		#lamda=1/(n_a*((hest(energy[-1])+tion(2, energy[-1])+tex(2, energy[-1]))))
		
		lamda=1/(n_a*((hest(0))))
		
		xd=np.append(xd,(rw()*(np.cos(theta[m-1])*np.cos(phi))-x[-1]))
		#xd=np.append(xd,(-rw()*(np.cos(theta[m-1]))+x[-1]))
		yd=np.append(yd,rw()*np.sin(theta[m])*(np.cos(phi))-y[-1])
		zd=np.append(zd,rw()*np.sin(phi)-z[-1])
		
		xd=np.append(xd,(rw()*(np.cos(theta[m-1]))*np.cos(phi)+x[-1]))
		#xd=np.append(xd,(-rw()*(np.cos(theta[m-1]))+x[-1]))
		yd=np.append(yd,rw()*np.sin(theta[m])*np.cos(phi)+y[-1])
		zd=np.append(zd,rw()*np.sin(theta[m])*np.sin(phi)+z[-1])
	
		delt=(electroncharge/me)*E[l]
		
		while ((x[-1]) < (xd[-1])) : 
	
			#x=[0]
		
			#y=[0]
				
			vx=np.append(vx,(vx[-1]+delt*h))
		
			vy=np.append(vy, vy[-1])
		
			vz=np.append(vz, vz[-1])
		
			x=np.append(x,vx[-1]*h+delt*h**2/2+x[-1])
		
			y=np.append(y,vy[-1]*h+y[-1])
			
			z=np.append(z,vz[-1]*h+z[-1])
			
			
			if ((y[-1]) < (yd[-1])): 
		
				break
		
			elif ((z[-1]) < (zd[-1])): 
		
				break
		
			elif ((x[-1]) > (xd[-2])):
			
				break
			
			elif ((y[-1]) > (yd[-2])):
			
				break
			
			elif ((z[-1]) > (zd[-2])):
			
				break
		
			#r1=np.sqrt((x[-1]-x2[-1])**2+(y[-1]-y2[-1])**2)
			
			#print(round(r1,9),"	", round(dr,9),"	 ", y[-1]+yd[-1])
						
			#j+=1
		
	
		energy=np.append(energy, 0.5*me*(vx[-1]**2+vy[-1]**2+vz[-1]**2)/electroncharge)
		
		energy=np.append(energy, me*c**2*(1/np.sqrt(1-(vx[-1]**2+vy[-1]**2)/c**2)-1)/electroncharge)

		#print(energy[-1],"		", x[-1]/xd[-1])
		
		#inelaio, energ =mainfunc(energy, e_ion, e_exc , 2, 0, 0, inelaio, 0, 0)

		'''
		if (energy[-1]> e_exc[2] and energy[-1] < e_ion[2]):
		
			energy=np.append(energy, energy[-1]-e_exc[2])
		
		elif(energy[-1] > e_ion[2]):
		
			energy=np.append(energy, energy[-1]-e_ion[2])
			inelaio=np.append(inelaio, inelaio[-1]+1)
			inelaio.pop(0)
		else:
		
			energy=np.append(energy, energy[-1])

		'''
		#energy=np.append(energy, energ)
		
		

		theta.append(np.random.uniform(0, 2*np.pi))
		phi=np.random.uniform(0, 2*np.pi)
		
		
		vx=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.cos(theta)*np.cos(phi)])
		
		vy=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(theta)*np.cos(phi)])
			
		vz=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(phi)])		

	
		#h=(me/(electroncharge*E[l]))*3*10**4				

		#print(xd[-1],"		", x[-1])
	
		#theta=random.random()*3.14/2
		
		
	
		#print(xd[-1])
	
		'''
		#delt=0
		dr=np.sqrt(xd[-1]**2+yd[-1]**2)
	
		if x2[-1]<0:
		
			x2[-1]=0
	
		r1=0
		
		j=0	
		
			
		#print(xl)
		
		
	
		x2=np.append(x2,xd[-1]+x[-2])
		y2=np.append(y2,yd[-1]+y[-2])

		#energy.append(E[l]*(r1*np.cos(theta[m])))
		
	
		#print(v)
	
		#theta.append(np.random.uniform(-np.pi, np.pi))
	
		#energy.append((1/2)*me*(vx[-1]**2+vy[-1]**2)/electroncharge)
		
		#energy.append(me*c**2*(1/np.sqrt(1-(v[-1]/c)**2))/electroncharge)
			
		#energy.append(E[l]*(xl))
	
		#energy=(1/2)*me*(vx[-1]**2+vy[-1]**2)/electroncharge



		#print(energy[-1])
		
		'''
		
		#print(energy[-1])
		
		#en=np.append(en, energy/electroncharge)
			
		
		#r2=r2+x2[-1]
		
		#x=np.append(x, xd[-1])
		#y=np.append(y, yd[-1])
			
		m+=1	
		i+=1
	
	#print(round(E[l]/100), (np.sqrt(np.mean(energy)*electroncharge*2/(me))*10**-7))
	
	vd.append(vel(np.mean(energy)*electroncharge)*10**-7)
	
	print(round(E[l]/100), vel(np.mean(energy)*electroncharge)*10**-7,"	", inelaio[-1])
	
	l+=1

#print(len(vd))

#vd[0]=0
	
#plt.plot(x,z)
#plt.plot(x, y)

plt.scatter(E/100, vd)

#plt.plot(E/100, vd)

plt.show()

