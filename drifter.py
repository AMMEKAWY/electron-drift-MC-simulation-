import random
import numpy as np
from matplotlib import pyplot as plt
from tioncscalc import tion
from hest import hest
from texcs import tex
#from randnumgen import rand
#from scipy.optimize import curve_fit
from iterr import *
import random

#He, Ne, Ar

e_exc=[19.8, 1.47, 11.6]					#y
e_ion=[24.6, 15.8, 15.8]				#l


E=np.linspace(10**6,1.*10**7, 3)

#E=np.linspace(1,9*10**6, 50)


d=5*10**-3/50
n=1

k=1.38*10**-23
p=760*133.3*n



t=300

rad=5*10**-3
vol=(np.pi*(rad)**2)*d

n_l=2.69*10**25

N_avo=6.022*10**23
R=8.3145
#n_a=n_l*p*(273/t)*vol

#n_a=p*vol/(np.sqrt(2)*k*t)

electroncharge=1.6*10**-19



me=9*10**-31
c=3*10**8
n_a=n_l*(273/t)*p
n_a=p/(k*t)

m_at=6.6335209*10**-26

e_at=(3/2)*k*t
v_at=np.sqrt(2*e_at/m_at)

def rw():

	
	#lamda=10**-50
	yx=np.random.uniform(0,1)
	ro=-np.log(1-yx)*(lamda)
	
	#print(ro)
	
	return ro

def vel(en):



	A=1/((en/(me*c**2))+1)**2
	
	return c*np.sqrt(1-A) 

	
	
def n_inelastic_per_meter(x,y):

	if (x==0):
	
		return 0
		
	else:

		#return (1-np.exp(-x/y))
		return np.exp(-y/x)
		#return 1
		
def n_p_ions(x,y,l):

	if (x==0):
	
		return 1
		
	else:

		#return ((1-np.exp(-(x)/y))*(1-np.exp(-x/l)))
		return ((1-np.exp(-x/y))*(np.exp(-l/x)))
		#return ((np.exp(-l/x)))
		
	
	
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
			
							
				if (n_inelastic_per_meter(ener[-1]-excE[ln],excE[ln])<0.5):
					
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

									
				if (n_inelastic_per_meter(ener[-1]-excE[ln],excE[ln])<0.5):
					
					
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

'''

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

#h=10**-15*6

#h=10**-8

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
'''


#E=[5*10**6]

#E=[0]

l=2

#hin=np.logspace(12, 2, 10, base=10)

#print(hin)

vd=[]

while l != len(E):


	
	counter=0
	
	ion=[0]
	time=[]
	mfp=[]
	xd=np.array([0])
	x_del=[0]
	alp=[]
	
	while counter != 30:

		energy=np.array([0.00000001])
		#vd=[np.sqrt(2*(energy[-1]*electroncharge)/me)]
	

		inelaio=[0]
		i=1
		xl=0
	
		#p=[0]
	
		v=np.array([0])
	
	
	
		h=10**-15
	
		theta=np.random.uniform(0, 2*np.pi)
		phi=np.random.uniform(0, np.pi/2)	
		
		x2=x=np.array([0])
		y2=y=np.array([0])
		z2=z=np.array([0])
		vx=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.cos(theta)*np.cos(phi)])
		vy=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(theta)*np.cos(phi)])
		vz=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(phi)])
		
		yd=np.array([0])
		zd=np.array([0])
	
		#vd=[np.sqrt(vx[-1]**2+vy[-1]**2+vz[-1]**2)]
	
		v=np.array([np.sqrt(2*(energy[-1]*electroncharge)/me)*np.cos(theta)*np.cos(phi)])
	
		#x=np.array([0])
		#y=np.array([0])

		m=0
	
		ttau=0
	
		#print(h)

		i=1
		ttau=0

		while x[-1] < d:
	
	
			lamda=1/(n_a*((hest(energy[-1])+tion(2, energy[-1])+tex(2, energy[-1]))))
			#lamda=1/(n_a*((hest(energy[-1])+tion(2, energy[-1]))))
			#lamda=1/(n_a*((hest(energy[-1]))))
			#lamda=1/(n_a*((np.pi*(180*10**-12)**2)))
			#ind=np.random.uniform(-1, 1)
		
			ind=0
		
			#print(x[-1])
			
			#print(lamda)
		
			xd=np.append(xd,(rw()))
			#xd=np.append(xd,(-rw()*(np.cos(theta[m-1]))+x[-1]))
			yd=np.append(yd,rw())
			zd=np.append(zd,rw())
		
			xd=np.append(xd,(-rw()))
			#xd=np.append(xd,(-rw()*(np.cos(theta[m-1]))+x[-1]))
			yd=np.append(yd,-rw())
			zd=np.append(zd,-rw())
	
		
			#print(lamda)
	
			delt=(electroncharge/me)*E[l]
		

		
			xf,yf,zf=(x[-1]),(y[-1]), (z[-1])
			
			xd1,yd1,zd1=(xd[-2]), (yd[-2]) ,(zd[-2])
			
			xd2,yd2,zd2= (xd[-1]) ,(yd[-1]) ,(zd[-1]) 
			vxf,vyf,vzf=(vx[-1]+ind*v_at),(vy[-1]+ind*v_at), (vz[-1]+ind*v_at)
			
			tt=(iterr(xf, yf, zf, xd1, yd1, zd1, xd2, yd2, zd2, vxf, vyf, vzf, delt, h))
		
			'''
			if tt[0]/d>1.05:
			
				counter=counter-1
				break		
			'''
			vx=np.append(vx,tt[3])
		
			vy=np.append(vy, tt[4])
		
			vz=np.append(vz, tt[5])
		
			x=np.append(x,tt[0])
		
			y=np.append(y,tt[1])
				
			z=np.append(z,tt[2])
	
			energy=np.append(energy, 0.5*me*(vx[-1]**2+vy[-1]**2+vz[-1]**2)/electroncharge)


		
			if energy[-1]>3000:
		
				print("out")

		
			inelaio, energ =mainfunc(energy, e_ion, e_exc , 2, 0, 0, inelaio, 0, 0)

			energy=np.append(energy, energ)
		
		

			theta=random.uniform(0, 2*np.pi)
			phi=np.random.uniform(0, np.pi)
		
		
			vx=np.append(vx, np.sqrt(2*(energy[-1]*electroncharge)/me)*np.cos(theta)*np.cos(phi))
		
			vy=np.append(vy, np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(theta)*np.cos(phi))
			
			vz=np.append(vz, np.sqrt(2*(energy[-1]*electroncharge)/me)*np.sin(phi))		

			
			v=np.append(v,np.sqrt(vx[-1]**2+vy[-1]**2+vz[-1]**2))
		
			ttau=ttau+(me*abs(v[-1]-v[-2])/(electroncharge*E[l]))
			
			h=me*abs(v[-1]-v[-2])/(electroncharge*E[0])*0.0001
			
			#alp.append(n_a*(tex(2, energy[-1])+tion(2, energy[-1])))
			
			if inelaio[-1]>0:
			
				ion.append(inelaio[-1])
		
			xl=ttau/i
			
			#print(round(x[-1]*100/d),"		",inelaio[-1])
		
			#print((inelaio[-1]/i)*1/(100*lamda))
				
			m+=1	
			i+=1
		
		counter+=1	

		
	
		time.append(xl)

		veld=(electroncharge*E[-1]*np.mean(time)*10**-7/me)
	
		vd.append(veld)
	
		xll=(vd[-1])*10**7*np.mean(time)/i
		
		#print(i)
		
		A=(1/xll)
	
		B=A*e_ion[2]/E[l]
	
		#print(vd[-1])
	
		alphat=((np.mean(ion))/(d))*(10**-2)
		
		alp.append(alphat)
	
		#alphat=((np.mean(ion)/i)/(xll))*(10**-2)
	
		#alphat=A*np.exp(-B)
	
		print(round(E[l]/100),"		",np.mean(alp),"	", counter,"	",np.mean(vd))
	
	l+=1

#print(np.mean(mfp)/np.mean(time)*10**-7)

#print(len(vd))

#vd[0]=0
	
#plt.plot(x,z)
#plt.plot(x, y)

#plt.scatter(E/100, vd)

#plt.plot(E/100, vd)

#plt.show()

