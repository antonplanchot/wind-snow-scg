import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d  # Fonction pour la 3D
import numpy as np
import netCDF4 as nc

data=nc.Dataset(r"C:\Users\Toshiba\Documents\stage L3\Documents\code_amory_5_juin\forcings_20100101-20200101.nc")

#extraire les variables meteo du fichier netcdf
U=data.variables['U2'][:,0,0]   #[m/s]

##variables temporelles
date=[]
an=2010
while an<2020:
    a=0
    if an==2012 or an==2016:
        nb_jour=366
    else:
        nb_jour=365
    while a<nb_jour*24:
        date.append(an+a/(24*nb_jour))
        a+=1
    an+=1
date.remove(2010)
date.append(2020)

Lmois=[31,28,31,30,31,30,31,31,30,31,30,31]
mois=[]
moi=1
while moi<13:
    a=0
    nb_jour=Lmois[moi-1]
    while a<nb_jour*24:
        mois.append(moi+a/(24*nb_jour))
        a+=1
    moi+=1

##declaration des constantes
rhos0=300
rhoice=920

vk=0.4   #constante de Von Karman   sans dimension

d=0.5    #dendricite
s=0.5   #sphericite

C_D=10**-3   #coefficient de drag pour le moment   sans dimension
z0=10**-4   #longueur de rugosite aerodynamique   m

#calcule u* a partir de U
ustar=U*vk/np.log(2/z0)

#indice d'erodabilite
iER=0.75*d-0.5*s+0.5

#vitesse de friction seuil standard
ustarT0=(C_D**(0.5))*(np.log(2.868)-np.log(1+iER))/0.085

##plot 2D

size=15
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)

x = np.arange(0,100, 1)
y = np.arange(251,451,2)
z=np.zeros((len(x),len(y)))
for i in range(len(x)):
    for j in range(len(y)):
        z[i,j]=0.25*x[i]*vk/(np.log(2/z0))-ustarT0*np.exp(rhoice/rhos0-rhoice/y[j])

alpha = [ '300', '350','400','450']
alpha2 = ['0','5','10', '15','20']


fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(z,origin='lower')


plt.plot(x,[4*ustarT0*np.exp(rhoice/rhos0-rhoice/k)/(vk/np.log(2/z0)) for k in y], color='r', label='u*=u*t',linewidth=5)
fig.colorbar(cax,label='u*-u*t [m/s]')
#plt.plot([24,24],[0,100],color='aqua',label='ρs=ρ0',linewidth=5)
plt.plot([99,99],[0,100],color='aqua',label='ρs=ρmax',linewidth=5)
ax.set_xticks(range(24,100,25))
ax.set_yticks(range(0,100,20))
ax.set_xticklabels(alpha)
ax.set_yticklabels(alpha2)
plt.xlabel('ρs [kg/m3]')
plt.ylabel('U [m/s]')
plt.text(40,55,'érosion',color="red",fontsize='17',weight='bold')
plt.text(35,15,"pas d'érosion",color="red",size='17',weight='bold')
plt.legend()
plt.grid()
plt.xlim(0,100)
plt.ylim(0,100)
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\f(u_rhos)2D.png",bbox_inches="tight")
plt.show()


# plt.plot([ustarT0*np.exp(rhoice/rhos0-rhoice/k) for k in y]) #ustarT
# plt.show()

##ustarT(rhos) pour différents rhos0
size=15
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)

plt.figure(figsize=(8,5))
#plt.subplot(2,1,1)
Lrhos0=[150,175,200,225,250,300,350]
for x in Lrhos0:
    Lrho=[k for k in range(x,451)]
    LustarT=[ustarT0*np.exp(rhoice/x-rhoice/y) for y in Lrho]
    if x==300:
        plt.plot(Lrho,LustarT,color='r')
    else:
        plt.plot(Lrho,LustarT,color='b')

size2='12'
plt.text(180, 2.5,"ρ0=150 kg/m3", color='b',fontsize=size2)
plt.text(285, 2.5,"175 kg/m3", color='b',fontsize=size2)
plt.text(345, 2,"200 kg/m3", color='b',fontsize=size2)
plt.text(385, 1.5,"225 kg/m3", color='b',fontsize=size2)
plt.text(385, 1,"250 kg/m3", color='b',fontsize=size2)
plt.text(365, 0.5,"300 kg/m3", color='r',fontsize=size2)
plt.text(415, 0.2,"350 kg/m3", color='b',fontsize=size2)

plt.xlabel("ρs [kg/m3]")
plt.ylabel("u*t [m/s]")
plt.ylim(0,3)
plt.grid()
# plt.subplot(2,1,2)
# plt.plot(date[:40000],ustar[:40000])
# plt.xlabel('temps [ans)]')
# plt.ylabel('u*t [m/s]')
# plt.grid()
plt.text(150,2.7,'(a)',fontsize=size+5)
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\ustarT(rhos)_pour_differents_rhos0.png",bbox_inches="tight")
plt.show()

plt.figure(figsize=(5,5))
plt.plot(mois,ustar[:365*24])
plt.xlabel('temps [mois]')
plt.ylabel('u*t [m/s]')
plt.gca().xaxis.set_ticks(range(2,13,2))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.text(1,1.3,'(b)',fontsize=size+5)
plt.show()