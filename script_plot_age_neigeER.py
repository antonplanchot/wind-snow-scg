import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import netCDF4 as nc
from script_erosion3C_2 import erosion3C

##charger fichier netcdf et extraire les variables meteo
data=nc.Dataset("C:/Users/Toshiba/Documents/stage L3/Documents/code_amory_5_juin/forcings_20100101-20200101.nc")

T=data.variables['T2'][:,0,0]   #[K]
P=data.variables['PRES'][:,0,0] #[hPa]
U=data.variables['U2'][:,0,0]   #[m/s]
SF=data.variables['RRR'][:,0,0]  #[mmwe]

##création des variables de temps
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

#erosion sur periode 2010_2020
rhos0=300
rhos0ind=0
rhosmax=450
tDR=24*3600
rhos,AC,EReff=erosion3C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
ACtot=np.sum(AC,axis=1)
ERefftot=list(np.sum(EReff,axis=1))

##série temporelle des ages maximum et quantité d'érosion brute
nb_date=len(U)
Lage_neige_erode=[]
Lqte_erod=[]
datebis=[]
print(ACtot[0])
for k in range(1,nb_date):
    if ACtot[k]<ACtot[k-1]:
        i=1
        test=True
        while i<k-1 and test:
            if ACtot[k-1-i]<=ACtot[k]:
                Lage_neige_erode.append(i/24)
                datebis.append(date[k])
                test=False
                Lqte_erod.append(ACtot[k-1]-ACtot[k])
            i+=1

#plot série temporelle age d'erosion
plt.figure(figsize=(12,12))
plt.subplot(2,1,1)
plt.plot(datebis, Lage_neige_erode,'*')
plt.ylabel('age maximum de la neige érodée [jours]')
plt.grid()
plt.subplot(2,1,2)
plt.plot(date,ACtot)
plt.ylabel('accumulation totale [mmwe]')
plt.grid()
plt.show()
plt.savefig("C:/Users/Toshiba/Documents/stage L3/Documents/figures/age_neige_erodee.png",bbox_inches="tight")

#plot quantité d'erosion en fonction de l'âge maximum de la neige errodée
plt.plot(Lage_neige_erode,Lqte_erod,'.',markersize=1)
plt.xlabel('age maximum de la neige érodée [jours]')
plt.ylabel('érosion pour une occurence [mmwe]')
plt.show()
plt.savefig("C:/Users/Toshiba/Documents/stage L3/Documents/figures/age_occurence_erosion.png",bbox_inches="tight")

plt.hist2d(Lage_neige_erode, Lqte_erod, bins=40, norm=colors.LogNorm())
plt.colorbar(label="nombre d'épisodes d'érosion")
plt.xlabel('âge maximum de la neige érodée [jours]')
plt.ylabel('érosion pour une occurence [mm w.e.]')
plt.show()

#plot histo nombre d'épisode d'érosion en fonction de l'âge maximum de la neige errodée
plt.hist(Lage_neige_erode,bins=50)
plt.xlabel('age maximum de la neige érodée [jours]')
plt.ylabel("nombre d'occurence d'érosion")
plt.show()

#calcul de la quantité d'érosion par intervalle d'âge maximum de la neige errodée
step=5
ran=range(0,45+1,step)
J=[k for k in ran]
J2=[0 for _ in ran]
J3=[k+step/2 for k in ran]
compt=0
for k in range(len(Lage_neige_erode)):
    i=0
    test=True
    while test:
        if J[i]<=Lage_neige_erode[compt]<J[i+1]:
            J2[i]+=Lqte_erod[k]
            test=False
        i+=1
    compt+=1

#plot histo quantité d'erosion en fonction de l'âge maximum de la neige errodée
for k in range(len(J3)):
    plt.plot([J3[k],J3[k]],[0,J2[k]],linewidth='5',color='b')
plt.grid()
plt.xlabel('age maximum de la neige érodée [jours]')
plt.ylabel('érosion de 2010 à 2020 [mmwe]')
plt.show()
plt.savefig("C:/Users/Toshiba/Documents/stage L3/Documents/figures/erosion(age_neige).png",bbox_inches="tight")
S=np.sum(J2)
print(100*J2[0]/S,100*(J2[0]+J2[1])/S)
