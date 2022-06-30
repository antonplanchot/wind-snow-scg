import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd

data=nc.Dataset("/home/planchoa/Documents/code_amory_5_juin/forcings_20100101-20200101.nc")

#extraire les variables meteo du fichier netcdf
U=data.variables['U2'][:,0,0]   #[m/s]
time=data.variables['time'][:] #[hours since 01/01/2010]
P=data.variables['PRES'][:,0,0] #[hPa]
T=data.variables['T2'][:,0,0]   #[K]
SF=data.variables['RRR'][:,0,0]  #[K]
L=data.variables['LWin'][:,0,0] # #LW...
nb_date=U.shape[0]

##création des variables de temps adequates
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
vk=0.4   #constante de Von Karman   sans dimension

d=0.5    #dendricite
s=0.5   #sphericite

C_D=10**-3   #coefficient de drag pour le moment   sans dimension
z0=10**-4   #longueur de rugosite aerodynamique   m
g=9.81   # acceleration gravitationnelle   m/s2
dt=3600   # pas de temps des variables meteo   s

#calcule u* a partir de U: hypothese du profil logarithmique de vent #on se place donc à 2m ?
ustar=U*vk/np.log(2/z0)

#indice d'erodabilite
iER=0.75*d-0.5*s+0.5

#vitesse de friction seuil standard
ustarT0=(C_D**(0.5))*(np.log(2.868)-np.log(1+iER))/0.085

#hauteur de saltation
hs=0.08436*ustar**(1.27)

#caracteristique neige
Lrhos0=[k for k in range(230,380,10)]  #densite de neige fraiche en surface  kg/m3
Lrhos0ind=[k for k in range(300,380,10)]  #densite de neige fraiche en surface  kg/m3 (n'intervient pas dans les calculs de ustarT et drhos)
Lrhosmax=[k for k in range(350,550+1,25)]   #densite maximum de la neige en surface  kg/m3
rhoice=920   #densite de glace  kg/m3
LtDR=[1,3,6,12,24,48]   #temps caracteristique de compaction de la neige fraiche  s
Lfact=[0.2,0.35,0.6,0.8,1,1.5,2,2.5,3]
Lforcage_vent=[k/100 for k in range(90,140,5)]


Lrhos_mean=[]
LACtot_mean=[]
LACfinal=[]
# Lseuil_AC=[]
Loccurence_ER=[]
LEReff_moy=[]
# seuil=50
SFcumul=[np.sum(SF[:k+1]) for k in range(len(SF))]

## calcul 1D de la masse erodee

for k in range(1): #len(Lrhos0ind)
#choix pour étude de sensibilité (changer 3 choses!)

    choice=33 #0:fact, #1:rhos0, #2:vent , #3: tDR, #4: rhosmax, #5: rhos0ind

    rhos0=Lrhos0[7] #k=7 pour rhos0=300
    rhos0ind=Lrhos0ind[7] #k=7 pour rhos0ind=300
    rhosmax=Lrhosmax[4] #4
    tDR=LtDR[4]*3600 #4
    fo_vent=Lforcage_vent[2] #2
    fact=Lfact[4] #4

    #derivee de la densite de la neige en fonction du temps, pour rhos entre rhos0 et rhosmax
    if choice==5:
        drhos=(rhosmax-rhos0ind)/tDR
    else:
        drhos=(rhosmax-rhos0)/tDR

    SF2=np.copy(SF)*fact
    U1=fo_vent*U
    ustar=U1*vk/np.log(2/z0)

    #initialisation des 3 couches : CF (fraiche) ,CI (intermédiaire) ,CC (ciment)
    rhos=np.array([[np.nan,np.nan,np.nan]]*nb_date,dtype=float)
    AC=np.array([[0,0,0]]*nb_date,dtype=float)
    EReff=np.array([[0,0,0]]*nb_date,dtype=float)

    AC[0,:]=np.array([0.,0.,SF2[0]])
    for i in range(0,nb_date):
        if i>1:
            rhos[i,:]=rhos[i-1,:]

            #contribution des precipitations a l'accumulation nette
            AC[i,:]=AC[i-1,:]+np.array([0.,0.,SF2[i]])
        if AC[i,2]>0:
            if choice!=5:
                rhos[i,2]=rhos0
            else:
                rhos[i,2]=rhos0ind

        ustarT=list(map(lambda x: ustarT0*np.exp(x),(rhoice/rhos0-rhoice/rhos[i])))

        # for m in range(3):
        #     if rhos[i,m]<rhos0:
        #         ustarT[m]=ustarT0

        test=True

        couche=2
        if ustar[i]>ustarT[couche]: #si CF erodable et existante
            qsalt=(ustar[i]**2-ustarT[couche]**2) / (3.25*ustar[i]*g*hs[i])
            #fraction massique de neige dans la couche de saltation
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk    # erosion potentielle
            ER=Ep * P[i]*100/(287*T[i]) * dt # erosion max

            if ER<AC[i,couche]:
                test=False
                rhos[i,couche]+=drhos*dt
                if rhos[i,couche]>rhosmax:
                    rhos[i,couche]=rhosmax
            else:  #toute la couche est erodée
                rhos[i,couche]=np.nan
            EReff[i,couche]=-min(ER,AC[i,couche])  # erosion effective
            AC[i,couche]+=EReff[i,couche]

        couche=1
        if test and ustar[i]>ustarT[couche]: #si CI erodable et existante
            qsalt=(ustar[i]**2-ustarT[couche]**2) / (3.25*ustar[i]*g*hs[i])
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk
            ER=Ep * P[i]*100/(287*T[i]) * dt

            if ER<AC[i,couche]:
                test=False
                rhos[i,couche]+=drhos*dt
                if rhos[i,couche]>rhosmax:
                    rhos[i,couche]=rhosmax
            else:
                rhos[i,couche]=np.nan
            EReff[i,couche]=-min(ER,AC[i,couche])
            AC[i,couche]+=EReff[i,couche]

        couche=0
        if test and ustar[i]>ustarT[couche] and rhos[i,couche]<rhosmax: #si CC erodable et existante
            qsalt=(ustar[i]**2-ustarT[couche]**2) / (3.25*ustar[i]*g*hs[i])
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk
            ER=Ep * P[i]*100/(287*T[i]) * dt

            if ER>=AC[i,couche]:
                rhos[i,couche]=np.nan
            EReff[i,couche]=-min(ER,AC[i,couche])
            AC[i,couche]+=EReff[i,couche]

        #repaquatage des couches
        if rhos[i,1]==rhosmax: #CI-->CC
            rhos[i,1]=np.nan
            AC[i,0]+=AC[i,1]
            AC[i,1]=0
            rhos[i,0]=rhosmax

        if rhos[i,2]==rhosmax and np.isnan(rhos[i,1]): #CF-->CC
            rhos[i,2]=np.nan
            AC[i,0]+=AC[i,2]
            AC[i,2]=0
            rhos[i,0]=rhosmax

        elif rhos[i,2]>rhos0: #CF-->CI
            if np.isnan(rhos[i,1]):
                rhos[i,1]=rhos[i,2]
            else:
                rhos[i,1]=(rhos[i,1]*AC[i,1]+rhos[i,2]*AC[i,2])/(AC[i,1]+AC[i,2])
            rhos[i,2]=np.nan
            AC[i,1]+=AC[i,2]
            AC[i,2]=0


    mrhos = np.ma.masked_array(rhos[:,1],np.isnan(rhos[:,1]))#pour masquer les nan pour le calcul de moyenne
    Lrhos_mean.append(np.mean(mrhos))

    ERefftot=list(np.sum(EReff,axis=1))
    Loccurence_ER.append((len(ERefftot)-ERefftot.count(0))/10)

    LEReff_moy.append(-np.mean(ERefftot))

    ACtot=np.sum(AC,axis=1)
    LACtot_mean.append(np.mean(ACtot))
    LACfinal.append(ACtot[-1])

    # seuilAC=0
    # for k in range(len(ACtot)):
    #     if k>0:
    #         if ACtot[k]>seuil and ACtot[k-1]<seuil:
    #             seuilAC+=1
    # Lseuil_AC.append(seuilAC/10)

        ## plot serie temporelle

    size=13
    parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
    plt.rcParams.update(parameters)

    n=4
    plt.figure(figsize=(12,12))
    string='3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(fact)+', fU='+str(fo_vent)

    plt.subplot(n,1,1)
    #plt.title(string)
    plt.plot(date,SFcumul,color='b',label='acc. brute totale')
    plt.plot(date,ACtot,color='k',label='acc. nette totale')
    plt.ylabel('accumulation totale \n [mmwe]')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(2010,2020,1))
    plt.gca().yaxis.set_ticks(range(0,2001,500))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')


    plt.subplot(n,1,2)
    plt.plot(date,AC[:,2],color='b',label='CF')
    plt.plot(date,AC[:,1],color='y',label='CI')
    plt.plot(date,AC[:,0],color='r',label='CD')
    plt.ylabel('accumulation nette \n [mmwe]')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(2010,2020,1))
    #plt.gca().yaxis.set_ticks(range(0,701,100))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,3)
    plt.plot(date,SF2,color='b',label='Précipitation')
    plt.plot(date,ERefftot,color='r',label='Erosion effective totale')
    plt.ylabel('[mmwe/h]')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(2010,2020,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,4)
    plt.plot(date,rhos[:,1],color='k',label='masse volumique neige CI')#marker=',',linestyle='None',
    plt.ylabel('ρ [kg/m3]')
    plt.xlabel('temps (années)')
    plt.legend()
    plt.xlim(2009.5,2020.5)
    plt.gca().xaxis.set_ticks(range(2010,2020,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.savefig('/home/planchoa/Documents/figures/'+string+'.png',bbox_inches="tight")
    plt.show()
    # plt.pause(1)
    # plt.close()

# size=12
# parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
# plt.rcParams.update(parameters)
# n=2
# plt.figure(figsize=(10,3.5))
#
# Lvar=[[Lfact,"facteur de multiplication des précipitations",'3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF=variable, fU='+str(fo_vent)],[Lrhos0,'ρ0 [kg/m3]','3C_Erlimite_rhos0=variable, rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(fact)+', fU='+str(fo_vent)],[Lforcage_vent,'facteur de multiplication du vent','3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(fact)+', fU=variable'],[LtDR,"temps caractéristique de densification [h]",'3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR=variable, fSF='+str(fact)+', fU='+str(fo_vent)],[Lrhosmax,'ρmax [kg/m3]','3C_Erlimite_rhos0='+str(rhos0)+', rhosmax=variable, tDR='+str(tDR/3600)+'h, fSF='+str(fact)+', fU='+str(fo_vent)],[Lrhos0ind,'ρ0ind [kg/m3]','3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(fact)+', fU='+str(fo_vent)+', rhos0ind=variable']]
# var=Lvar[choice][0]
# strvar=Lvar[choice][1]
# string=Lvar[choice][2]
#
# marker='*'
#
# plt.subplot(1,2,1)
# #plt.title(string)
# plt.plot(var,Lrhos_mean,marker,label='masse volumique CI moyenne ')
# plt.grid()
# plt.ylabel('ρ [kg/m3]')
# plt.xlabel(strvar)
# plt.legend()
#
# plt.subplot(1,2,2)
# #plt.plot(var,LACtot_mean,'k'+marker,label='AC totale moyenne (moyenne temporelle)')
# #plt.plot(var,[x/10 for x in LACfinal],color='orange', marker=marker,linestyle='None',label='AC totale annuelle')
# plt.plot(var,[x/10 for x in LACfinal],color='k', marker=marker,linestyle='None',label='accumulation nette annuelle')
# plt.grid()
# plt.ylabel('accumulation [mm w.e. an-1]')
# plt.xlabel(strvar)
# plt.legend()
#
# # plt.subplot(2,2,3)
# # plt.plot(var,Loccurence_ER,'r'+marker,label="occurences d'épisodes d'érosion")
# # plt.grid()
# # plt.ylabel("occurences/an")
# # plt.xlabel(strvar)
# # plt.legend()
# #
# # plt.subplot(2,2,4)
# # plt.plot(var,LEReff_moy,'r'+marker,label="érosion effective moyenne par occurence")
# # plt.grid()
# # plt.ylabel('mm w.e./occurence')
# # plt.xlabel(strvar)
# # plt.legend()
#
# plt.savefig('/home/planchoa/Documents/figures/'+string+'.png',bbox_inches="tight")
# plt.show()


## sur 1 ans : 2010
size=13
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)

#
# n=4
# plt.figure(figsize=(12,12))
# string='3C_Erlimite_2010_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(fact)
#
# plt.subplot(n,1,1)
# #plt.title(string)
# #plt.plot(mois,ACtot[:365*24],color='k',label='accumulation totale')
# plt.plot(mois,AC[:,2][:365*24],color='b',label='CF')
# plt.plot(mois,AC[:,1][:365*24],color='y',label='CI')
# plt.plot(mois,AC[:,0][:365*24],color='r',label='CD')
# plt.ylabel('accumulation nette \n [mmwe]')
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,14,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
# plt.xticks(color='w')
#
# plt.subplot(n,1,2)
# plt.plot(mois,SF2[:365*24],color='b',label='Précipitation')
# plt.plot(mois,ERefftot[:365*24],color='r',label='Erosion effective totale')
# plt.ylabel('[mmwe/h]')
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,14,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
# plt.xticks(color='w')
#
# plt.subplot(n,1,3)
# plt.plot(mois,rhos[:,1][:365*24],color='k',label='masse volumique neige CI')#marker=',',linestyle='None',
# plt.ylabel('ρ [kg/m3]')
# plt.xlim(0.4,13.6)
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,14,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
# plt.xticks(color='w')
#
# plt.subplot(n,1,4)
# plt.plot(mois,U[:365*24],color='g',label='vitesse du vent U')#marker=',',linestyle='None',
# plt.ylabel('U [m/s]')
# plt.xlabel('temps (mois)')
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,13,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
#
# plt.savefig('/home/planchoa/Documents/figures/'+string+'.png',bbox_inches="tight")
#
# plt.show()


# ##plot direction vent
# # L2=[]
# # ERefftot2=[]
# # pas=24*14+1 #nombre impair [h] (pas de lissage
# # moitpas=pas//2
# # for k in range(moitpas,nb_date-moitpas):
# #     L2.append(np.mean(L[k-moitpas:k+moitpas+1]))
# #     ERefftot2.append(-np.mean(ERefftot[k-moitpas:k+moitpas+1]))
#
#
# # Lpond2=np.copy(Lpond)
# # Lpond2.sort()
#
# # n=3
# # plt.figure(figsize=(15,15))
# # plt.subplot(n,1,1)
# # plt.plot(Lpond)
# # plt.grid()
# # plt.subplot(n,1,2)
# # plt.plot(date,L, label='direction vent')
# # plt.plot(date[moitpas:-moitpas],L2, label='direction vent lissée')
# # plt.grid()
# # plt.subplot(n,1,3)
# # plt.plot(date,[-x for x in ERefftot])
# # plt.plot(date[moitpas:-moitpas],ERefftot2)
# # plt.grid()
# # plt.show()



# S=np.sum(ERefftot[82277+moitpas:82277+4950-moitpas])
# Lpond= [x*y for x,y in zip(WDIRliss,ERefftot[82277+moitpas:82277+4950-moitpas])]
# moy=np.sum(Lpond)/S #moyenne pondérée
#
# step=5
# ran=range(50,285+1,step)
# #ran=range(45,285+1,30)
# J=[k for k in ran]
# J2=[0 for _ in ran]
# compt=0
# for k in range(82277+moitpas,82277+4950-moitpas):
#     i=0
#     test=True
#     while test:
#         if J[i]<=WDIRliss[compt]<J[i+1]:
#             J2[i]+=-ERefftot[k]
#             test=False
#         i+=1
#     compt+=1
#
# for k in range(len(J2)):
#     if J2[k]==0:
#         J2[k]=np.nan
#
# J3=[k+step/2 for k in ran]
# plt.figure(figsize=(10,5))
# #plt.axvline(x=moy,color='gray')
# plt.plot(J3,J2,'.')
# #plt.text(110, 12.5, 'direction moyenne='+str(round(moy))+'°')
# plt.xlabel('direction du vent [°]')
# plt.ylabel('érosion effective [mmwe]')
# plt.title('érosion éolienne par intervalle de '+str(step)+'° de direction du vent du 22/5/2019 au 14/12/2019 (capteur 1, pas de 25h)')
# plt.savefig('/home/planchoa/Documents/figures/erosion(direction_vent).png',bbox_inches="tight")
# plt.show()


# # plt.plot(WDIR,[-x for x in ERefftot[82277:82277+4950]],'.',markersize='1')
# # plt.show()

##quantifier âge neige érodée



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
plt.savefig('/home/planchoa/Documents/figures/age_neige_erodee.png',bbox_inches="tight")


plt.plot(Lage_neige_erode,Lqte_erod,'.',markersize=1)
plt.xlabel('age maximum de la neige érodée [jours]')
plt.ylabel('érosion pour une occurence [mmwe]')
plt.show()
plt.savefig('/home/planchoa/Documents/figures/age_occurence_erosion.png',bbox_inches="tight")


plt.hist(Lage_neige_erode,bins=50)
plt.show()


step=5
ran=range(0,40+1,step)
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

for k in range(len(J3)):
    plt.plot([J3[k],J3[k]],[0,J2[k]],linewidth='5',color='b')
plt.grid()
plt.xlabel('age maximum de la neige érodée [jours]')
plt.ylabel('érosion de 2010 à 2020 [mmwe]')
plt.show()
plt.savefig('/home/planchoa/Documents/figures/erosion(age_neige).png',bbox_inches="tight")




from matplotlib import colors
fig, ax = plt.subplots(figsize=(5, 5), sharex=True, sharey=True,
                        tight_layout=True)


# As well as define normalization of the colors
hist=ax.hist2d(Lage_neige_erode, Lqte_erod, bins=40, norm=colors.LogNorm())
plt.colorbar()
plt.show()

## moyenne des vents en 2010

def month2day(month):
    bisextile=False
    L=[31,28+bisextile,31,30,31,30,31,31,30,31,30,31]
    print(np.sum(L[:month]))


plt.plot(mois,U[:365*24],color='g',label='vitesse du vent U')#marker=',',linestyle='None',
plt.ylabel('U [m/s]')
plt.xlabel('temps (mois)')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,13,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.show()

month2day(6)
month2day(9)

print(np.mean(U[181*24:273*24]))















