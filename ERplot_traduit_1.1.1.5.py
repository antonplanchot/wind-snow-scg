trigger=True
if trigger:
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
    fact=1
    SF=data.variables['RRR'][:,0,0]*fact  #[K] #WARNING coefficient multiplicateur a adapter (5 de base)
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

for k in range(1):
    #caracteristique neige
    Lrhos0=[150,200,225,250,275,300,350,400,450]  #densite de neige fraiche en surface  kg/m3
    rhos0=Lrhos0[5]#5
    Lrhosmax=[400,450,500,550]   #densite maximum de la neige en surface  kg/m3
    rhosmax=Lrhosmax[1]
    rhoice=920   #densite de glace  kg/m3
    LtDR=[1,3,6,12,24,48]   #temps caracteristique de compaction de la neige fraiche  s
    tDR=LtDR[4]*3600
    Lforcage_vent=[k/100 for k in range(90,140,5)]
    fo_vent=Lforcage_vent[2] #2

    #derivee de la densite de la neige en fonction du temps, pour rhos entre rhos0 et rhosmax
    drhos=(rhosmax-rhos0)/tDR

    ## calcul 1D de la masse erodee

    #initialisation de Acc
    Acc=0
    rhos=[rhos0]

    AC,ER,EReff=[],[],[]
    for i in range(0,nb_date):

        #contribution  des precipitations a l'accumulation nette
        AC.append(Acc+SF[i])
        #rétroaction négative de l'érosion sur elle meme
        if i>0:
            if Acc==0 and SF[i]==0:#pas de neige
                rhos.append(np.nan)
            elif Acc==0 and SF[i]>0:#que neige fraiche
                rhos.append(rhos0)
            elif Acc>0 and SF[i]==0:#que neige pas fraiche
                rhos.append(rhos[-1])
            else:
                rhos.append((Acc*rhos[-1]+SF[i]*rhos0)/(AC[-1]))

        ustarT=ustarT0*np.exp(rhoice/rhos0-rhoice/rhos[-1])

        # si erosion
        if ustar[i]>ustarT and not np.isnan(rhos[-1]) and rhos[-1]<rhosmax:
            qsalt=(ustar[i]**2-ustarT**2) / (3.25*ustar[i]*g*hs[i])   #kg/kg #cf article MAR
            #fraction massique de neige dans la couche de saltation


            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk   #m/s kg/kg #??
            # erosion potentielle


            ER.append(Ep * P[i]*100/(287*T[i]) * dt)   #m/s * kg/kg * kg/m3 * s = kg/m2 #Ep*
            # erosion max = accu nette #??


            EReff.append(-min(ER[i],AC[i]))   # erosion effective (on met un "-" pour le plot...)
            #on prend ER sauf quand il n'y a pas assez de neige accumulée.

            test=i
        # si pas erosion
        else:
            ER.append(0)
            EReff.append(0)

        AC[-1] += EReff[-1]   # accumulation nette ponctuelle
        Acc=AC[-1]

        if not np.isnan(rhos[-1]) and Acc>0 and rhos[-1]<rhosmax and test==i:
            #neige presente, pas densifiée au max, vient de subir érosion
            rhos[-1]+=drhos*dt
            if rhos[-1]>rhosmax:
                rhos[-1]=rhosmax

    ## plot serie temporelle

    size=13
    parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
    plt.rcParams.update(parameters)

    n=3
    plt.figure(figsize=(9,12))

    string='1C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(fact)+', fU='+str(fo_vent)
    plt.subplot(n,1,1)
    #plt.title(string)
    plt.plot(date,AC,color='k',label='Accumulation de neige')
    plt.ylabel('Accumulation nette \n [mmwe]')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(2010,2020+1,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,2)
    plt.plot(date,SF,color='b',label='Précipitation')
    plt.plot(date,EReff,color='r',label='Erosion effective')
    plt.ylabel('[mmwe/h]')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(2010,2020+1,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,3)
    plt.plot(date,rhos,color='k',label='masse volumique neige')#marker=',',linestyle='None',
    plt.ylabel('ρ [kg/m3]')
    plt.xlabel('temps (années)')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(2010,2020+1,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.savefig('/home/planchoa/Documents/figures/'+string+'.png',bbox_inches="tight")

    plt.show()
    # plt.pause(1.5)
    # plt.close()




#
# n=4
# plt.figure(figsize=(15,15))
#
# plt.subplot(n,1,1)
# plt.title('année 2010: rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fact_prec='+str(fact))
# plt.plot(mois,AC[:365*24],color='k',label='Net accumulation')
# plt.ylabel('[mmwe]')
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,14,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
#
# plt.subplot(n,1,2)
# plt.plot(mois,SF[:365*24],color='b',label='Precipitation')
# plt.plot(mois,EReff[:365*24],color='r',label='Effective erosion')
# plt.ylabel('[mmwe/h]')
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,14,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
#
# plt.subplot(n,1,3)
# plt.plot(mois,rhos[:365*24],color='k',label='densite moyenne de la neige')#marker=',',linestyle='None',
# plt.ylabel('[kg/m3]')
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,14,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
#
# plt.subplot(n,1,4)
# plt.plot(mois,U[:365*24],color='g',label='vitesse du vent')#marker=',',linestyle='None',
# plt.ylabel('[m/s]')
# plt.legend()
# plt.gca().xaxis.set_ticks(range(1,14,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
#
# plt.savefig('/home/planchoa/Documents/1couche_2010_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fact_prec='+str(fact)+'.png')
#
# plt.show()


# ## double plot
# fig, ax1 = plt.subplots()
# color = 'tab:blue'
# ax1.set_xlabel('time (day)')
# ax1.set_ylabel('rhos', color=color)
# ax1.plot(date, rhos,color=color)
# ax1.tick_params(axis='y', labelcolor=color)
#
# ax2 = ax1.twinx()
#
# color = 'tab:red'
# ax2.set_ylabel('AC', color=color)
# ax2.plot(date, AC, color=color)
# ax2.tick_params(axis='y', labelcolor=color)
#
# plt.gca().xaxis.set_ticks(range(2010,2020,1))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
# fig.tight_layout()
# plt.show()








