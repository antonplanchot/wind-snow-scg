import numpy as np

def erosion1C(T,P,U2,SF2,rhos0,rhos0ind,rhosmax,tDR):
    '''
    Calcul 1D de la masse erodee pour 1 couche au cours du temps
    Parameters
    ----------
    T : numpy array (len,)
        temperature [K]
    P : numpy array (len,)
        pression [hPa]
    U2 : numpy array (len,)
        vitesse du vent à 2m [m/s]
    SF2 : numpy array (len,)
        précipitation horaire [mmwe]
    rhos0 : float
        masse volumique de la neige fraiche en surface [kg/m3] (valeur standard:300)
    rhos0ind : float
        masse volumique de la neige fraiche en surface n'intervennant pas dans calcul de ustarT [kg/m3]
    rhosmax : float
        masse volumique maximum de la neige [kg/m3] (valeur standard:450)
    tDR : float
        temps caracteristique de densification par érosion de la neige [h] (valeur standard:24)

    Returns
    -------
    rhos : list (len)
        masse volumique des 3 couches au cours du temps [kg/m3]
    AC : list (len)
        épaisseure des 3 couches au cours du temps [mmwe]
    EReff : list (len)
        érosion subie par les 3 couches au cours du temps [kg/m3]
    '''
    #constantes
    vk=0.4   #constante de Von Karman   [sans dimension]
    C_D=10**-3   #coefficient de drag pour le moment   [sans dimension]
    z0=10**-4   #longueur de rugosite aerodynamique   [m]
    g=9.81   # acceleration gravitationnelle   [m/s2]
    dt=3600   # pas de temps des variables meteo   [s]
    
    d=0.5    #dendricite  [sans dimension]
    s=0.5   #sphericite  [sans dimension]
    iER=0.75*d-0.5*s+0.5   #indice d'erodabilite  [sans dimension]
    #vitesse standard de friction seuil [m/s]
    ustarT0=(C_D**(0.5))*(np.log(2.868)-np.log(1+iER))/0.085
    
    rhoice=920  #masse volumique de la glace  [kg/m3]

    nb_date=len(U2)
    ustar=U2*vk/np.log(2/z0)    #vitesse de friction seuil [m/s]

    if rhos0ind==0: #pas de test de sensibilité à rhos0ind
        rhos02=rhos0
    else:
        rhos02=rhos0ind
    #taux de densification de la neige
    drhos=(rhosmax-rhos02)/tDR

    #initialisation de Acc
    Acc=0
    rhos=[rhos0]*(SF2[0]>0)+[np.nan]*(SF2[0]==0)

    AC,ER,EReff=[],[],[]
    for i in range(0,nb_date):

        #contribution  des precipitations a l'accumulation nette
        AC.append(Acc+SF2[i])
        #rétroaction négative de l'érosion sur elle meme
        if i>0:
            if Acc==0 and SF2[i]==0:#pas de neige
                rhos.append(np.nan)
            elif Acc==0 and SF2[i]>0:#que neige fraiche
                rhos.append(rhos0)
            elif Acc>0 and SF2[i]==0:#que neige pas fraiche
                rhos.append(rhos[-1])
            else:
                rhos.append((Acc*rhos[-1]+SF2[i]*rhos02)/(AC[-1]))

        ustarT=ustarT0*np.exp(rhoice/rhos0-rhoice/rhos[-1]) #nan si rhos[-1]=nan

        # si erosion (et couche non vide)
        if ustar[i]>ustarT and rhos[-1]<rhosmax:
            qsalt=(ustar[i]**2-ustarT**2) / (3.25*g*0.08436*ustar[i]**2.27) #kg/kg 
            #fraction massique de neige dans la couche de saltation m/s
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk   #erosion potentielle m/s kg/kg
            ER.append(Ep * P[i]*100/(287*T[i]) * dt) # erosion max 
            EReff.append(-min(ER[i],AC[i]))   # erosion effective
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

    return rhos,AC,EReff