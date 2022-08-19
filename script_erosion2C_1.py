import numpy as np

def erosion2C(T,P,U2,SF2,rhos0,rhos0ind,rhosmax,tDR):
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
        masse volumique de la couche au cours du temps [kg/m3]
    AC : list (len)
        épaisseure de la couche au cours du temps [mmwe]
    EReff : list (len)
        érosion subie par la couche au cours du temps [kg/m3]
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
    rhos=[rhos02]*(SF2[0]>0)+[np.nan]*(SF2[0]==0)
    
    AC=np.array([[0,0]]*nb_date,dtype=float)
    EReff=[]
    
    AC[0,:]=np.array([0.,SF2[0]])
    for i in range(0,nb_date):

        #rétroaction négative de l'érosion sur elle meme
        if i>0:
            #contribution  des precipitations a l'accumulation nette
            AC[i]=np.array([AC[i-1,0],AC[i-1,1]+SF2[i]])
            if AC[i-1,1]==0 and SF2[i]==0:#pas de neige
                rhos.append(np.nan)
            elif AC[i-1,1]==0 and SF2[i]>0:#que neige fraiche
                rhos.append(rhos02)
            elif AC[i-1,1]>0 and SF2[i]==0:#que neige pas fraiche
                rhos.append(rhos[-1])
            else:
                rhos.append((AC[i-1,1]*rhos[-1]+SF2[i]*rhos02)/(AC[i,1]))

        ustarT=ustarT0*np.exp(rhoice/rhos0-rhoice/rhos[-1]) #nan si rhos[-1]=nan
        test=False
        # si erosion (et couche non vide)
        if ustar[i]>ustarT and rhos[-1]<rhosmax:
            qsalt=(ustar[i]**2-ustarT**2) / (3.25*g*0.08436*ustar[i]**2.27) #kg/kg 
            #fraction massique de neige dans la couche de saltation m/s
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk   #erosion potentielle m/s kg/kg
            ER=Ep * P[i]*100/(287*T[i]) * dt # erosion max 
            EReff.append(-min(ER,AC[i,1]))   # erosion effective
            #on prend ER sauf quand il n'y a pas assez de neige accumulée.
            test=True
        # si pas erosion
        else:
            EReff.append(0)

        AC[i,1] += EReff[-1]   # accumulation nette ponctuelle

        if AC[i,1]>0 and test:
            #neige presente, vient de subir érosion
            rhos[i]+=drhos*dt
            if rhos[i]>=rhosmax:
                AC[i,0]+=AC[i][1]
                AC[i,1]=0
                rhos[-1]=np.nan

    return rhos,AC,EReff