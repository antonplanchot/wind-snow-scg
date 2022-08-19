import numpy as np

def erosion4C(T,P,U2,SF2,rhos0,rhos0ind,rhosmax,tDR):
    '''
    Calcul 1D de la masse erodee pour 4 couches au cours du temps
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
    rhos : numpy array (len,4)
        masse volumique des 4 couches au cours du temps [kg/m3]
    AC : numpy array (len,4)
        épaisseure des 4 couches au cours du temps [mmwe]
    EReff : numpy array (len,4)
        érosion subie par les 4 couches au cours du temps [kg/m3]
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

    #initialisation des 3 couches : CF (fraiche) ,CI (intermédiaire) ,CC (ciment)
    rhos=np.array([[np.nan,np.nan,np.nan,np.nan]]*nb_date,dtype=float)
    AC=np.array([[0,0,0,0]]*nb_date,dtype=float)
    EReff=np.array([[0,0,0,0]]*nb_date,dtype=float)

    AC[0,:]=np.array([0.,0.,0.,SF2[0]])
    for i in range(nb_date):
        if i>1:
            rhos[i,:]=rhos[i-1,:]

            #contribution des precipitations a l'accumulation brute
            AC[i,:]=AC[i-1,:]+np.array([0.,0.,0.,SF2[i]])
        if AC[i,3]>0:
            rhos[i,3]=rhos02

        ustarT=list(map(lambda x: ustarT0*np.exp(x),(rhoice/rhos0-rhoice/rhos[i])))

        test=True
        
        couche=3
        fracdtF=0
        if ustar[i]>ustarT[couche]: #si CF erodable et existante
            qsalt=(ustar[i]**2-ustarT[couche]**2) / (3.25*g*0.08436*ustar[i]**2.27)
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
                fracdtF=AC[i,couche]/ER #proportion de dt écoulée
            EReff[i,couche]=-min(ER,AC[i,couche])  # erosion effective
            AC[i,couche]+=EReff[i,couche]
        
        couche=2
        fracdtFI=0
        if test and ustar[i]>ustarT[couche]: #si CFI erodable et existante
            qsalt=(ustar[i]**2-ustarT[couche]**2) / (3.25*g*0.08436*ustar[i]**2.27)
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk
            ER=Ep * P[i]*100/(287*T[i]) * (1-fracdtF) * dt 
            #on réduit le temps d'érosion si la couche précédente a été érodée

            if ER<AC[i,couche]:
                test=False
                rhos[i,couche]+=drhos*dt
                if rhos[i,couche]>rhosmax:
                    rhos[i,couche]=rhosmax
            else:
                rhos[i,couche]=np.nan
                fracdtFI=AC[i,couche]/ER
            EReff[i,couche]=-min(ER,AC[i,couche])
            AC[i,couche]+=EReff[i,couche]
            
        couche=1
        fracdtI=0
        if test and ustar[i]>ustarT[couche]: #si CI erodable et existante
            qsalt=(ustar[i]**2-ustarT[couche]**2) / (3.25*g*0.08436*ustar[i]**2.27)
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk
            ER=Ep * P[i]*100/(287*T[i]) * (1-fracdtF) * (1-fracdtFI) * dt 
            #on réduit le temps d'érosion si la couche précédente a été érodée

            if ER<AC[i,couche]:
                test=False
                rhos[i,couche]+=drhos*dt
                if rhos[i,couche]>rhosmax:
                    rhos[i,couche]=rhosmax
            else:
                rhos[i,couche]=np.nan
                fracdtI=AC[i,couche]/ER
            EReff[i,couche]=-min(ER,AC[i,couche])
            AC[i,couche]+=EReff[i,couche]

        couche=0
        if test and ustar[i]>ustarT[couche] and rhos[i,couche]<rhosmax: #si CC erodable et existante
            qsalt=(ustar[i]**2-ustarT[couche]**2) / (3.25*g*0.08436*ustar[i]**2.27)
            Ep = qsalt * C_D * ustar[i] * np.log(10/z0) / vk
            ER=Ep * P[i]*100/(287*T[i]) * (1-fracdtF) * (1-fracdtFI) * (1-fracdtI) * dt

            if ER>=AC[i,couche]:
                rhos[i,couche]=np.nan
            EReff[i,couche]=-min(ER,AC[i,couche])
            AC[i,couche]+=EReff[i,couche]

        #réorganisation des couches
        if rhos[i,1]==rhosmax: #CI-->CC
            rhos[i,1]=np.nan
            AC[i,0]+=AC[i,1]
            AC[i,1]=0
            rhos[i,0]=rhosmax
        #si cette opération est effectuée, CI a été densifiée donc CF et CFI sont vides
        
        else:
            if rhos[i,3]==rhosmax and np.isnan(rhos[i,1]) and np.isnan(rhos[i,2]): #CF-->CC
            #(CI et CFI sont vides)
                rhos[i,3]=np.nan
                AC[i,0]+=AC[i,3]
                AC[i,3]=0
                rhos[i,0]=rhosmax
    
            else:
                if rhos[i,3]>rhos02: 
                    if np.isnan(rhos[i,1]): #CF-->CI 
                    #CI est vide (donc CFI aussi)
                        rhos[i,1]=rhos[i,3]
                        AC[i,1]+=AC[i,3]
                    else:
                        if np.isnan(rhos[i,2]): #CF-->CFI
                        #CI est non vide, CFI est vide
                            rhos[i,2]=rhos[i,3]
                        else: #CF-->CFI (CFI alimenté que si CI non vide)
                            rhos[i,2]=(rhos[i,2]*AC[i,2]+rhos[i,3]*AC[i,3])/(AC[i,2]+AC[i,3])
                        AC[i,2]+=AC[i,3]
                    rhos[i,3]=np.nan
                    AC[i,3]=0
                    
                if rhos[i,2]>rhos[i,1]: #CFI-->CI 
                #(CFI a été densifiée, par érosion ou par fusion avec CF)
                    rhos[i,1]=(rhos[i,1]*AC[i,1]+rhos[i,2]*AC[i,2])/(AC[i,1]+AC[i,2])
                    rhos[i,2]=np.nan
                    AC[i,1]+=AC[i,2]
                    AC[i,2]=0

    return rhos,AC,EReff