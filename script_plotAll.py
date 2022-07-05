import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from script_erosion3C_1 import erosion3C
from script_erosion1C_1 import erosion1C
from script_plot1C import plotplurian1C, plotsensibilite1C, plotannuel1C
from script_plot3C import plotplurian3C, plotsensibilite3C, plotannuel3C

'''
plots:
-pluriannuelle: ## plot serie temporelle pluriannuelle
-sensibilité:
'''


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

##preparation des plots
LfactU=[k/100 for k in range(0,180,5)] #facteur multiplicatif de la vitesse du vent
LfactSF=[0.2,0.35,0.6,0.8,1,1.5,2,2.5,3] #facteur multiplicatif des précipitations
Lrhos0=[k for k in range(230,380,10)]  #densite de neige fraiche en surface  kg/m3
Lrhos0ind=[0]+[k for k in range(300,380,10)]  #densite de neige fraiche en surface n'intervennant pas dans calcul de ustarT kg/m3
Lrhosmax=[k for k in range(350,550+1,25)]   #densite maximum de la neige en surface  kg/m3
LtDR=[1,3,6,12,24,48]   #temps caracteristique de compaction de la neige fraiche

Lsensibilite=[[True],LfactU,LfactSF,Lrhos0,Lrhos0ind,Lrhosmax,LtDR]
        

#______________________________________________________________________________


## plot serie temporelle pluriannuelle
#choisir le modèle à utiliser
nb_couche=3 #1 ou 3

Lplurian=[[erosion1C,plotplurian1C],[erosion3C,plotplurian3C]]

#paramètre d'entrée à faire varier
choice=0 #aucun:0 #factU:1 #factSF:2 #rhos0:3 #rhos0ind:4 #rhosmax:5 #tDR:6

size=13
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)
for k in range(len(Lsensibilite[choice])-(choice==4)):
    factU=LfactU[(k-LfactU.index(1))*(choice==1)+LfactU.index(1)]
    factSF=LfactSF[(k-LfactSF.index(1))*(choice==2)+LfactSF.index(1)] 
    rhos0=Lrhos0[(k-Lrhos0.index(300))*(choice==3)+Lrhos0.index(300)]
    rhos0ind=Lrhos0ind[(k+1)*(choice==4)]
    rhosmax=Lrhosmax[(k-Lrhosmax.index(450))*(choice==5)+Lrhosmax.index(450)]
    tDR=LtDR[(k-LtDR.index(24))*(choice==6)+LtDR.index(24)]*3600

    U2=U*factU
    SF2=SF*factSF
    rhos,AC,EReff=Lplurian[nb_couche//2][0](T,P,U2,SF2,rhos0,rhos0ind,rhosmax,tDR)
    string='_rhos0='+str(rhos0)+'_rhosmax='+str(rhosmax)+'_tDR='+str(tDR/3600)+'h_fSF='+str(factSF)+'_fU='+str(factU)
    Lplurian[nb_couche//2][1](choice,date,string,SF2,rhos,AC,EReff)


#______________________________________________________________________________


##plot étude de sensibilité 3C
#paramètre d'entrée à faire varier
choice=1 #factU:1 #factSF:2 #rhos0:3 #rhos0ind:4 #rhosmax:5 #tDR:6

Lrhos_mean=[]
#LACtot_mean=[]
LACfinal=[]
Loccurence_ER=[]
LEReff_moy=[]
for k in range(len(Lsensibilite[choice])-(choice==4)):
    factU=LfactU[(k-LfactU.index(1))*(choice==1)+LfactU.index(1)]
    factSF=LfactSF[(k-LfactSF.index(1))*(choice==2)+LfactSF.index(1)] 
    rhos0=Lrhos0[(k-Lrhos0.index(300))*(choice==3)+Lrhos0.index(300)]
    rhos0ind=Lrhos0ind[(k+1)*(choice==4)]
    rhosmax=Lrhosmax[(k-Lrhosmax.index(450))*(choice==5)+Lrhosmax.index(450)]
    tDR=LtDR[(k-LtDR.index(24))*(choice==6)+LtDR.index(24)]*3600

    U2=U*factU
    SF2=SF*factSF
    rhos,AC,EReff=erosion3C(T,P,U2,SF2,rhos0,rhos0ind,rhosmax,tDR)
   
    mrhos = np.ma.masked_array(rhos[:,1],np.isnan(rhos[:,1]))#pour masquer les nan pour le calcul de moyenne
    Lrhos_mean.append(np.mean(mrhos))

    ACtot=np.sum(AC,axis=1)
    #LACtot_mean.append(np.mean(ACtot))
    LACfinal.append(ACtot[-1])
    ERefftot=list(np.sum(EReff,axis=1))
    nb_ER=len(ERefftot)-ERefftot.count(0)
    Loccurence_ER.append(nb_ER/10)
    LEReff_moy.append(-np.sum(ERefftot)/nb_ER)

Lvar=[[LfactU,'facteur de multiplication du vent','3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU=variable'],[LfactSF,"facteur de multiplication des précipitations",'3C_erl_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF=variable, fU='+str(factU)],[Lrhos0,'ρ0 [kg/m3]','3C_Erlimite_rhos0=variable, rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU='+str(factU)],[Lrhos0ind,'ρ0ind [kg/m3]','3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU='+str(factU)+', rhos0ind=variable'],[LtDR,"temps caractéristique de densification [h]",'3C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR=variable, fSF='+str(factSF)+', fU='+str(factU)],[Lrhosmax,'ρmax [kg/m3]','3C_Erlimite_rhos0='+str(rhos0)+', rhosmax=variable, tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU='+str(factU)]]
var=Lvar[choice-1][0]
strvar=Lvar[choice-1][1]
string=Lvar[choice-1][2]
plotsensibilite3C(var,strvar,string,Lrhos_mean,LACfinal,Loccurence_ER,LEReff_moy)



##plot étude de sensibilité 1C
#paramètre d'entrée à faire varier
choice=1 #factU:1 #factSF:2 #rhos0:3 #rhos0ind:4 #rhosmax:5 #tDR:6

Lrhos_mean=[]
#LACtot_mean=[]
LACfinal=[]
Loccurence_ER=[]
LEReff_moy=[]
for k in range(len(Lsensibilite[choice])-(choice==4)):
    factU=LfactU[(k-LfactU.index(1))*(choice==1)+LfactU.index(1)]
    factSF=LfactSF[(k-LfactSF.index(1))*(choice==2)+LfactSF.index(1)] 
    rhos0=Lrhos0[(k-Lrhos0.index(300))*(choice==3)+Lrhos0.index(300)]
    rhos0ind=Lrhos0ind[(k+1)*(choice==4)]
    rhosmax=Lrhosmax[(k-Lrhosmax.index(450))*(choice==5)+Lrhosmax.index(450)]
    tDR=LtDR[(k-LtDR.index(24))*(choice==6)+LtDR.index(24)]*3600

    U2=U*factU
    SF2=SF*factSF
    rhos,AC,EReff=erosion3C(T,P,U2,SF2,rhos0,rhos0ind,rhosmax,tDR)
   
    mrhos = np.ma.masked_array(rhos,np.isnan(rhos))#pour masquer les nan pour le calcul de moyenne
    Lrhos_mean.append(np.mean(mrhos))

    #LACtot_mean.append(np.mean(ACtot))
    LACfinal.append(AC[-1])
    
    ERefftot=list(np.sum(EReff,axis=1))
    nb_ER=len(EReff)-EReff.count(0)
    Loccurence_ER.append(nb_ER/10)
    LEReff_moy.append(-np.sum(EReff)/nb_ER)

Lvar=[[LfactU,'facteur de multiplication du vent','1C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU=variable'],[LfactSF,"facteur de multiplication des précipitations",'1C_erl_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF=variable, fU='+str(factU)],[Lrhos0,'ρ0 [kg/m3]','1C_Erlimite_rhos0=variable, rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU='+str(factU)],[Lrhos0ind,'ρ0ind [kg/m3]','1C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU='+str(factU)+', rhos0ind=variable'],[LtDR,"temps caractéristique de densification [h]",'1C_Erlimite_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR=variable, fSF='+str(factSF)+', fU='+str(factU)],[Lrhosmax,'ρmax [kg/m3]','1C_Erlimite_rhos0='+str(rhos0)+', rhosmax=variable, tDR='+str(tDR/3600)+'h, fSF='+str(factSF)+', fU='+str(factU)]]
var=Lvar[choice-1][0]
strvar=Lvar[choice-1][1]
string=Lvar[choice-1][2]
plotsensibilite1C(var,strvar,string,Lrhos_mean,LACfinal,Loccurence_ER,LEReff_moy)

#______________________________________________________________________________

## plot serie temporelle annuelle (2010)
#choisir le modèle à utiliser
nb_couche=3 #1 ou 3

Lan=[[erosion1C,plotannuel1C],[erosion3C,plotannuel3C]]

#paramètre d'entrée à faire varier
choice=0 #aucun:0 #factU:1 #factSF:2 #rhos0:3 #rhos0ind:4 #rhosmax:5 #tDR:6

size=13
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)
for k in range(len(Lsensibilite[choice])-(choice==4)):
    factU=LfactU[(k-LfactU.index(1))*(choice==1)+LfactU.index(1)]
    factSF=LfactSF[(k-LfactSF.index(1))*(choice==2)+LfactSF.index(1)] 
    rhos0=Lrhos0[(k-Lrhos0.index(300))*(choice==3)+Lrhos0.index(300)]
    rhos0ind=Lrhos0ind[(k+1)*(choice==4)]
    rhosmax=Lrhosmax[(k-Lrhosmax.index(450))*(choice==5)+Lrhosmax.index(450)]
    tDR=LtDR[(k-LtDR.index(24))*(choice==6)+LtDR.index(24)]*3600

    U2=U*factU
    SF2=SF*factSF
    rhos,AC,EReff=Lan[nb_couche//2][0](T,P,U2,SF2,rhos0,rhos0ind,rhosmax,tDR)
    string='_rhos0='+str(rhos0)+'_rhosmax='+str(rhosmax)+'_tDR='+str(tDR/3600)+'h_fSF='+str(factSF)+'_fU='+str(factU)
    Lan[nb_couche//2][1](choice,mois,string,U2,SF2,rhos,AC,EReff)