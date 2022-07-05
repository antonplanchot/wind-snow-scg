import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd
from datetime import datetime
from script_erosion3C_1 import erosion3C
'''
3. Negative bias on WS_MAX, WS_AVG, WS_MAX_2, WS_AVG_2 at times during the monsoon season due to
rime ice accretion on the propellers.
4. WDIR_2 offset by 180 degrees due to apparent internal wiring error.
5. Erroneous data from wind sensor 1 (WS_AVG, WS_MAX, WDIR) beginning 15 DECEMBER 2019 due to sensor and/or mount
damage/failure.
6. No data from wind sensor 2 (WS_AVG_2, WS_MAX_2, WDIR_2) beginning 5 JANUARY 2020 due to either sensor or cable
failure.
7. No data from 1 SEPTEMBER 2020 through 7 OCTOBER 2020 due to low battery voltage presumably caused by burial of solar
panels by deep snow. Intermittent data from 7 OCTOBER 2020 through 11 OCTOBER 2020 as snow began to ablate.'''

##chargement des fichiers et variables
dataERA5=nc.Dataset("C:/Users/Toshiba/Documents/stage L3/Documents/code_amory_5_juin/forcings_20100101-20200101.nc")
T=dataERA5.variables['T2'][:,0,0]   #[K]
P=dataERA5.variables['PRES'][:,0,0] #[hPa]
U=dataERA5.variables['U2'][:,0,0]   #[m/s]
SF=dataERA5.variables['RRR'][:,0,0]  #[mmwe]

dataAWS=pd.read_csv("C:/Users/Toshiba/Documents/stage L3/Documents/code_amory_5_juin/South_Col_20210630.csv")
dataAWS=np.array(dataAWS)

def conv_date(time):
    ti= nc.num2date(time[:], units = 'hours since 2010-01-01 01:00:00.0')
    date_conv=[da.strftime('%Y/%m/%d/%H') for da in ti]
    return(date_conv)
#mesure d'indice 82277 : '5/22/2019 6:00'
print(conv_date([82277]))
m=4950
#mesure d'indice 82277 + 4950 : '12/14/2019 12:00'
print(conv_date([82277+m]))

WS_MAX=dataAWS[:,4][:m]
WS_AVG=dataAWS[:,5][:m]
WDIR=dataAWS[:,6][:m]
WS_MAX_2=dataAWS[:,7][:m]
WS_AVG_2=dataAWS[:,8][:m]
WDIR_2=dataAWS[:,9][:m]
date=dataAWS[:,0]

#prise en compt biais de mesure capteur2
for k in range(len(WDIR_2)):
    WDIR_2[k]+=180
    if WDIR_2[k]>=360:
        WDIR_2[k]-=360

Ltime=[datetime.strptime(x, '%m/%d/%Y %H:%M') for x in date][:m]
#Ltime = matplotlib.dates.date2num(Ltime)

#moyenne glissante sur les variables mesurées
WDIRliss,WDIR_2liss=[],[]
WS_MAXliss,WS_AVGliss=[],[]
pas=24*3+1 #nombre impair [h] (pas de lissage)
moitpas=pas//2
for k in range(moitpas,len(WDIR)-moitpas):
    WDIRliss.append(np.mean( WDIR[k-moitpas:k+moitpas+1]))
    WDIR_2liss.append(np.mean( WDIR_2[k-moitpas:k+moitpas+1]))
    WS_MAXliss.append(np.mean(  WS_MAX[k-moitpas:k+moitpas+1]))
    WS_AVGliss.append(np.mean(  WS_AVG[k-moitpas:k+moitpas+1]))

#suppresion des dernières valeurs de direction du vent du capteur 2
WDIR_2liss[-310:]=310*[np.nan]

##plot caractéristiques mesures météo
nb=2
size=13
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)
plt.figure(figsize=(11,9))
# plt.subplot(nb,1,1)
# plt.title("mesures de vents sur station automatique au col sud de l'Everest (du 5/22/2019 6:00 au 12/14/2019 13:00)")
# plt.plot(Ltime,WDIR,'.',markersize='2',label='direction du vent (capteur1)')
# plt.plot(Ltime,WDIR_2,'.',markersize='2',label='(capteur2, correction de 180°)')
# plt.grid()
# plt.legend()

plt.subplot(nb,1,1)
plt.plot(Ltime[moitpas:-moitpas],WDIRliss,label='capteur 1')
plt.plot(Ltime[moitpas:-moitpas],WDIR_2liss,label='capteur 2')
plt.ylabel('direction lissée (sur '+str(pas)+'h) du vent [°]')
plt.grid()
plt.legend()
plt.xticks(color='w')

# plt.subplot(nb,1,3)
# plt.plot(Ltime,WS_MAX, label='vitesse maximum pendant 5s du vent sur 1h ')
# plt.plot(Ltime,WS_MAX_2)
# plt.grid()
# plt.legend()

plt.subplot(nb,1,2)
plt.plot(Ltime[moitpas:-moitpas],WS_AVGliss,label='moyenne (sur 1h) lissée (sur '+str(pas)+'h)')
# plt.plot(Ltime,WS_AVG,label='vitesse moyenne du vent sur 1h')
plt.ylabel('vitesse moyenne du vent [m/s]')
# plt.plot(Ltime,WS_AVG_2)
plt.grid()
plt.legend()

# plt.subplot(nb,1,5)
# plt.plot(Ltime[moitpas:-moitpas],[y/x if x>0 else np.nan for x,y in zip(WS_AVGliss,WS_MAXliss)],label='WS_MAX/WS_AVG')
# plt.grid()
# plt.legend()

plt.savefig("C:/Users/Toshiba/Documents/stage L3/Documents/figures/mesures_vent_SC.png",bbox_inches="tight")
plt.show()

#plot vitesse du vent en fonction de sa direction
plt.plot(WDIR,WS_AVG,'.',markersize='1')
plt.xlabel('direction du vent [°]')
plt.ylabel('vitesse moyenne du vent [°]')
plt.show()

#somme de l'érosion par intervalle de direction du vent
rhos0=300
rhos0ind=0
rhosmax=450
tDR=24*3600
rhos,AC,EReff=erosion3C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
ERefftot=list(np.sum(EReff,axis=1))

step=5
ran=range(50,285+1,step)
#ran=range(45,285+1,30)
J=[k for k in ran]
J2=[0 for _ in ran]
J3=[k+step/2 for k in ran]
compt=0
for k in range(82277+moitpas,82277+4950-moitpas):
    i=0
    test=True
    while test:
        if J[i]<=WDIRliss[compt]<J[i+1]:
            J2[i]+=-ERefftot[k]
            test=False
        i+=1
    compt+=1

#plot erosion en fonction de la direction du vent (par intervalle)
plt.figure(figsize=(10,5))
plt.plot(J3,[x if x!=0 else np.nan for x in J2],'.')
plt.xlabel('direction du vent [°]')
plt.ylabel('érosion effective [mmwe]')
plt.title('érosion éolienne par intervalle de '+str(step)+'° de direction du vent du 22/5/2019 au 14/12/2019 (capteur 1, pas de 25h)')
plt.savefig("C:/Users/Toshiba/Documents/stage L3/Documents/figures/erosion(direction_vent).png",bbox_inches="tight")
plt.show()

## double plot erosion et vitesse du vent en fonction de la direction du vent
fig, ax1 = plt.subplots()
fig.set_size_inches(7, 5)
color = 'tab:blue'
ax1.set_xlabel('direction du vent [°]')
ax1.set_ylabel('vitesse moyenne du vent [m.s-1]', color=color)
ax1.plot(WDIR,WS_AVG,'.',markersize='1',color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('érosion effective [mmwe]', color=color)
ax2.plot(J3,[x if x!=0 else np.nan for x in J2],'.',color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.grid()
fig.suptitle("direction du vent et de l'érosion au SCG du 22/5/2019 au 14/12/2019")
plt.savefig("C:/Users/Toshiba/Documents/stage L3/Documents/figures/VV_ER(DV)_SC.png",bbox_inches="tight")
plt.show()

## rose des vents por l'érosion
#(http://youarealegend.blogspot.com/2008/09/windrose.html)
from windrose import WindroseAxes
size=18
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)

wd=WDIR
ws=[-x for x in ERefftot[82277:82277+4950]]
def new_axes():
    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
    rect = [0.1, 0.1, 0.8, 0.8]
    ax = WindroseAxes(fig, rect)
    fig.add_axes(ax)
    return ax

#A stacked histogram with normed (displayed in percent) (normed=True/False) results
#attention, c'est juste un histogramme du vent ! l'intendité de l'érosion n'est en fait pas prise en compte...
ax = new_axes()
ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
#plt.title('Erosion cumulée modélisée [mm w.e] selon les directions mesurées du vent \n du 22/5/2019 au 14/12/2019 au Col Sud')
plt.text(0.4*np.pi,45,'Erosion [mmwe]',fontsize='20')
plt.savefig("C:/Users/Toshiba/Documents/stage L3/Documents/figures/rose_des_erosions.png",bbox_inches="tight")
plt.show()


