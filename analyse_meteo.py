import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates
from datetime import datetime

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


data=pd.read_csv("/home/planchoa/Documents/code_amory_5_juin/South_Col_20210630.csv")
data=np.array(data)


#82277 : '5/22/2019 6:00'

m=4950
#82277 + 4950 : '12/14/2019 13:00'


WS_MAX=data[:,4][:m]
WS_AVG=data[:,5][:m]
WDIR=data[:,6][:m]
WS_MAX_2=data[:,7][:m]
WS_AVG_2=data[:,8][:m]
WDIR_2=data[:,9][:m]
date=data[:,0]

for k in range(len(WDIR_2)):
    WDIR_2[k]+=180
    if WDIR_2[k]>=360:
        WDIR_2[k]-=360

Ltime=[datetime.strptime(x, '%m/%d/%Y %H:%M') for x in date][:m]
#Ltime = matplotlib.dates.date2num(Ltime)


WDIRliss,WDIR_2liss=[],[]
WS_MAXliss,WS_AVGliss=[],[]
pas=24*3+1 #nombre impair [h] (pas de lissage)
moitpas=pas//2
for k in range(moitpas,len(WDIR)-moitpas):
    WDIRliss.append(np.mean( WDIR[k-moitpas:k+moitpas+1]))
    WDIR_2liss.append(np.mean( WDIR_2[k-moitpas:k+moitpas+1]))
    WS_MAXliss.append(np.mean(  WS_MAX[k-moitpas:k+moitpas+1]))
    WS_AVGliss.append(np.mean(  WS_AVG[k-moitpas:k+moitpas+1]))

WDIR_2liss[-310:]=310*[np.nan]

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

plt.savefig('/home/planchoa/Documents/figures/mesures_vent_SC.png',bbox_inches="tight")
plt.show()





plt.plot(WDIR,WS_AVG,'.',markersize='1')
plt.xlabel('direction du vent [°]')
plt.ylabel('vitesse moyenne du vent [°]')
plt.show()



## double plot
fig, ax1 = plt.subplots()
#fig.set_size_inches(7, 7)
color = 'tab:blue'
ax1.set_xlabel('direction du vent [°]')
ax1.set_ylabel('vitesse moyenne du vent [m.s-1]', color=color)
ax1.plot(WDIR,WS_AVG,'.',markersize='1',color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

color = 'tab:red'
ax2.set_ylabel('érosion effective [mmwe]', color=color)
ax2.plot(J3,J2,'.',color=color)
ax2.tick_params(axis='y', labelcolor=color)

#plt.gca().xaxis.set_ticks(range(0,350,10))
# plt.gca().xaxis.grid()
# plt.gca().yaxis.grid()
fig.tight_layout()
plt.grid()
fig.suptitle("direction du vent et de l'érosion au SCG du 22/5/2019 au 14/12/2019")
plt.savefig('/home/planchoa/Documents/figures/VV_ER(DV)_SC.png',bbox_inches="tight")

plt.show()








## rose des vents
# import plotly.graph_objects as go
#
# fig = go.Figure()
#
# fig.add_trace(go.Barpolar(r=[0,0]+J2+[0,0],name='érosion',))
#
# fig.show()

##autre

from windrose import WindroseAxes
import matplotlib.cm as cm
from numpy.random import random
from numpy import arange

#Create wind speed and direction variables
# ws = J2
# wd = J3

size=18
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)

wd=WDIR
ws=[-x for x in ERefftot[82277:82277+4950]]
#A quick way to create new windrose axes...
def new_axes():
    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
    rect = [0.1, 0.1, 0.8, 0.8]
    ax = WindroseAxes(fig, rect)# axisbg='w')
    fig.add_axes(ax)
    return ax

#...and adjust the legend box
# def set_legend(ax):
#     l = ax.legend(axespad=-0.10)
#     plt.setp(l.get_texts(), fontsize=8)

#A stacked histogram with normed (displayed in percent) results :
ax = new_axes()
ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
# set_legend(ax)
#plt.title('Erosion cumulée modélisée [mm w.e] selon les directions mesurées du vent \n du 22/5/2019 au 14/12/2019 au Col Sud')
plt.text(0.4*np.pi,45,'Erosion [mmwe]',fontsize='20')
plt.savefig('/home/planchoa/Documents/figures/rose_des_erosions.png',bbox_inches="tight")
plt.show()





