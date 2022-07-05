import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from script_erosion3C_1 import erosion3C

U=np.array([])
P=np.array([])
T=np.array([])
SF=np.array([])
Lyear=['19500101-19600101','19600101-19700101','19700101-19800101','19800101-19900101','19900101-20000101','20000101-20100101','20100101-20200101']
deb=0
print(Lyear[deb][:8],'-->','20200101')
for y in range(deb,7):
    data=nc.Dataset("forcings/forcings_"+Lyear[y]+".nc")
    U=np.concatenate((U,data.variables['U2'][:,0,0]))
    P=np.concatenate((P,data.variables['PRES'][:,0,0]))
    T=np.concatenate((T,data.variables['T2'][:,0,0]))
    SF=np.concatenate((SF,data.variables['RRR'][:,0,0]))
nb_date=len(U)


##création des variables de temps adequates

date=[]
year=int(Lyear[deb][:4])
an=year
while an<2020:
    a=0
    if an%4==0:
        nb_jour=366
    else:
        nb_jour=365
    while a<nb_jour*24:
        date.append(an+a/(24*nb_jour))
        a+=1
    an+=1
if year!=1950:
    date.remove(year)
date.append(2020)
#début 01/01/1950 00:00
#fin 01/01/2020 00:00


bisextile=False
Lmois=[31,28+bisextile,31,30,31,30,31,31,30,31,30,31]
mois=[]
moi=1
while moi<13:
    a=0
    nb_jour=Lmois[moi-1]
    while a<nb_jour*24:
        mois.append(moi+a/(24*nb_jour))
        a+=1
    moi+=1


## plot serie temporelle
def plotplurian3Clong(date,string,SF2,rhos,AC,EReff):
    SFcumul=[np.sum(SF2[:k+1]) for k in range(len(SF2))]
    ACtot=np.sum(AC,axis=1)
    ERefftot=list(np.sum(EReff,axis=1))

    n=4
    plt.figure(figsize=(20,12))
    print(string)
    plt.subplot(n,1,1)
    plt.plot(date,SFcumul,color='b',label='acc. brute totale')
    plt.plot(date,ACtot,color='k',label='acc. nette totale')
    plt.ylabel('accumulation totale \n [mmwe]')
    plt.legend()
    plt.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,2)
    plt.plot(date,AC[:,2],color='b',label='CF')
    plt.plot(date,AC[:,1],color='y',label='CI')
    plt.plot(date,AC[:,0],color='r',label='CD')
    plt.ylabel('accumulation nette \n [mmwe]')
    plt.legend()
    plt.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,3)
    plt.plot(date,SF2,color='b',label='Précipitation')
    plt.plot(date,ERefftot,color='r',label='Erosion effective totale')
    plt.ylabel('[mmwe/h]')
    plt.legend()
    plt.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,4)
    plt.plot(date,rhos[:,1],color='k',label='masse volumique neige CI')#marker=',',linestyle='None',
    plt.ylabel('ρ [kg/m3]')
    plt.xlabel('temps [années]')
    plt.legend()
    plt.grid()
    plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\3C_1950_2020"+string+".png")
    plt.show()


rhos0=300
rhos0ind=0
rhosmax=450
tDR=24*3600
rhos,AC,EReff=erosion3C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
SFcumul=[np.sum(SF[:k+1]) for k in range(len(SF))]
ERefftot=list(np.sum(EReff,axis=1))
string='_rhos0='+str(rhos0)+'_rhosmax='+str(rhosmax)+'_tDR='+str(tDR/3600)+'h_fSF=1_fU=1'
plotplurian3Clong(date,string,SF,rhos,AC,EReff)


## analyse différences d'accumulation entre années

#échelle annuelle
sum_an_ERefftot=[]
sum_an_SF=[]
sum_ja_SF=[]
moy_an_U=[]
moy_ja_U=[]
var_an_U=[]
var_ja_U=[]
lenyeartot=0
for k in range(year,2020):#on va jusqu'au 31/12/2019 23:00
    lenyear=(365+(k%4==0))*24
    sum_an_ERefftot.append(-np.sum(ERefftot[lenyeartot:lenyeartot+lenyear]))
    sum_an_SF.append(np.sum(SF[lenyeartot:lenyeartot+lenyear]))
    sum_ja_SF.append(np.sum(SF[lenyeartot+6*30*24:lenyeartot+lenyear-4*30*24]))
    moy_an_U.append(np.mean(U[lenyeartot:lenyeartot+lenyear]))
    moy_ja_U.append(np.mean(U[lenyeartot+6*30*24:lenyeartot+lenyear-4*30*24]))
    var_an_U.append(np.var(U[lenyeartot:lenyeartot+lenyear]))
    var_ja_U.append(np.var(U[lenyeartot+6*30*24:lenyeartot+lenyear-4*30*24]))
    lenyeartot+=lenyear

AC_an=[x-y for x,y in zip(sum_an_SF,sum_an_ERefftot)]

xtime=range(year,2020)

plt.plot(xtime,AC_an)
plt.show()

# plt.plot(sum_an_SF,AC_an,'*')
# plt.show()#tendance faible
# r=np.corrcoef(sum_an_SF,AC_an)
# print(r[0,1])


##accumulation en fonction des précipitations totales de juillet/aout
r=round(np.corrcoef(sum_ja_SF,AC_an)[0,1],3)
a,b=np.polyfit(sum_ja_SF,AC_an,1)
plt.plot(sum_ja_SF,AC_an,'*')
plt.plot([min(sum_ja_SF),max(sum_ja_SF)],[a*min(sum_ja_SF)+b,a*max(sum_ja_SF)+b])
plt.xlabel('somme de précipitation sur juillet-aout [m/s]')
plt.ylabel('accumulation annuelle [mmwe]')
plt.text(70,30,'r='+str(r))
plt.title('1950-2020')
plt.show()#tendance
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\sum_SF_ja.png",bbox_inches="tight")

# plt.plot(moy_an_U,AC_an,'*')
# plt.show()#no tendance

##accumulation en fonction de la moyenne du vent de juillet/aout
r=round(np.corrcoef(moy_ja_U,AC_an)[0,1],3)
a,b=np.polyfit(moy_ja_U,AC_an,1)
plt.plot(moy_ja_U,AC_an,'*')
plt.plot([min(moy_ja_U),max(moy_ja_U)],[a*min(moy_ja_U)+b,a*max(moy_ja_U)+b])
plt.xlabel('moyenne du vent sur juillet-aout [m/s]')
plt.ylabel('accumulation annuelle [mmwe]')
plt.text(2,50,'r='+str(r))
plt.title('1950-2020')
plt.show()#tendance
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\moy_vent_ja.png",bbox_inches="tight")


# plt.plot(var_an_U,AC_an,'*')
# plt.show()#no tendance

##accumulation en fonction de la variance du vent de juillet/aout
r=round(np.corrcoef(var_ja_U,AC_an)[0,1],3)
a,b=np.polyfit(var_ja_U,AC_an,1)
plt.plot(var_ja_U,AC_an,'*')
plt.plot([min(var_ja_U),max(var_ja_U)],[a*min(var_ja_U)+b,a*max(var_ja_U)+b])
plt.xlabel('variance du vent sur juillet-aout')
plt.ylabel('accumulation annuelle [mmwe]')
plt.text(2,50,'r='+str(r))
plt.title('1950-2020')
plt.show()#tendance
#plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\var_vent_ja.png",bbox_inches="tight")

# plt.plot(moy_an_U,sum_an_SF,'*')
# plt.show()#no tendance


##trouver période discriminante

#jour début i=58
LrsumSF=[]
LrmoyU=[]
LrvarU=[]
for i in range(150):
    a,b,c=[],[],[]
    lenyeartot=0
    for k in range(year,2020):#on va jusqu'au 31/12/2019 23:00
        lenyear=(365+(k%4==0))*24
        a.append(np.sum(SF[lenyeartot+(4*30+i)*24:lenyeartot+lenyear-4*30*24]))
        b.append(np.mean(U[lenyeartot+(4*30+i)*24:lenyeartot+lenyear-4*30*24]))
        c.append(np.var(U[lenyeartot+(4*30+i)*24:lenyeartot+lenyear-4*30*24]))
        lenyeartot+=lenyear
    LrsumSF.append(np.corrcoef(a,AC_an)[0,1])
    LrmoyU.append(np.corrcoef(b,AC_an)[0,1])
    LrvarU.append(np.corrcoef(c,AC_an)[0,1])


plt.plot(LrvarU,label='r(var(U),Acc. annuelle) : coef. de corrélation')
plt.plot(LrmoyU,label='r(moy(U),Acc. annuelle)')
plt.plot(LrsumSF,label='r(somme(SF),Acc. annuelle)')
plt.xlabel('jour de début de la période échantillonée [jour après le 1er mai] ')
plt.legend()
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\r_deb.png",bbox_inches="tight")
plt.show()


#jour fin i=90
LrsumSF=[]
LrmoyU=[]
LrvarU=[]
for i in range(150):
    a,b,c=[],[],[]
    lenyeartot=0
    for k in range(year,2020):#on va jusqu'au 31/12/2019 23:00
        lenyear=(365+(k%4==0))*24
        a.append(np.sum(SF[lenyeartot+(4*30+58)*24:lenyeartot+lenyear-(6*30-i)*24]))
        b.append(np.mean(U[lenyeartot+(4*30+58)*24:lenyeartot+lenyear-(6*30-i)*24]))
        c.append(np.var(U[lenyeartot+(4*30+58)*24:lenyeartot+lenyear-(6*30-i)*24]))
        lenyeartot+=lenyear
    LrsumSF.append(np.corrcoef(a,AC_an)[0,1])
    LrmoyU.append(np.corrcoef(b,AC_an)[0,1])
    LrvarU.append(np.corrcoef(c,AC_an)[0,1])

plt.plot(LrvarU,label='r(var(U),Acc. annuelle) : coef. de corrélation')
plt.plot(LrmoyU,label='r(moy(U),Acc. annuelle)')
plt.plot(LrsumSF,label='r(somme(SF),Acc. annuelle)')
plt.xlabel('jour de fin de la période échantillonée [jour après le 1er juillet] ')
plt.legend()
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\r_fin.png",bbox_inches="tight")
plt.show()

## sur 1 ans : comparaison de 2 années extremes: 1988 et 1992
size=13
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)


n=8
plt.figure(figsize=(12,12))
string='3C_Erlimite_1988/1992_rhos0='+str(rhos0)+', rhosmax='+str(rhosmax)+', tDR='+str(tDR/3600)
deb1=(36*365+9)*24
fin1=deb1+365*24
deb2=(31*365+8)*24
fin2=deb2+365*24
plt.subplot(4,2,1)
#plt.title(string)
plt.plot(mois,[x-SFcumul[deb1]for x in SFcumul][deb1:fin1],color='k',label='précipitation totale')
plt.plot(mois,AC[:,2][deb1:fin1],color='b',label='CF 1986')
plt.plot(mois,AC[:,1][deb1:fin1],color='y',label='CI')
plt.plot(mois,[x-AC[:,0][deb1] for x in AC[:,0]][deb1:fin1],color='r',label='CD')
plt.ylabel('accumulation nette \n [mmwe]')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,14,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.xticks(color='w')

plt.subplot(4,2,3)
plt.plot(mois,SF[deb1:fin1],color='b',label='Précipitation')
plt.plot(mois,ERefftot[deb1:fin1],color='r',label='Erosion effective totale')
plt.ylabel('[mmwe/h]')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,14,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.xticks(color='w')

plt.subplot(4,2,5)
plt.plot(mois,rhos[:,1][deb1:fin1],color='k',label='masse volumique neige CI')#marker=',',linestyle='None',
plt.ylabel('ρ [kg/m3]')
plt.xlim(0.4,13.6)
plt.legend()
plt.gca().xaxis.set_ticks(range(1,14,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.xticks(color='w')

plt.subplot(4,2,7)
plt.plot(mois,U[deb1:fin1],color='g',label='vitesse du vent U')#marker=',',linestyle='None',
plt.ylabel('U [m/s]')
plt.xlabel('temps (mois)')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,13,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()

#------------------
plt.subplot(4,2,2)
#plt.title(string)
plt.plot(mois,[x-SFcumul[deb2]for x in SFcumul][deb2:fin2],color='k',label='précipitation totale')
plt.plot(mois,AC[:,2][deb2:fin2],color='b',label='CF 1981')
plt.plot(mois,AC[:,1][deb2:fin2],color='y',label='CI')
plt.plot(mois,[x-AC[:,0][deb2] for x in AC[:,0]][deb2:fin2],color='r',label='CD')
plt.ylabel('accumulation nette \n [mmwe]')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,14,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.xticks(color='w')

plt.subplot(4,2,4)
plt.plot(mois,SF[deb2:fin2],color='b',label='Précipitation')
plt.plot(mois,ERefftot[deb2:fin2],color='r',label='Erosion effective totale')
plt.ylabel('[mmwe/h]')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,14,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.xticks(color='w')

plt.subplot(4,2,6)
plt.plot(mois,rhos[:,1][deb2:fin2],color='k',label='masse volumique neige CI')#marker=',',linestyle='None',
plt.ylabel('ρ [kg/m3]')
plt.xlim(0.4,13.6)
plt.legend()
plt.gca().xaxis.set_ticks(range(1,14,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.xticks(color='w')

plt.subplot(4,2,8)
plt.plot(mois,U[deb2:fin2],color='g',label='vitesse du vent U')#marker=',',linestyle='None',
plt.ylabel('U [m/s]')
plt.xlabel('temps (mois)')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,13,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()

#plt.savefig('/home/planchoa/Documents/figures/'+string+'.png',bbox_inches="tight")

plt.show()

## moyenne des vents en 1950 de jullet a aout
def month2day(month):
    bisextile=False
    L=[31,28+bisextile,31,30,31,30,31,31,30,31,30,31]
    print(np.sum(L[:month]))

month2day(6)
month2day(9)
print(np.mean(U[181*24:273*24]))

#plot vent 1950
plt.plot(mois,U[:365*24],color='g',label='vitesse du vent U')#marker=',',linestyle='None',
plt.ylabel('U [m/s]')
plt.xlabel('temps (mois)')
plt.legend()
plt.gca().xaxis.set_ticks(range(1,13,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.show()


