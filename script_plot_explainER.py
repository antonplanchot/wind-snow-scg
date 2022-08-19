import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from script_erosion1C_1 import erosion1C
from script_erosion2C_1 import erosion2C
from script_erosion3C_2 import erosion3C
from script_erosion4C_1 import erosion4C

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
    SFcumul=[]
    a=0
    for x in SF2:
        a+=x
        SFcumul.append(a)
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

    plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\3C_1950_2020"+string+".png")
    plt.show()


## plot serie temporelle simplifiée
def plotplurian3Clongsimp(date,string,SF2,rhos,AC,EReff):
    SFcumul=[]
    a=0
    for x in SF2:
        a+=x
        SFcumul.append(a)
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

    plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\3C_1950_2020"+string+".png")
    plt.show()


## analyse différences d'accumulation entre années (3C)
rhos0=300
rhos0ind=0
rhosmax=450
tDR=24*3600
rhos,AC,EReff=erosion3C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
SFcumul=[]
a=0
for x in SF:
    a+=x
    SFcumul.append(a)
ERefftot=list(np.sum(EReff,axis=1))
#string='_rhos0='+str(rhos0)+'_rhosmax='+str(rhosmax)+'_tDR='+str(tDR/3600)+'h_fSF=1_fU=1'
#plotplurian3Clongsimp(date,string,SF,rhos,AC,EReff)


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
    var_an_U.append(np.std(U[lenyeartot:lenyeartot+lenyear]))
    var_ja_U.append(np.std(U[lenyeartot+6*30*24:lenyeartot+lenyear-4*30*24]))
    lenyeartot+=lenyear

AC_an=[x-y for x,y in zip(sum_an_SF,sum_an_ERefftot)]

xtime=range(year,2020)

plt.plot(xtime,AC_an)
plt.show()

plt.plot(sum_an_SF,AC_an,'*')
plt.show()#tendance faible
r=np.corrcoef(sum_an_SF,AC_an)
print(r[0,1])


##accumulation en fonction des précipitations totales de juillet/aout
r=round(np.corrcoef(sum_an_SF,AC_an)[0,1],3)
a,b=np.polyfit(sum_an_SF,AC_an,1)
plt.plot(sum_an_SF,AC_an,'*')
plt.plot([min(sum_an_SF),max(sum_an_SF)],[a*min(sum_an_SF)+b,a*max(sum_an_SF)+b])
plt.xlabel('précipitations annuelles [mm w.e.]')
plt.ylabel('accumulation nette annuelle [mm w.e.]')
plt.text(230,30,'r='+str(r))
#plt.title('1950-2020')
plt.show()#tendance
#plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\sum_SF_ja.png",bbox_inches="tight")

# plt.plot(moy_an_U,AC_an,'*')
# plt.show()#no tendance

##accumulation en fonction de la moyenne du vent de juillet/aout
r=round(np.corrcoef(moy_ja_U,AC_an)[0,1],3)
a,b=np.polyfit(moy_ja_U,AC_an,1)
plt.plot(moy_ja_U,AC_an,'*')
plt.plot([min(moy_ja_U),max(moy_ja_U)],[a*min(moy_ja_U)+b,a*max(moy_ja_U)+b])
plt.xlabel('moyenne du vent sur juillet-août [m s-1]')
plt.ylabel('accumulation nette annuelle [mm w.e.]')
plt.text(2,50,'r='+str(r))
#plt.title('1950-2020')
plt.show()#tendance
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\moy_vent_ja.png",bbox_inches="tight")


# plt.plot(var_an_U,AC_an,'*')
# plt.show()#no tendance

##accumulation en fonction de la variance du vent de juillet/aout
r=round(np.corrcoef(var_ja_U,AC_an)[0,1],3)
a,b=np.polyfit(var_ja_U,AC_an,1)
plt.plot(var_ja_U,AC_an,'*')
plt.plot([min(var_ja_U),max(var_ja_U)],[a*min(var_ja_U)+b,a*max(var_ja_U)+b])
plt.xlabel('écart-type du vent sur juillet-août')
plt.ylabel('accumulation nette annuelle [mm w.e.]')
plt.text(2.4,50,'r='+str(r))
#plt.title('1950-2020')
plt.show()#tendance
#plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\var_vent_ja.png",bbox_inches="tight")

# plt.plot(moy_an_U,sum_an_SF,'*')
# plt.show()#no tendance


r=round(np.corrcoef(moy_ja_U,var_ja_U)[0,1],3)
a,b=np.polyfit(moy_ja_U,var_ja_U,1)
plt.plot(moy_ja_U,var_ja_U,'*')
plt.plot([min(moy_ja_U),max(moy_ja_U)],[a*min(moy_ja_U)+b,a*max(moy_ja_U)+b])
plt.xlabel('moyenne de la vitesse du vent sur juillet-août [m s-1]')
plt.ylabel('écart-type vit. du vent sur juillet-août')
plt.text(3,2.4,'r='+str(r))
#plt.title('1950-2020')
plt.show()#tendance

#____________________________________________________
##calcul de l'érosion pour tout les modèles
rhos0=300
rhos0ind=0
rhosmax=450
tDR=24*3600
rhos1,AC1,EReff1=erosion1C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
rhos2,AC2,EReff2=erosion2C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
rhos,AC,EReff=erosion3C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
rhos4,AC4,EReff4=erosion4C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
SFcumul=[]
a=0
for x in SF:
    a+=x
    SFcumul.append(a)

#____________________________________________________
##série temporelle d'érosion sur 70 ans
Acc1=[]
Acc=[]
Acc4=[]
prec=[]
frac=[]
for k in range(0,70):
    Acc.append(AC[(k+1)*365*24,0]-AC[k*365*24,0])
    Acc4.append(AC4[(k+1)*365*24,0]-AC4[k*365*24,0])
    prec.append(SFcumul[(k+1)*365*24]-SFcumul[k*365*24])
    frac.append(100*Acc[-1]/prec[-1])
print(np.mean(Acc),np.mean(Acc4))
print(np.std(Acc))

c=0
for x in Acc:
    if x>93:
        print(x,1950+c)
    c+=1
Lyear=[k for k in range(1950,2020)]
m=np.mean(frac)
plt.plot(Lyear,Acc)
plt.show()
print(np.std(frac))

plt.figure(figsize=(13,4))
size=13
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)
plt.subplot(1,2,1)
plt.bar(Lyear,prec,.7,label='précipitation brute')
plt.bar(Lyear,Acc,.7,label='précipitation nette')
plt.xlabel('temps [années]')
plt.ylabel("prcipitation annuelle [en mm w.e.]")
plt.legend()

plt.subplot(1,2,2)  
plt.plot(Lyear,frac)
plt.axhline(y=m, color='grey', linestyle='--',label="efficacité moyenne d'accumulation: 25%")
plt.grid()
plt.xlabel('temps [années]')
plt.ylabel("efficacité de l'accumulation [en %]")
plt.legend()
plt.show()


string='_rhos0='+str(rhos0)+'_rhosmax='+str(rhosmax)+'_tDR='+str(tDR/3600)+'h_fSF=1_fU=1'
plotplurian3Clong(date,string,SF,rhos,AC,EReff)

#______________________________________________________________
##statistiques des maxima annuels d'épaisseur de couches (1C, 3C, 4C)
import numpy as np
#1 couche
M=[]
for i in range(10):
    max=0
    for j in range(365*24):
        if AC1[i*365*24+j]>max:
            max=AC1[i*365*24+j]
    M.append(max)
print(np.mean(M),np.std(M),np.max(M))
    

#3 couches
M=[]
for i in range(10):
    max=0
    for j in range(365*24):
        if AC[i*365*24+j,1]>max:
            max=AC[i*365*24+j,1]
    M.append(max)
print(np.mean(M),np.std(M),np.max(M))

#4 couches
M=[]
for i in range(10):
    max=0
    for j in range(365*24):
        if AC4[i*365*24+j,1]>max:
            max=AC4[i*365*24+j,1]
    M.append(max)
print(np.mean(M),np.std(M),np.max(M)) 

 
#________________________________________________________________

an=1950
k=0
Pm=[0]*365*24
Um=[0]*365*24
ACm1=[0]*365*24
ACm2=[0]*365*24
ACm=[0]*365*24
ACm4=[0]*365*24
a1,a2,a3,a4=0,0,0,0

ranep=range(0,150+1,10)
L=[k for k in ranep]
L1=[0 for k in ranep]
L2=[0 for k in ranep]
L3=[0 for k in ranep]
L4=[0 for k in ranep]

L1dd=[0 for k in ranep]
L2dd=[0 for k in ranep]
L3dd=[0 for k in ranep]
L4dd=[0 for k in ranep]
while an<2020:
    bisext=0
    i=0
    test=False
    while i<365:
        if an//4!=0 or i!=59 or test:
            h=0
            while h<24:
                if an>1950 and i==0 and h==0:
                    a2,a3,a4=AC2[k,0],AC[k,0],AC4[k,0]
                Pm[i*24+h]+=P[k]
                Um[i*24+h]+=U[k]
                ACm1[i*24+h]+=AC1[k]
                ACm2[i*24+h]+=AC2[k,1]+AC2[k,0]-a2
                ACm[i*24+h]+=AC[k,1]+AC[k,2]+AC[k,0]-a3
                ACm4[i*24+h]+=AC4[k,1]+AC4[k,2]+AC4[k,3]+AC4[k,0]-a4
                if k>0:
                    for j in range(len(ranep)-1):
                        #densif
                        if ranep[j]<AC1[k]<=ranep[j+1] and EReff1[k+1]:
                            L1[j]+=1
                        if ranep[j]<AC2[k,1]<=ranep[j+1] and EReff2[k+1]:
                            L2[j]+=1
                        if ranep[j]<AC[k,1]<=ranep[j+1] and EReff[k+1,1]:
                            L3[j]+=1
                        if ranep[j]<AC4[k,1]<=ranep[j+1] and EReff4[k+1,1]:
                            L4[j]+=1
                        if ranep[j]<AC4[k,2]<=ranep[j+1] and EReff4[k+1,2]:
                            L4[j]+=1
                        if ranep[j]<AC[k,2]<=ranep[j+1] and EReff[k+1,2]:
                            L3[j]+=1
                        if ranep[j]<AC4[k,3]<=ranep[j+1] and EReff4[k+1,3]:
                            L4[j]+=1
                            
                        # #dédensif
                        if ranep[j]<AC1[k]<=ranep[j+1] and (not EReff1[k] and AC1[k]-AC1[k-1]!=0):
                            L1dd[j]+=1
                        if ranep[j]<AC2[k,1]<=ranep[j+1] and (not EReff2[k] and AC2[k,1]-AC2[k-1,1]!=0):
                            L2dd[j]+=1
                        if ranep[j]<AC[k,1]<=ranep[j+1] and (not EReff[k,1] and AC[k,1]-AC[k-1,1]!=0):
                            L3dd[j]+=1
                        if ranep[j]<AC4[k,1]<=ranep[j+1] and (not EReff4[k,1] and AC4[k,1]-AC4[k-1,1]!=0):
                            L4dd[j]+=1
                        if ranep[j]<AC4[k,2]<=ranep[j+1] and (not EReff4[k,2] and AC4[k+1,2]-AC4[k-1,2]!=0):
                            L4dd[j]+=1
                            
                        ##densif union dédensif
                        # if ranep[j]<AC1[k]<=ranep[j+1] and (EReff1[k] or AC1[k]-AC1[k-1]!=0):
                        #     L1[j]+=1
                        # if ranep[j]<AC2[k,1]<=ranep[j+1] and (EReff2[k] or AC2[k,1]-AC2[k-1,1]!=0):
                        #     L2[j]+=1
                        # if ranep[j]<AC[k,1]<=ranep[j+1] and (EReff[k,1] or AC[k,1]-AC[k-1,1]!=0):
                        #     L3[j]+=1
                        # if ranep[j]<AC4[k,1]<=ranep[j+1] and (EReff4[k,1] or AC4[k,1]-AC4[k-1,1]!=0):
                        #     L4[j]+=1
                        # if ranep[j]<AC4[k,2]<=ranep[j+1] and (EReff4[k,2] or AC4[k,2]-AC4[k-1,2]!=0):
                        #     L4[j]+=1
                k+=1
                h+=1
            i+=1
        else:
            k+=24
            test=True
            bisext+=1
    # plt.subplot(1,2,1)
    # plt.plot(mois,AC1[(an-1950)*(365+bisext)*24:(an+1-1950)*(365+bisext)*24], color='grey')
    # plt.subplot(1,2,2)
    # plt.plot(mois,AC4[(an-1950)*(365+bisext)*24:(an+1-1950)*(365+bisext)*24,1], color='red')
    an+=1
Pm=[k/70 for k in Pm]
Um=[k/70 for k in Um]
ACm1=[k/70 for k in ACm1]
ACm2=[k/70 for k in ACm2]
ACm=[k/70 for k in ACm]
ACm4=[k/70 for k in ACm4]

div=70
L1=[k/div for k in L1]
L2=[k/div for k in L2]
L3=[k/div for k in L3]
L4=[k/div for k in L4]

L1dd=[k/div for k in L1dd]
L2dd=[k/div for k in L2dd]
L3dd=[k/div for k in L3dd]
L4dd=[k/div for k in L4dd]

#___________________________________________________
#plot profondeur de densification (pour 2c, 3C et 4C)
size=11
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)
larg=2
stop=-1
f, (ax, ax2) = plt.subplots(2, 1,figsize=(6,4), sharex=True)

#plt.bar([k-3 for k in ranep],L1,larg,label="1 couche")
ax.bar([5+k-2 for k in ranep][:stop],L2[:stop],larg,label="2 couches")
ax.bar([5+k for k in ranep][:stop],L3[:stop],larg,label="3 couches")
ax.bar([5+k+2 for k in ranep][:stop],L4[:stop],larg,color='darkred',label="4 couches")
ax2.bar([5+k-2 for k in ranep][:stop],L2[:stop],larg)
ax2.bar([5+k for k in ranep][:stop],L3[:stop],larg)
ax2.bar([5+k+2 for k in ranep][:stop],L4[:stop],larg,color='darkred')
#plt.bar([5+k+4 for k in ranep],L4bis,larg,label="4 couches")
plt.gca().xaxis.set_ticks(range(0,150+1,10))
plt.gca().yaxis.set_ticks(range(0,10+1,2))
plt.gca().yaxis.grid()
ax2.set_xlabel('épaisseur de la couche [mm w.e.]')
ax.set_ylabel("nombre de densification [par an]", loc='top')
ax.legend(title = "Modèles à:")
#plt.ylim(0,15)
ax.text(-17,108,'(b)',fontsize=16)

ax.set_ylim(50, 110)  # outliers only
ax2.set_ylim(0, 12) 
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
plt.show()



#________________________________________________
#plot profondeur de dédensification (pour 3C et 4C)
size=12
parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
plt.rcParams.update(parameters)
larg=2
stop=-6
f, (ax, ax2) = plt.subplots(2, 1,figsize=(6,5), sharex=True)

#plt.bar([k-3 for k in ranep],L1,larg,label="1 couche")
#plt.bar([5+k-2 for k in ranep][:stop],L2[:stop],larg,label="2 couches")
ax.bar([5+k for k in ranep][:stop],L3dd[:stop],larg,color='darkorange',label="3 couches")
ax.bar([5+k+2 for k in ranep][:stop],L4dd[:stop],larg,color='darkred',label="4 couches")
ax2.bar([5+k for k in ranep][:stop],L3dd[:stop],larg,color='darkorange',label="3 couches")
ax2.bar([5+k+2 for k in ranep][:stop],L4dd[:stop],larg,color='darkred',label="4 couches")
#plt.bar([5+k+4 for k in ranep][:stop],L4bis[:stop],larg,label="4 couches")
plt.gca().xaxis.set_ticks(range(0,150+1,10))
plt.gca().yaxis.set_ticks(range(0,6+1,1))
plt.gca().yaxis.grid()
ax2.set_xlabel('épaisseur de la couche [mm w.e.]')
ax.set_ylabel("nombre de dédensification [par an]", loc='top')
ax.set_ylabel("nombre de dédensification [par an]", loc='top')
ax.legend(title = "Modèles à:")
#plt.ylim(0,15)

ax.set_ylim(7, 35)  # outliers only
ax2.set_ylim(0, 6) 
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
plt.show()


#_____________________________________________________
##accumulation annuelle de neige pour plusieurs modèles
#print(np.mean( [x*(y+5) for x,y in zip(L1,ranep)]),np.mean(L2),np.mean(L3),np.mean(L4))
#ACm_jour=[np.mean(ACm[k*24:(k+1)*24]) for k in range(365)]
plt.figure(figsize=([4,4]))
plt.plot(mois[1:],ACm1[1:],color='green',label='1 couche')
plt.plot(mois[1:],ACm2[1:],color='#1f77b4',label='2 couches')
plt.plot(mois[1:],ACm[1:],color='#ff7f0e',label='3 couches')
plt.plot(mois[1:],ACm4[1:],color='darkred',label='4 couches')
plt.gca().xaxis.set_ticks(range(1,13,1))
plt.gca().xaxis.grid()
plt.gca().yaxis.grid()
plt.xlabel('temps [mois]')
plt.ylabel('épaisseur totale de neige [mm w.e.]')
plt.text(-2,80,'(a)',fontsize=16)
plt.legend()
plt.show()


#____________________________________________________
##détermination période disciminante pour variabilité interannuelle 
#jour début i=58
LrsumSF=[]
LrmoyU=[]
LrvarU=[]
for i in range(105):
    a,b,c=[],[],[]
    lenyeartot=0
    for k in range(year,2020):#on va jusqu'au 31/12/2019 23:00
        lenyear=(365+(k%4==0))*24
        a.append(np.sum(SF[lenyeartot+(4*30+i)*24:lenyeartot+lenyear-(4*30+2)*24]))
        b.append(np.mean(U[lenyeartot+(4*30+i)*24:lenyeartot+lenyear-(4*30+2)*24]))
        c.append(np.var(U[lenyeartot+(4*30+i)*24:lenyeartot+lenyear-(4*30+2)*24]))
        lenyeartot+=lenyear
    LrsumSF.append(np.corrcoef(a,AC_an)[0,1]**2)
    LrmoyU.append(np.corrcoef(b,AC_an)[0,1]**2)
    LrvarU.append(np.corrcoef(c,AC_an)[0,1]**2)


plt.plot(LrvarU,label='r(std(U), AAN)^2')
plt.plot(LrmoyU,label='r(moy(U), AAN)^2')
plt.plot(LrsumSF,label='r(somme(P), AAN)^2')
plt.ylabel("coef. de corr. linéaire au carré")
plt.xlabel("date de début de l'échantillonnage [jour après le 1er mai] \n (fin le 31 août)")
plt.legend()
plt.grid()
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\r_deb.png",bbox_inches="tight")
plt.show()


#jour fin i=90
LrsumSF=[]
LrmoyU=[]
LrvarU=[]
for i in range(105):
    a,b,c=[],[],[]
    lenyeartot=0
    for k in range(year,2020):#on va jusqu'au 31/12/2019 23:00
        lenyear=(365+(k%4==0))*24
        a.append(np.sum(SF[lenyeartot+(4*30+58)*24:lenyeartot+lenyear-(6*30+2-i)*24]))
        b.append(np.mean(U[lenyeartot+(4*30+58)*24:lenyeartot+lenyear-(6*30-i)*24]))
        c.append(np.var(U[lenyeartot+(4*30+58)*24:lenyeartot+lenyear-(6*30-i)*24]))
        lenyeartot+=lenyear
    LrsumSF.append(np.corrcoef(a,AC_an)[0,1]**2)
    LrmoyU.append(np.corrcoef(b,AC_an)[0,1]**2)
    LrvarU.append(np.corrcoef(c,AC_an)[0,1]**2)

plt.plot(LrvarU,label='r(std(U), AAN)^2')
plt.plot(LrmoyU,label='r(moy(U), AAN)^2')
plt.plot(LrsumSF,label='r(somme(P), AAN)^2')
plt.xlabel("date de fin de l'échantillonnage [jour après le 1er juillet] \n (début le 28 juin) ")
plt.ylabel("coef. de corr. linéaire au carré")
plt.legend()
plt.grid()
plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\r_fin.png",bbox_inches="tight")
plt.show()


#______________________________________________________________
## sur 1 ans : comparaison de 2 années extremes: 1988 et 1992 (3C)
rhos0=300
rhos0ind=0
rhosmax=450
tDR=24*3600
rhos,AC,EReff=erosion3C(T,P,U,SF,rhos0,rhos0ind,rhosmax,tDR)
SFcumul=[]
a=0
for x in SF:
    a+=x
    SFcumul.append(a)
ERefftot=list(np.sum(EReff,axis=1))

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


#______________________________________________________
## moyenne des vents en 1950 de juillet à aout
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


