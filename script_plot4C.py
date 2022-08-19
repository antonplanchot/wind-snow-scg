import numpy as np
import matplotlib.pyplot as plt

##preparation plot pluriannuel
def plotplurian4C(choice,date,string,SF2,rhos,AC,EReff):
    SFcumul=[np.sum(SF2[:k+1]) for k in range(len(SF2))]
    ACtot=np.sum(AC,axis=1)
    ERefftot=list(np.sum(EReff,axis=1))

    n=4
    plt.figure(figsize=(12,12))
    print(string)
    plt.subplot(n,1,1)
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
    plt.plot(date,AC[:,3],color='b',label='CF')
    plt.plot(date,AC[:,2],color='maroon',label='CFI')
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
    plt.plot(date,rhos[:,2],color='maroon',label='masse volumique neige CFI')
    plt.plot(date,rhos[:,1],color='k',label='masse volumique neige CI')#marker=',',linestyle='None',
    plt.ylabel('ρ [kg/m3]')
    plt.xlabel('temps [années]')
    plt.legend()
    plt.xlim(2009.5,2020.5)
    plt.gca().xaxis.set_ticks(range(2010,2020,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\3C_2010_2020"+string+".png",bbox_inches="tight")
    plt.show()
    if choice!=0:
        plt.pause(1)
        plt.close()
    print(SFcumul[-1]-ACtot[-1])

##préparation plot sensibilité
def plotsensibilite4C(var,strvar,string,Lrhos_mean,LACfinal,Loccurence_ER,LEReff_moy):
    size=13
    parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
    plt.rcParams.update(parameters)
    
    plt.figure(figsize=(10,3.5))
    marker='*'
    
    plt.subplot(1,2,1)
    plt.plot(var,Lrhos_mean,marker,label='masse volumique CI \n (moyenne temporelle)')
    plt.grid()
    plt.ylabel('ρ [kg/m3]')
    plt.xlabel(strvar)
    plt.legend()
    
    plt.subplot(1,2,2)
    #plt.plot(var,LACtot_mean,'k'+marker,color='orange',label='AC totale moyenne (moyenne temporelle)')
    plt.plot(var,[x/10 for x in LACfinal],color='k', marker=marker,linestyle='None',label='acc. nette annuelle')
    plt.grid()
    plt.ylabel('accumulation [mm w.e. an-1]')
    plt.xlabel(strvar)
    plt.legend()
    
    # plt.subplot(2,2,3)
    # plt.plot(var,Loccurence_ER,'r'+marker,label="occurences d'épisodes d'érosion")
    # plt.grid()
    # plt.ylabel("occurences/an")
    # plt.xlabel(strvar)
    # plt.legend()
    #
    # plt.subplot(2,2,4)
    # plt.plot(var,LEReff_moy,'r'+marker,label="érosion effective moyenne par occurence")
    # plt.grid()
    # plt.ylabel('mm w.e./occurence')
    # plt.xlabel(strvar)
    # plt.legend()
    
    plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\3C_2010_2020"+string+".png")
    plt.show()


##preparation plot annuel
def plotannuel4C(choice,mois,string,U2,SF2,rhos,AC,EReff):
    print(type(SF2))
    SFcumul=[np.sum(SF2[:k+1]) for k in range(len(SF2))]
    ERefftot=list(np.sum(EReff,axis=1))

    n=4
    size=13
    parameters = {'axes.labelsize': size,'axes.titlesize': size,'legend.fontsize': size,'xtick.labelsize': size,'ytick.labelsize': size}
    plt.rcParams.update(parameters)
    plt.figure(figsize=(12,9))
    print(string)

    plt.subplot(n,1,1)
    #plt.title(string)
    plt.plot(mois,SFcumul[:365*24],color='k',label='accumulation brute')
    plt.plot(mois,AC[:,3][:365*24],color='b',label='CF')
    plt.plot(mois,AC[:,2][:365*24],color='maroon',label='CFI')
    plt.plot(mois,AC[:,1][:365*24],color='y',label='CI')
    plt.plot(mois,AC[:,0][:365*24],color='r',label='CD')
    plt.ylabel('accumulation nette \n [mmwe]')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(1,14,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,2)
    plt.plot(mois,SF2[:365*24],color='b',label='Précipitation')
    plt.plot(mois,ERefftot[:365*24],color='r',label='Erosion effective totale')
    plt.ylabel('[mmwe/h]')
    plt.legend()
    plt.gca().xaxis.set_ticks(range(1,14,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')

    plt.subplot(n,1,3)
    plt.plot(mois,rhos[:,2][:365*24],color='maroon',label='masse volumique neige CFI')
    plt.plot(mois,rhos[:,1][:365*24],color='k',label='masse volumique neige CI')#marker=',',linestyle='None',
    plt.ylabel('ρ [kg/m3]')
    plt.xlim(0.4,13.6)
    plt.legend()
    plt.gca().xaxis.set_ticks(range(1,14,1))
    plt.gca().xaxis.grid()
    plt.gca().yaxis.grid()
    plt.xticks(color='w')

    # plt.subplot(n,1,4)
    # plt.plot(mois,U2[:365*24],color='g',label='vitesse du vent U')#marker=',',linestyle='None',
    # plt.ylabel('U [m/s]')
    # plt.xlabel('temps [mois]')
    # plt.legend()
    # plt.gca().xaxis.set_ticks(range(1,13,1))
    # plt.gca().xaxis.grid()
    # plt.gca().yaxis.grid()

    plt.savefig(r"C:\Users\Toshiba\Documents\stage L3\Documents\figures\3C_2010"+string+".png")
    plt.show()
    if choice!=0:
        plt.pause(1)
        plt.close()