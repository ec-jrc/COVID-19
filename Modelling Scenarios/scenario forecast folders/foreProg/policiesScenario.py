# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 07:14:13 2020

@author: annunal
"""
import json
import numpy as np
import matplotlib.pyplot as plt
import os,random
from shutil import copyfile
from datetime import datetime
from datetime import timedelta
#from ECDC_model import ECDC_model
from analisi import analisi
from plotUtils import plotFig
from polyScenarioNat import polyNational
#def getCountries():
#    url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
#    if  os.path.exists('tmp.txt'):
#        os.remove('tmp.txt')
#    cmd='curl -o tmp.txt '+url
#    print (cmd)
#    os.system(cmd)
#    f=open('tmp.txt','r',encoding="utf-8")
#    rows=f.read().split('\n')
#    f.close()
#    return rows

#def getrows():
#    url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
#    if  os.path.exists('tmp.txt'):
#        os.remove('tmp.txt')
#    cmd='curl -o tmp.txt '+url
#    print (cmd)
#    os.system(cmd)
#    f=open('tmp.txt','r',encoding="utf-8")
#    rows=f.read().split('\n')
#    f.close()
#    return rows

def getlist(fname):
    f=open(fname,'r',encoding="utf-8")
    rows=f.read().split('\n')
    f.close()
    return rows

def getCountry(region, rows):
    
    for r in rows:
      if r =='':continue
      r=r.replace(',"','"')
      #print (r)
      p=r.split(',')
      if region==p[3]:
        cou=p[2]
    
    return cou

def getPopICU(dire0,country,reg=''):
    if reg=='':
        f=open(dire0+'Pop_ICUs.csv','r',encoding="utf-8")
        rows=f.read().split('\n')
        f.close()
        for row in rows:
            if row.split('\t')[0]==country:
                pop=float(row.split('\t')[1])
                icu=float(row.split('\t')[2])
                if pop==0: pop=1
                return pop,icu,icu/pop
    else:
        f=open(dire0+'Pop_ICUs_region.csv','r',encoding="utf-8")
        rows=f.read().split('\n')
        f.close()
        for row in rows:
            if row.split('\t')[0]==country and row.split('\t')[1]==reg:
                
                icu=float(row.split('\t')[2])
                return 0,icu,0
        return 0,0,getPopICU(dire0,cou)[2]        
    print ('country/reg '+country+ ''+reg+' not found in '+dire0+'Pop_ICUs.csv')
    return 0,0,0


    
def mainScenario(dire0,direIn,DAYFORE, c0=100,c1=100,FILTERCOUNTRY={},prefix='',calibrationPeriod=30):


    #analisi(0,0,dire0,DAYFORE)

    if not os.path.isdir(dire0+DAYFORE):
        os.makedirs(dire0+DAYFORE)

    #STUDY='BYREGION'
    STUDY='BYCOUNTRY'
    random.seed()
    #rows=getrows()
    FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES=dire0+DAYFORE+'/jrc-covid-19-all-days-by-country.csv'


    if not os.path.exists(FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES):
        #os.remove(FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES)
        url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
        cmd='curl -o "'+ FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES+ '" '+url
        print (cmd)
        os.system(cmd)

    FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS=dire0+DAYFORE+'/jrc-covid-19-all-days-by-regions.csv'

    if not os.path.exists(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS):
        #os.remove(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS)
        url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
        cmd='curl -o "'+ FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS+ '" '+url
        print (cmd)
        os.system(cmd)

    rows=getlist(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS)

    cous=[]
    EUData=[]
    REGSData=[]

    #countries=getCountries()
    countries=getlist(FILE_EPIDEMIOLOGY_WEBCRITECH_COUNTRIES)

    for r in countries:
        if r.startswith('Date') or r=='': continue
        cc=r.split(',')[2]
        if not cc in cous:
            cous.append(cc)

    # read countries data
    for c in cous:
        TIME=[]
        CUMUL=[]
        ICUs=[]
        NEWPOS=[]

        cum0=0
        for r in countries:
            if r.startswith('Date') or r=='': continue
            cc=r.split(',')[2]
            if cc==c:
                TIME.append(datetime.strptime(r.split(',')[0],'%Y-%m-%d'))
                if r.split(',')[5]=='':
                    cum=0
                else:
                    cum=int(r.split(',')[5])
                CUMUL.append(cum)
                icus=r.split(',')[10]
                if icus=='':
                    icus=0
                else:
                    icus=int(icus)
                ICUs.append(icus)
                NEWPOS.append(cum-cum0)
                cum0=cum
        EUData.append([c,TIME,CUMUL,NEWPOS,ICUs])


    # read regions data
    regions={}
    mcMax=0
    iculockMax=1e9
    iculockMin=0
    tunlock=-1e6
    tlock=-1e6
    waitTime=7

    for r in rows:
        if r.startswith('Date') or r=='': continue
        tv=datetime.strptime(r.split(',')[0],'%Y-%m-%d')
        cc=r.split(',')[2]
        reg=r.split(',')[3]
        if not reg in regions:
            regions[reg]=[cc,[],[],[]]
            cum0=0
        else:
            cum0=regions[reg][2][-1]
        cums=r.split(',')[6]
        if cums=='':
            cum=0
        else:
            cum=float(cums)
        #if reg=='Lombardia':
        #    print(tv,cum0,cum,cum-cum0)
        regions[reg][1].append(tv)
        regions[reg][2].append(cum)
        regions[reg][3].append(cum-cum0)
        

    for case in range(c0,c1+1):

        if case==0:
            lockMax=1e6  # means never
            lockMin=1e6  #1000
            StartTimeControl=0
            R0Target=0.7
            mcMax=0
        elif case==1:
            lockMax=20  
            lockMin=20  #1000
            StartTimeControl=30
            R0Target=0.7
            mcMax=0
        elif case==2:
            lockMax=100  # means never
            lockMin=100  #1000
            StartTimeControl=30
            R0Target=0.7
            mcMax=0
        elif case==3:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=0.7
            mcMax=0
        elif case==4:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=0.95
            mcMax=0
        elif case==5:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=1.2
            mcMax=0
        elif case==6:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=1.5
            mcMax=0
        elif case==7:
            lockMax=400  # means never
            lockMin=0  #1000
            StartTimeControl=30
            R0Target=0.95
            mcMax=0
    
        elif case==100:
            lockMax=1e6  # means never
            lockMin=1e6  #1000
            StartTimeControl=0
            R0Target=0.95
            mcMax=0
        elif case==101:
            lockMax=20  
            lockMin=20  #1000
            StartTimeControl=30
            R0Target=0.95
            mcMax=0
        elif case==102:
            lockMax=100  # means never
            lockMin=100  #1000
            StartTimeControl=30
            R0Target=0.95
            mcMax=0
        elif case==103:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=0.95
            mcMax=0

        elif case==104:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=0.7
            mcMax=0    
        elif case==105:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=1.2
            mcMax=0    
        elif case==106:
            lockMax=400  # means never
            lockMin=400  #1000
            StartTimeControl=30
            R0Target=1.5
            mcMax=0    
        elif case==107:
            lockMax=400  # means never
            lockMin=0  #1000
            StartTimeControl=30
            R0Target=0.95
            mcMax=0    
        elif case==110:
            lockMax=-1
            lockMin=-1
            StartTimeControl=30
            R0Target=0.7



    #---!-----------STUDY LOOP---------
        

        for STUDY in ['BYCOUNTRYNAT','BYCOUNTRYNAT_ANYREG','BYCOUNTRY','BYREGION']:
        #for STUDY in ['BYREGION']:
            try:
                if c0 !=-1:
                    dire1=dire0+'/'+DAYFORE+'/'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+"_"+STUDY
                else:
                    dire1=dire0+'/'+DAYFORE+'/imposed_'+prefix+'_'+STUDY
            
                if not os.path.exists(dire1):
                    os.mkdir(dire1)
                else:
                    if os.path.exists(dire1+'/000/Italy_newInfe.jpg'):
                        continue
                flog=open(dire1+"/'logCalc.txt","w")

                print (" * * " +STUDY +" * * ")
        #-------!----  MC LOOP ------------
                for mc in range(0,mcMax+1):
                    print ("mc="+format(mc))
                    if mc>0:
                        fmc1=random.randrange(-50,50)/100 #(90,110)/100.0
                        fmc2=random.randrange(-200,200)/100
                    else:
                        fmc1=0.0
                        fmc2=0.0
            
                    dire=dire1+'/'+"{:03d}".format(mc)+"/"
                    if not os.path.exists(dire):
                        os.mkdir(dire)
            
                    lockImplementation=20  # implemenattion of measures
                
                    print (direIn)
                    if STUDY=="BYCOUNTRY":
                        #fname=dire0+'/'+DAYFORE+'/FITTING_ALLCOUNTRIES_SIR_(-30 -1).json'
                        fname=direIn+ '/FITTING_ALLCOUNTRIES/data/FITTING_ALLCOUNTRIES_SIR_(-'+format(calibrationPeriod)+' -1).json'
                        #fname=dire0+'/'+DAYFORE+'/FITTING_ALLCOUNTRIES_SIR_(-30 -1)_13092020.json'
                        if not os.path.exists(fname):
                            fname=direIn+ '/FITTING_SELECTEDCOUNTRIES/data/FITTING_SELECTEDCOUNTRIES_SIR_(-'+format(calibrationPeriod)+' -1).json'
                    else:
                        #fname=dire0+'/'+DAYFORE+'/FITTING_ALLREGIONS_SIR_(-30 -1).json'
                        fname=direIn+ '/FITTING_ALLREGIONS/data/FITTING_ALLREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).json'
                        #fname=dire0+'/'+DAYFORE+'/FITTING_ALLREGIONS_SIR_(-30 -1)_13092020.json'
                        if not os.path.exists(fname):
                            fname=direIn+ '/FITTING_SELECTEDREGIONS/data/FITTING_SELECTEDREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).json'
                    #f=open(fname)
                    data=json.load(open(fname))

                    ecList = "Austria, Belgium, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg, Malta, Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia, Spain, Sweden".split(', ')
                    #ecList = "Greece".split(', ')
                    eucmParticip = "Iceland, Montenegro, North Macedonia, Norway, Serbia, Turkey, United Kingdom, Switzerland".split(', ')

                    if STUDY=='BYCOUNTRYNAT' or STUDY=='BYCOUNTRYNAT_ANYREG':
                        polyNational(direIn,dire,data, STUDY,ecList,lockMax,lockMin,R0Target,StartTimeControl,fmc1,fmc2)
                        continue

                    EU={}
                    lock=False

                    print(data.items)
                    for (reg,v) in data.items():

                        if STUDY=='BYREGION':
                            #cou=getCountry(reg,rows)
                            try:
                                cou=regions[reg][0]
                            except:
                                continue
                        else:
                            cou=reg
                
                        if cou=='Hungary':
                           cou=cou
                        else:
                            continue
                
                        print (cou,reg)
                        if not (cou in ecList or cou in eucmParticip):
                            continue
                        if "Russian Fed. Jewish" in reg:
                            reg=reg
                    
                        #if reg=='Lombardia' or reg=='Italy':
                        #    reg=reg
                        #else:
                        #    continue
                    
                    
                        r0=data[reg]['r0']+fmc1
                        if r0<0.1: r0=0.1
                        Trecov=data[reg]['Trecov']+fmc2
                        flog.write(format(mc)+','+dire+','+format(fmc1)+','+format(fmc2)+','+format(r0)+','+format(Trecov)+'\n')
                        N=data[reg]['population']
                        R=data[reg]['_ivp']['R']
                        I=data[reg]['_ivp']['I']
                        time0=data[reg]['time'][0]

                        if len(FILTERCOUNTRY)>0:
                            if cou in FILTERCOUNTRY:
                                lockdownDate=datetime.strptime(FILTERCOUNTRY[cou][0],'%Y-%m-%d')
                                date0=datetime.strptime(time0,'%Y/%m/%d')
                                StartTimeControl=(lockdownDate-date0).days
                            
                                R0Target=FILTERCOUNTRY[cou][1]
                                lockMax=FILTERCOUNTRY[cou][2]
                                lockMin=FILTERCOUNTRY[cou][3]
                                lockICUPerc=FILTERCOUNTRY[cou][4]
                                unlockICUPerc=FILTERCOUNTRY[cou][5]
                                waitTime=FILTERCOUNTRY[cou][6]
                            else:
                                continue
                    
                            print (cou,reg,I,Trecov,time0)
                        if STUDY=='BYCOUNTRY':
                            if lockMax==-1 and lockMin==-1:
                                'get ICU max by country'
                                p,icu,den=getPopICU(dire0,cou)
                                iculockMax=icu*lockICUPerc/100
                                iculockMin=icu*unlockICUPerc/100
                                lockMax=1e9
                                lockMin=-1
                        else:
                            if lockMax==-1 and lockMin==-1:
                                'get ICU max by region'
                                p,icu,den=getPopICU(dire0,cou,reg)
                                if p==0 and icu==0:
                                    icu=den*N
                                iculockMax=icu*lockICUPerc/100
                                iculockMin=icu*unlockICUPerc/100
                                lockMax=1e9
                                lockMin=-1
                        print(cou,reg)
                        I0=I
                        daysLock=0
                        if waitTime=='':
                            waitTime=14
                        S=N-I-R
                        II=[]
                        RR=[]
                        SS=[]
                        CUMPOS=[]
                        NEWPOS=[]
                        LOCK=[]
                        REPR=[]
                        CUM_14days=[]
                        II.append(I)
                        RR.append(R)
                        SS.append(S)
                        REPR.append(r0)
                        CUMPOS.append(N-S)
                        #NEWPOS.append(0)   # it is only to have the same number of data
                        CP0=N-S
                        LOCK.append(0)


                        #print(cou,reg,r0,Trecov,N)
    
                        dt=0.1
                        r00=r0
                        R0Lock=r00
                        R0UnLock=r00
                        firstLock=True
                        tlock=0
                        tunlock=-1e6
                        lock=False
                        unlocking=False
                        unlocked=False
    
                        xdata=np.linspace(0,30*6*int(1/dt),30*6*int(1/dt))*dt  # 6 months forecast
                        first = True
                        #f=open(dire+cou+'.csv','w')
                        for i in range(len(xdata)):
                           # if i==1191:
                           #     i=i
                    
                            dSdt=-r0/Trecov*S*I/N
                            dIdt=r0/Trecov*S*I/N - 1/Trecov*I
                            dRdt=1/Trecov*I
        
                            S +=dSdt*dt
                            I +=dIdt*dt
                            R +=dRdt*dt
                            CP =N-S
        
                            II.append(I)
                            RR.append(R)
                            SS.append(S)
                            CUMPOS.append(CP)
                            npo=(CP-CP0)/dt
                            NEWPOS.append(npo)
                            CP0=CP
                            #print(format(xdata[i])+','+format(S)+','+format(I)+','+format(R)+','+format(CP)+','+format(r0)+','+format(npo))
        
                            #  I     =  CURRENT POSITIVE !!!
                            #  -dSdt =  NEW POSITIVE
                            #  N-S   =  CUMULATIVE POSITIVE
                            if first:
                                cumIncidence_14days=(CP-CUMPOS[0])/0.1*14*100000.0/N
                                first=False
                            else:
                                if len(CUMPOS)>15/dt:
                                    cumIncidence_14days=(CP-CUMPOS[len(CUMPOS)-int(14/dt)])*100000.0/N           
                            #if cou=='France':
                            print(i,cumIncidence_14days,lockMax,lockMin)
                            CUM_14days.append(cumIncidence_14days)
                            ICUEstimate=npo*0.09
                            if xdata[i]>StartTimeControl:
                                if (cumIncidence_14days>lockMax or (ICUEstimate>iculockMax and xdata[i]-tunlock>waitTime)) and not lock :  # and (I0<lockMax):
                                
                                    tlock=xdata[i]
                                    lock=True
                                    R0Lock=r0
                                if (cumIncidence_14days<lockMin or (ICUEstimate<iculockMin and xdata[i]-tlock>waitTime)) and lock :
                                    tunlock=xdata[i]
                                    lock=False
                                    R0UnLock=r0
                            def transF(t,tlock,delta,R_0_start,R_0_end):
                                if t<tlock:
                                    return R_0_start
                                elif t>tlock+delta*2:
                                    return R_0_end
                                else:
                                    ff=10/delta
                                    return (R_0_start-R_0_end)/(1+np.exp(-ff*(-t+tlock+delta/2)))+R_0_end
        
                            if lock: # and xdata[i]>daysAfterLock+tlock:
                                #r0=min(0.7,r00)
                                if firstLock:
                                #    deltaDelay=14   # 1 week delay due to disease propagation
                                    firstLock=False
                                #else:
                                #    deltaDelay=0
                                deltaDelay=7
                                r0=transF(xdata[i],tlock+deltaDelay,lockImplementation,R0Lock,R0Target)
                            else:
                                deltaDelay=2   # immediate releasing of people

                                r0=transF(xdata[i],tunlock+deltaDelay,lockImplementation,R0UnLock,r00)
                            #print (xdata[i],r0,tlock,tunlock)
                            REPR.append(r0)
                            if lock:
                                LOCK.append(1)
                            else:
                                LOCK.append(0)
        
                            #print (i,r0, R0Lock,R0UnLock)

                        #f.close()
                        if not cou in EU:
                            EU[cou]=[]

                        EU[cou].append([reg,II,SS,CUMPOS,N,LOCK,NEWPOS,REPR,CUM_14days])
                
                        if False:
                            fig, ax = plt.subplots()
                            ax.set_title('Cumulative Infected '+ cou,fontsize=10,fontweight="bold")    
                            #plt.plot(myoptim.model.time.values, myoptim.model.predictions["I"].values, 'o')
                            #myoptim.model.predictions["I"].plot(ax=ax)
        
                         #   print (np.shape(xdata),np.shape(II))
        
                            plt.plot(xdata,II)

                            plt.show()
                    
    
                    rr=open(dire+'cumPos_byRegion.csv','w', encoding="utf-8")
                    cc=open(dire+'cumPos_byCountry.csv','w', encoding="utf-8")
                    Nrr=open(dire+'newPos_byRegion.csv','w', encoding="utf-8")
                    Ncc=open(dire+'newPos_byCountry.csv','w', encoding="utf-8")

                    lt=open(dire+'lockTime2.csv','w', encoding="utf-8")
                    lp=open(dire+'lockPop2.csv','w', encoding="utf-8")
                    ls=open(dire+'lockStatus2.csv','w', encoding="utf-8")
                    rt=open(dire+'rt2.csv','w', encoding="utf-8")

                    rr.write('days,Date')
                    cc.write('days,Date')
                    lp.write('days,Date')
                    ls.write('days,Date')
                    Nrr.write('days,Date')
                    Ncc.write('days,Date')
                    rt.write('days,Date')

                    #lt.write('Date')
                    for c in EU:
                        cc.write(','+c)
                        Ncc.write(','+c)
                        lp.write(','+c)
                        #lt.write(',Date,'+c)
                        for regs in EU[c]:
                            if STUDY=='BYREGION':
                                rr.write(','+regs[0]+'@'+c)
                                ls.write(','+regs[0]+'@'+c)
                                Nrr.write(','+regs[0]+'@'+c)
                            else:
                                rr.write(','+regs[0])
                                ls.write(','+regs[0])

                                Nrr.write(','+regs[0])



                    lt.write('Country,TotDays,TotPopLock,Population\n')
                    lp.write(',TOT\n')
                    cc.write(',TOT\n')
                    rr.write(',TOT\n')
                    Ncc.write(',TOT\n')
                    Nrr.write(',TOT\n')
                    rt.write('\n')
        
                    for i in range(len(xdata)):
                        ti=datetime.strptime(time0,'%Y/%m/%d')+timedelta(days=xdata[i])

                        cc.write(format(xdata[i])+','+format(ti))
                        rr.write(format(xdata[i])+','+format(ti))
                        lp.write(format(xdata[i])+','+format(ti))
                        ls.write(format(xdata[i])+','+format(ti))
                        rt.write(format(xdata[i])+','+format(ti))
                        if i>0:
                            Ncc.write(format(xdata[i])+','+format(ti))
                            Nrr.write(format(xdata[i])+','+format(ti))
    
                        tot=0
                        totCum=0
    
                        TotdaysLock=0
                        TotpopLock=0
                        for cou in EU:
                            totcou=0
                            totcouCum=0
                            newCum=0
                            Lock=0
                            popLock=0
                            for reg,II,SS,CUMPOS,N,LOCK,NEWPOS,REPR in EU[cou]:
                                # [reg,II,SS,CPCP,N,LOCK,NEWPOS,REPR])
                                #   0   1  2  3   4  5
                                #II=reg[1]
                                totcou +=II[i]
                                tot+=II[i]
                                #rr.write(','+format(int(reg[1][i])))
                                #SS=reg[2]
                                #CUMPOS=reg[3]
                                #NEWPOS=reg[6]
                                #LOCK=reg[5]
                                totcouCum +=CUMPOS[i]
                                totCum +=totcouCum
                                rr.write(','+format(int(CUMPOS[i])))
                                ls.write(','+format(LOCK[i]))
                                rt.write(','+format(REPR[i]))
                                if i>0:
                                    Nrr.write(','+format(int(NEWPOS[i])))  #format(int((CUMPOS[i]-CUMPOS[i-1])/(xdata[i]-xdata[i-1]))))
                                    newCum +=int(NEWPOS[i])
                                #N=reg[4]
                                #LOCK=reg[5]
            
                                if LOCK[i]==1:
                                    popLock +=N

                            cc.write(','+format(int(totcouCum)))  
                            if i>0:
                                Ncc.write(','+format(newCum))
                            totcouCum0=totcouCum
                            lp.write(','+format(popLock))
                            TotpopLock += popLock

                        rr.write(','+format(totCum)+'\n')
                        ls.write('\n')
                        cc.write(','+format(totCum)+'\n')
                        lp.write(','+format(TotpopLock)+'\n')
                        rt.write('\n')
                        if i>0:
                            Nrr.write(',\n')
                            Ncc.write(',\n')


                    for cou in EU:
                        TotdaysLock=0
                        TotpopLock=0
                
                        for i in range(len(xdata)):
                            atleastOnelocked=False
                            TotpopLock=0
                            totN=0
                            for reg in EU[cou]:
                                N=reg[4]
                                totN +=N
                                LOCK=reg[5]
                                if LOCK[i]==1:
                                    atleastOnelocked=True
                                    TotpopLock +=N
                            #if atleastOnelocked:
                                    TotdaysLock +=dt*N
        

                        lt.write(cou+','+format(int(TotdaysLock))+','+format(TotpopLock)+','+format(totN)+'\n')
                    cc.close()
                    rr.close()
                    lt.close()
                    lp.close()
                    ls.close()
                    Ncc.close()
                    Nrr.close()
                    rt.close()

                #dire=dire0+'/'+DAYFORE+'/'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+"_"+STUDY+'/'+"{:03d}".format(mc)
            
                TIME=[]
                for i in range(len(xdata)):
                  ti=datetime.strptime(time0,'%Y/%m/%d')+timedelta(days=xdata[i])
                  TIME.append(ti)
                fr=open(dire+'/times_icu_bycou.csv','w', encoding="utf-8")
                fr.write('Country, dateMin, DateMax, ICUMax\n')
                fr.close()

                for cou in EU:
                    NEWPOS=[]
                    ICUS=[]
                    ICUMin=[]
                    ICUMax=[]
                    MAXICU=[]
                    if ((cou=='France' or cou=='Denmark' or cou=='Greece' or cou=='Slovenia' or cou=='Iceland') and STUDY=='BYREGION'):
                      continue
                    if cou=='Iceland':
                        cou=cou
                    Pop,maxICU,den=getPopICU(dire0,cou)
                    for i in range(len(xdata)):
                        newCum=0
                        icu=0;icumin=0;icumax=0
                        for reg in EU[cou]:
                            NP=reg[6]
                            newCum +=int(NP[i])
                            icu +=int(NP[i]*0.09)
                            icumin +=int(NP[i]*0.05)
                            icumax +=int(NP[i]*0.15)
        
                        NEWPOS.append(newCum)
                        ICUS.append(icu)
                        ICUMin.append(icumin)
                        ICUMax.append(icumax)
                        MAXICU.append(maxICU)
                    for rec in EUData:
                        if rec[0]==cou:
                            TIME_DATA=rec[1]
                            NEWPOS_DATA=rec[3]
                            ICUS_DATA=rec[4]

                    fname=dire+ "/"+cou.replace(' ','_') +"_newInfe.jpg"
                    plotFig(dates1=TIME_DATA,q1=NEWPOS_DATA,lab1='Data '+cou, mark1=4,mark2=1,dates2=TIME,q2=NEWPOS,lab2='Forecast '+cou,
                        ymin=0,ymax=np.max(NEWPOS_DATA)*2,
                        #title=cou+': New Positive Quantities',saveFile='')
                        title=cou+': New Positive Quantities',saveFile=fname)
                    
                    fname=dire+ "/"+cou.replace(' ','_') +"_np14.jpg"
                    plotFig(dates1=TIME,q1=CUM_14days,lab1='NP14_100000 '+cou,
                        #title=cou+': New Positive Quantities',saveFile='')
                        title=cou+': New Pos 14 days per 100000 pop',saveFile=fname)
                    
                    if (cou=='France' or cou=='Denmark' or cou=='Greece' or cou=='Slovenia'  or cou=='Iceland') and STUDY=='BYCOUNTRY':
                        if c0 !=-1:
                            direREG=dire0+'/'+DAYFORE+'/'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+'_BYREGION/'+"{:03d}".format(mc)
                        else:
                            direREG=dire0+'/'+DAYFORE+'/imposed_'+prefix+'_'+'BYREGION/'+"{:03d}".format(mc)
                        if not os.path.exists(direREG):
                          os.makedirs(direREG)
                        fname1=direREG + "/"+cou.replace(' ','_') +"_newInfe.jpg"
                        print('****',fname,fname1)
                        copyfile(fname,fname1)
                    

                 
                    ymax1=max(max(np.max(ICUS_DATA)*2,400),maxICU*1.5)
                    imin=0; imax=0
                    imin2=0;imax2=0
                    for i in range(len(xdata)):
                        if ICUMax[i]>MAXICU[i] and imax==0:
                            #if ICUMax[i]>ICUMax[imax]:
                                imax=i
                        if ICUMin[i]>MAXICU[i] and imin==0:
                            #if ICUMin[i]>ICUMin[imin]:
                                imin=i
                        if ICUMax[i]>MAXICU[i]*2 and imax2==0:
                            #if ICUMax[i]>ICUMax[imax2]:
                                imax2=i
                        if ICUMin[i]>MAXICU[i]*2 and imin2==0:
                            #if ICUMin[i]>ICUMax[imin2]:
                                imin2=i
                    
                    try:
                        fname=dire+ "/"+cou.replace(' ','_') +"_icus.jpg"
                        if imax!=0 and imin!=0:
                            riga=cou+','+format(TIME[imax])+','+format(TIME[imin])+','+format(MAXICU[imax])
                            if imax2!=0 and imin2!=0:
                              riga+=','+format(TIME[imax2])+','+format(TIME[imin2])
                            riga +='\n'
                            fr=open(dire+'/times_icu_bycou.csv','a', encoding="utf-8")
                            fr.write(riga)
                            fr.close()
                        plotFig(dates1=TIME_DATA,q1=ICUS_DATA,lab1='Data '+cou, mark1=4,mark2=1,dates2=TIME,q2=ICUS,lab2='Forecast '+cou,
                            q3=ICUMin,lab3='Min',q4=ICUMax,lab4='Max',mark3=1,mark4=1, q5=MAXICU,lab5='MAX ICUs',dates3=TIME,dates4=TIME,dates5=TIME,
                            ymin=0,ymax=ymax1,
                            title=cou+': ICU occupancy',saveFile=fname)
                    except:
                        print('error in '+cou)

                    if (cou=='France' or cou=='Denmark' or cou=='Greece' or cou=='Slovenia' or cou=='Iceland') and STUDY=='BYCOUNTRY':
                        if c0 !=-1:
                            direREG=dire0+'/'+DAYFORE+'/'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+'_BYREGION/'+"{:03d}".format(mc)
                        else:
                            direREG=dire0+'/'+DAYFORE+'/imposed_'+prefix+'_'+'BYREGION/'+"{:03d}".format(mc)
                        if not os.path.exists(direREG):
                          os.makedirs(direREG)
                        fname1=direREG + "/"+cou.replace(' ','_') +"_icus.jpg"
                        print('****',fname,fname1)
                        copyfile(fname,fname1)
                    
                    
                    
                if STUDY=='BYREGION'  and False:
                    for cou in EU:
                        if not os.path.exists(dire+ "/"+cou):
                            os.makedirs(dire+ "/"+cou)
                        for reg in EU[cou]:
                            NEWPOS=reg[6]
            #                NP=reg[6]
            #                for i in range(len(xdata)):
            #                    NEWPOS.append(NP[i])
                            if reg[0] in regions:
                                TIME_DATA=regions[reg[0]][1]
                                NEWPOS_DATA=regions[reg[0]][3]
                                        #ICUS_DATA=rec[4]
                            else:
                                TIME_DATA=[]
                                NEWPOS_DATA=[]
                    
                            fname=dire+ "/"+cou.replace(' ','_')+'/'+reg[0] +"_newInfe.jpg"
                            plotFig(dates1=TIME_DATA,q1=NEWPOS_DATA,lab1='Data '+cou, mark1=4,mark2=1,dates2=TIME,q2=NEWPOS,lab2='Forecast '+cou,
                                ymin=0,ymax=np.max(NEWPOS_DATA)*2,
                                #title=cou+': New Positive Quantities',saveFile='')
                                title=cou+' '+ reg[0]+': New Positive Quantities',saveFile=fname)
        
            except Exception as e:
                print(e)
    #-------!----  MC LOOP ------------
            #end of mc loop
    #---!-----------STUDY LOOP---------
        if mcMax>1: 
             analisi(case,case,dire0,DAYFORE)
        flog.close()
