from analisi import getrows
import numpy as np
import matplotlib.pyplot as plt
import os,random
from shutil import copyfile
from datetime import datetime
from datetime import timedelta

def transF(t,tlock,delta,R_0_start,R_0_end):
    if t<tlock:
        return R_0_start
    elif t>tlock+delta*2:
        return R_0_end
    else:
        ff=10/delta
        return (R_0_start-R_0_end)/(1+np.exp(-ff*(-t+tlock+delta/2)))+R_0_end

def polyNational(direIn,dire,data, STUDY,ecList,lockMax,lockMin,R0Target,StartTimeControl,fmc1,fmc2):

    regions=getrows(direIn)
    EU={}
    for cou in ecList:

        iculockMax=1e9
        iculockMin=0
        tunlock=-1e6
        tlock=-1e6
        waitTime=7
        lock=False
        print ('---------------'+cou+'--------------------------')
        regData=[]
        for (reg,v) in data.items():
            try:
                #print (reg)
                couReg=regions[reg][0]
            except:
                couReg=''
                continue
            #print(couReg)
            if cou==couReg:
                r0=data[reg]['r0']+fmc1
                if r0<0.1: r0=0.1
                Trecov=data[reg]['Trecov']+fmc2
                #flog.write(format(mc)+','+dire+','+format(fmc1)+','+format(fmc2)+','+format(r0)+','+format(Trecov)+'\n')
                N=data[reg]['population']
                R=data[reg]['_ivp']['R']
                I=data[reg]['_ivp']['I']
                S=N-I-R
                time0=data[reg]['time'][0]
                I0=I
                daysLock=0
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
                r00=r0
                R0Lock=r00
                R0UnLock=r00

                regData.append([cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0Lock,R0UnLock])

                print(cou,reg)


                #print(cou,reg,r0,Trecov,N)
    
        dt=0.1
        lockImplementation=20  # implemenattion of measures
        tlock=0
        tunlock=-1e6
        lock=False
        unlocking=False
        unlocked=False
        CPNAT=[]
        NEWPONAT=[]
        firstLock=True

        xdata=np.linspace(0,30*6*int(1/dt),30*6*int(1/dt))*dt  # 6 months forecast
        first = True
        #f=open(dire+cou+'.csv','w')
        if len(regData)>0:
    #-------|--- for time
            for i in range(len(xdata)):
                print(i)
                cpnat=0
                nponat=0
                natPop=0
    #-----------|--- for regions
                for k in range(len(regData)):
                    regd=regData[k]
                    cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0Lock,R0UnLock=regd
                    ' ricalcolo r0 based on lock/unlock'
                    if lock: # and xdata[i]>daysAfterLock+tlock:
                        if firstLock:
                            firstLock=False
                        deltaDelay=7
                        r01=transF(xdata[i],tlock+deltaDelay,lockImplementation,R0Lock,R0Target)
                    else:
                        deltaDelay=2   # immediate releasing of people

                        r01=transF(xdata[i],tunlock+deltaDelay,lockImplementation,R0UnLock,r00)
                    print(reg,r0,r01,R0Lock,R0UnLock)
                    r0=r01
                    dSdt=-r0/Trecov*S*I/N
                    dIdt=r0/Trecov*S*I/N - 1/Trecov*I
                    dRdt=1/Trecov*I
        
                    S +=dSdt*dt
                    I +=dIdt*dt
                    R +=dRdt*dt
                    CP =N-S
                    cpnat +=CP
                    natPop +=N

                    II.append(I)
                    RR.append(R)
                    SS.append(S)
                    CUMPOS.append(CP)
                    npo=(CP-CP0)/dt
                    nponat +=npo
                    NEWPOS.append(npo)
                    CP0=CP
                    regd=cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0Lock,R0UnLock
                    regData[k]=regd
                
                    #print(format(xdata[i])+','+format(S)+','+format(I)+','+format(R)+','+format(CP)+','+format(r0)+','+format(npo))
        
                    #  I     =  CURRENT POSITIVE !!!
                    #  -dSdt =  NEW POSITIVE
                    #  N-S   =  CUMULATIVE POSITIVE
    #-----------|--- end of for regions            
                CPNAT.append(cpnat)
                NEWPONAT.append(nponat)
            
                if len(CUMPOS)<15/dt:
                    cumIncidence_14days=0
                else:
                    if STUDY=='BYCOUNTRYNAT':
                        cumIncidence_14days=(cpnat-CPNAT[len(CPNAT)-int(14/dt)])*100000.0/natPop          
                    elif STUDY=='BYCOUNTRYNAT_ANYREG':
                        # cerca il cum piu'alto nelle regioni
                        cumIncidence_14days=0
                        for k in range(len(regData)):
                            regd=regData[k]
                            cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0lock,R0unlock=regd
                            spi=(CUMPOS[len(CUMPOS)-1]-CUMPOS[len(CUMPOS)-int(14/dt)])*100000.0/N           
                            if spi>cumIncidence_14days:
                                cumIncidence_14days=spi
                    CUM_14days.append(cumIncidence_14days)
                    ICUEstimate=nponat*0.09
                    if xdata[i]>StartTimeControl:
                        if (cumIncidence_14days>lockMax or (ICUEstimate>iculockMax and xdata[i]-tunlock>waitTime)) and not lock :  # and (I0<lockMax):
                                
                            tlock=xdata[i]
                            lock=True
                        
                            for k in range(len(regData)):
                                regd=regData[k]
                                cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0Lock,R0UnLock=regd
                                R0Lock=r0
                                regd=cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0Lock,R0UnLock
                                regData[k]=regd
                        if (cumIncidence_14days<lockMin or (ICUEstimate<iculockMin and xdata[i]-tlock>waitTime)) and lock :
                            tunlock=xdata[i]
                            lock=False
                            for k in range(len(regData)):
                                regd=regData[k]
                                cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0Lock,R0UnLock=regd
                                R0UnLock=r0
                                regd=cou,reg,r0,Trecov,N,S,R,I,time0,II,RR,SS,CUMPOS,NEWPOS,LOCK,REPR,CUM_14days,r00,CP0,R0Lock,R0UnLock
                                regData[k]=regd
        
                    #print (xdata[i],r0,tlock,tunlock)
                REPR.append(r0)
                if lock:
                    LOCK.append(1)
                else:
                    LOCK.append(0)
        
                    #print (i,r0, R0Lock,R0UnLock)
    #-------|--- for time end
                #f.close()
            if not cou in EU:
                EU[cou]=[]

            EU[cou].append([cou,CPNAT,N,LOCK,NEWPONAT]) #,CUM_14days])
                
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
    try:
        for c in EU:
            cc.write(','+c)
            Ncc.write(','+c)
            lp.write(','+c)
            #lt.write(',Date,'+c)
            for regs in EU[c]:
                reg,CPNAT,N,LOCK,NEWPONAT=regs
                if STUDY=='BYREGION':
                    rr.write(','+regs[0]+'@'+c)
                    ls.write(','+regs[0]+'@'+c)
                    Nrr.write(','+regs[0]+'@'+c)
                else:
                    rr.write(','+reg)
                    ls.write(','+reg)

                    Nrr.write(','+reg)



        lt.write('Country,TotDays,TotPopLock,Population\n')
        lp.write(',TOT\n')
        cc.write(',TOT\n')
        rr.write(',TOT\n')
        Ncc.write(',TOT\n')
        Nrr.write(',TOT\n')
        rt.write('\n')
    except Exception as e:
        print(e)
    try:
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
                
                for reg,CUMPOS,N,LOCK,NEWPOS in EU[cou]:
                    # [reg,II,SS,CPCP,N,LOCK,NEWPOS,REPR])
                    #   0   1  2  3   4  5
                    #II=reg[1]
                    #totcou +=II[i]
                    #tot+=II[i]
                    #rr.write(','+format(int(reg[1][i])))
                    #SS=reg[2]
                    #CUMPOS=reg[3]
                    #NEWPOS=reg[6]
                    #LOCK=reg[5]
                    totcouCum +=CUMPOS[i]
                    totCum +=totcouCum
                    rr.write(','+format(int(CUMPOS[i])))
                    ls.write(','+format(LOCK[i]))
                    #rt.write(','+format(REPR[i]))
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
    except Exception as e:
        print(e)

    for cou in EU:
        TotdaysLock=0
        TotpopLock=0
                
        for i in range(len(xdata)):
            atleastOnelocked=False
            TotpopLock=0
            totN=0
            for reg,CPNAT,N,LOCK,NEWPONAT in EU[cou]:
                    
                totN +=N
                    
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
            
        #TIME=[]
        #for i in range(len(xdata)):
        #    ti=datetime.strptime(time0,'%Y/%m/%d')+timedelta(days=xdata[i])
        #    TIME.append(ti)
        #fr=open(dire+'/times_icu_bycou.csv','w', encoding="utf-8")
        #fr.write('Country, dateMin, DateMax, ICUMax\n')
        #fr.close()

