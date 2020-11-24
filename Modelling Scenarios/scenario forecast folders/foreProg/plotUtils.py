import os, sys, getopt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
from matplotlib.colors import ListedColormap
from matplotlib.dates import date2num, num2date
from matplotlib import ticker
import pandas as pd
import datetime
import numpy as np
import os 
from scipy import stats as sps
from scipy.interpolate import interp1d



def plotFig( **kwargs):
    
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    days = mdates.DayLocator()
    hours = mdates.HourLocator()
    dayFmt = mdates.DateFormatter('%d/%m')

    fig, ax = plt.subplots()
    # format the ticks
    
    
    colori=['b','r','o','g','m']
    title=kwargs.pop('title','')
    lab1=kwargs.pop('lab1','' )
    lab2=kwargs.pop('lab2','' )
    lab3=kwargs.pop('lab3','' )
    lab4=kwargs.pop('lab4','' )
    lab5=kwargs.pop('lab5','' )

    mark1=kwargs.pop('mark1',2 )
    mark2=kwargs.pop('mark2',2 )
    mark3=kwargs.pop('mark3',2 )
    mark4=kwargs.pop('mark4',2 )
    mark5=kwargs.pop('mark5',2 )

    dates1=kwargs.pop('dates1','')
    dates2=kwargs.pop('dates2','')
    dates3=kwargs.pop('dates3','')
    dates4=kwargs.pop('dates4','')
    dates5=kwargs.pop('dates5','')
    q1=kwargs.pop('q1','')
    q2=kwargs.pop('q2','')
    q3=kwargs.pop('q3','')
    q4=kwargs.pop('q4','')
    q5=kwargs.pop('q5','')
    
    ylab=kwargs.pop('ylab','Individuals')
    xlab=kwargs.pop('xlab','')
    outfile=kwargs.pop('out','test.jpg')
    xdate=kwargs.pop('xdate',True)
    xlog=kwargs.pop('xlog',False)
    ylog=kwargs.pop('ylog',False)
    
    xmin=kwargs.pop('xmin',-1)
    xmax=kwargs.pop('xmax',-1)
    ymin=kwargs.pop('ymin',-1)
    ymax=kwargs.pop('ymax',-1)
    
    areafill=kwargs.pop('areafill',False)
    setLine=kwargs.pop('setLine','')
    
    saveFile=kwargs.pop('saveFile','')
        
    if lab1 !='':
        ax.plot(dates1,q1,color='blue', marker='o',linewidth=2,markersize=mark1,label=lab1)
    if lab2 !='':
        ax.plot(dates2,q2,color='red', marker='o',linewidth=1,markersize=mark2,label=lab2)
    if lab3 !='':
        ax.plot(dates3,q3,color='orange', marker='o',linewidth=1,markersize=mark3,label=lab3)
    if lab4 !='':
        ax.plot(dates4,q4,color='green', marker='o',linewidth=1,markersize=mark4,label=lab4)
    if lab5 !='':
        ax.plot(dates5,q5,color='black',linestyle='dotted', marker='',linewidth=2,markersize=mark5,label=lab5)

    if not (ymin==-1 and ymax==-1):
        plt.ylim(ymin,ymax)
    
    if not (xmin==-1 and xmax==-1):
        plt.xlim(xmin,xmax)    
       
    if xdate:
        ax.xaxis.set_major_locator(days)
        ax.xaxis.set_major_formatter(dayFmt)
        ax.format_xdata = mdates.DateFormatter('%d %m')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, 30))

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
        
    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

    if areafill and lab2 !='':
        q2i=np.interp(date2num(dates1),date2num(dates2),q2)
        vmax=np.maximum(q1,q2i)
        vmin=np.minimum(q1,q2i)
        if lab3 != '':
            q3i=np.interp(date2num(dates1),date2num(dates3),q3)
            vmax=np.maximum(vmax,q3i)
            vmin=np.minimum(vmin,q3i)

        if lab4 != '':
            q4i=np.interp(date2num(dates1),date2num(dates4),q4)
            vmax=np.maximum(vmax,q4i)
            vmin=np.minimum(vmin,q4i)

        if lab5 != '':
            q5i=np.interp(date2num(dates1),date2num(dates5),q5)
            vmax=np.maximum(vmax,q5i)
            vmin=np.minimum(vmin,q5i)
            
        vmax *=1.2
        vmin *=0.8
        ax.fill_between(dates1,
                    vmin,vmax,
                    
                    color='k',
                    alpha=.1,
                    lw=0,
                    zorder=3)

        


    ax.set_title(title,fontsize=10,fontweight="bold")
        
    plt.ylabel(ylab)
    if xlab !='':
        plt.xlabel(xlab)
    #plt.legend(bbox_to_anchor=(0.5,-0.22), loc='center', ncol=4,fontsize=10,frameon=False)
    ax.legend()
    
        
    xmin,xmax=ax.get_xlim()
    ymin,ymax=ax.get_ylim()
    
    if setLine !='':
        line1=[];xline1=[]
        xline1.append(xmin);line1.append(setLine)
        xline1.append(xmax);line1.append(setLine)
        
        ax.plot(xline1,line1,color='red',linewidth=2, linestyle='dotted')
        ax.set_xlim(xmin,xmax)
    
    #ig,ax=plt.plot(dates,infected, fatalities,currPos,icus)
    h=200;w=400
    h=h/45
    w=w/45
    fig.set_size_inches(w,h)

    if 'Sweden'in saveFile:
        saveFile=saveFile
    print ('creating ',saveFile)
    if saveFile=='':
        plt.show()
    else:
        plt.savefig(saveFile,transparent=False)
    plt.clf()
    plt.close('all')
    return xmin,xmax,ymin,ymax
 
def createFigures(cou,dirout,factdir,rows,rundate): 
    print ('\n****************************')
    print ('***           '+cou+'        ***')
    print ('*****************************')
    if not os.path.exists(dirout):
        os.makedirs(dirout)
    now=rundate #datetime.datetime.now().strftime("%Y%m%d")    
    print('rundate=',rundate)    
   # if not os.path.exists(dirout+'\\'+now):
   #     os.makedirs(dirout+'\\'+now)
        
     #1.  leggere i dati
  
    if cou=='Czech_Republic': cou='Czech Republic' #PB
    infected=[];fatalities=[];recovered=[];dates=[];icus=[];hospitalized=[];currPos=[]
    for r in rows:

        if r=='': continue
        p=r.split(",")
        if p[2]==cou.strip():            
            #rint (r)
            if p[5]!='0':
                #rint (r)
                #print (p[2]+' '+p[0])
                d=datetime.datetime.strptime(p[0],'%Y-%m-%d')
                dates.append(d)
                for ij in range(5,11):
                   if p[ij]=='': p[ij]=0
                infected.append(int(p[5]))
                fatalities.append(int(p[6]))
                recovered.append(int(p[7]))
                icus.append(int(p[10]))
                hospitalized.append(int(p[9]))
                currPos.append(int(p[5])-int(p[6])-int(p[7]))

    vinfected=[];vfatalities=[];vCurrPos=[];vicus=[];vdates=[]
    ninfected=[];nfatalities=[];nicus=[];ndates=[]
    for k in range(len(dates)):
        if k>0:
            vdates.append(dates[k])
            vinfected.append((infected[k]-infected[k-1])/gd(dates[k],dates[k-1]))
            vfatalities.append((fatalities[k]-fatalities[k-1])/gd(dates[k],dates[k-1]))
            vCurrPos.append((currPos[k]-currPos[k-1])/gd(dates[k],dates[k-1]))
            vicus.append((icus[k]-icus[k-1])/gd(dates[k],dates[k-1]))
        if k>=7:
            #ndates.append(dates[k])
            ninfected.append(infected[k]-infected[k-7])
            nfatalities.append(fatalities[k]-fatalities[k-7])
            nicus.append(icus[k]-icus[k-7])

    for k in range(len(ninfected)):
         if ninfected[k]<=0:ninfected[k]=1
         if nfatalities[k]<=0:nfatalities[k]=1
         if nicus[k]<=0:nicus[k]=1

    #calculation of r
    rv_cases=[];rv_fatalities=[];rt_icus=[]
    xrv_cases=[];xrv_fatalities=[];xrt_icus=[]

    rv_casesS=[];rv_fatalitiesS=[];rt_icusS=[];dates_rS=[]
    xrv_casesS=[];xrv_fatalitiesS=[];xrt_icusS=[];dates_rS=[]

    for k in range(8,len(vinfected)):
        #print (k,len(ninfected))
        if vinfected[k]>0 and vinfected[k-7]>0:
            #print(np.log(ninfected[k-1]))
            #print(np.log(ninfected[k-1-7])) 
            #print(vdates[k],vdates[k-7])
            rv_cases.append ((np.log(vinfected[k])-np.log(vinfected[k-7])) /gd(vdates[k],vdates[k-7])*7.0+1.0)
            xrv_cases.append(vdates[k])
           #print (vdates[k])

        if vfatalities[k]>0 and vfatalities[k-7]>0:
            rv_fatalities.append ((np.log(vfatalities[k])-np.log(vfatalities[k-7]))/gd(vdates[k],vdates[k-7])*7.0+1.0)
            xrv_fatalities.append(vdates[k])

    for k in range(8,len(icus)):
        if icus[k]>0 and icus[k-7]>0:
            rt_icus.append((np.log(icus[k])-np.log(icus[k-7]))/gd(dates[k],dates[k-7])*7.0+1.0)
            xrt_icus.append(dates[k])
            
    for k in range(len(rv_cases)):
        rv_casesS.append(solve2(k-14,k,rv_cases,xrv_cases,xrv_cases[k]))
        xrv_casesS.append(xrv_cases[k])

    for k in range(len(rv_fatalities)):    
        rv_fatalitiesS.append(solve2(k-14,k,rv_fatalities,xrv_fatalities,xrv_fatalities[k]))
        xrv_fatalitiesS.append(xrv_fatalities[k])
        
    for k in range(len(rt_icus)):    
        rt_icusS.append(solve2(k-14,k,rt_icus,dates,dates[k]))
        xrt_icusS.append(xrt_icus[k])
                            
    #if cou=='Czech_Republic': cou=cou.replace("_"," ")
    xr0,r0,rl,rh, csv=calcR0(cou)
    fname=dirout+ "\\" + now +"_"+cou+"_RtKs.jpg"
    plotRt(cou,xr0,r0,rl,rh,fname,'K. Systrom')    
    
    #readFolder='E:\\CV\\FACTSHEETS\\scripts_py3\\test\\NAT'
    readFolder=factdir+'\\NAT'
    
    xr0C,r0C,rlC,rhC=calcR0_CRAN(cou,readFolder)
    fname=dirout+ "\\" + now +"_"+cou+"_RtCRAN.jpg"
    #print (xr0C,r0C)
    plotRt(cou,xr0C,r0C,rlC,rhC,fname,'CRAN')    
    
    xr0J,r0J,rlJ,rhJ=calcR0_JRC(cou)
    fname=dirout+ "\\" + now +"_"+cou+"_RtJRC.jpg"
 #   print (xr0J,r0J)
    plotRt(cou,xr0J,r0J,rlJ,rhJ,fname,'JRC')    
    
    xr0RKI,r0RKI=calcR0_RKI(cou)
    fname=dirout+ "\\" + now +"_"+cou+"_RtRKI.jpg"
 #   print (xr0J,r0J)
    plotRt(cou,xr0RKI,r0RKI,r0RKI,r0RKI,fname,'RKI')    
    
    # f=open(dirout+"\\"+now+"_"+cou+"_r0.csv","w")
    # f.write(csv.replace("\n",""))
    # f.close()
    if cou=='Czech Republic': cou='Czech_Republic'   
    print(cou)    
    #print(r0)
    print('rundate2=',rundate)
    print(now,dirout)
    fname=dirout+ "\\" + now +"_"+cou+"_Cumulative.jpg"
    plotFig(dates1=dates,q1=infected,lab1='Positive',dates2=dates,q2=fatalities,lab2='Fatalities',
             dates3=dates,q3=currPos,lab3='Current Positive',dates4=dates,q4=icus,lab4='ICUs',
             title=cou+': Cumulative Quantities',saveFile=fname)

    fname=dirout+ "\\" + now+"_"+cou+ "_Daily.jpg"    
    plotFig(dates1=vdates,q1=vinfected,lab1='Daily Positive',dates2=vdates,q2=vfatalities,lab2='Daily Fatalities',
             dates3=vdates,q3=vCurrPos,lab3='Daily Current Positive',dates4=vdates,q4=vicus,lab4='Daily ICUs',
             title=cou+': Daily Quantities',saveFile=fname)
    
    fname=dirout+ "\\" + now +"_"+cou+ "_Epidemic.jpg"    
    plotFig(dates1=infected[7:],q1=ninfected,lab1='New Pos. last week',dates2=fatalities[7:],q2=nfatalities,lab2='New Fatalities last week',
            dates3=icus[7:],q3=nicus,lab3='New ICUs last week',mark2=0,mark3=0,mark4=0,mark5=0,
            title=cou+': Epidemic Status',xlog=True,ylog=True,xdate=False,xlab='Overall number of cases/fatalities/icus',saveFile=fname)

#    fname=dirout+ "\\"+ now+"_"+cou+ "_ReprNumber.jpg"    
#    plotFig(dates1=xrv_cases,q1=rv_cases,lab1='rv_cases',dates2=xrv_fatalities,q2=rv_fatalities,lab2='rv_fatalities',
#            dates3=xrt_icus,q3=rt_icus,lab3='rt_icus',
#            title=cou+': Reproduction Number',ymin=0.0,ymax=4.5, setLine=1,saveFile=fname)
    
#    fname=dirout+ "\\" + now+"_"+cou+ "_ReprNumberS.jpg"    
#    plotFig(dates1=xrv_casesS,q1=rv_casesS,lab1='rv_cases (S)',dates2=xrv_fatalitiesS,q2=rv_fatalitiesS,lab2='rv_fatalities (S)',
#            dates3=xrt_icusS,q3=rt_icus,lab3='rt_icus (S)',
#            title=cou+': Reproduction Number (S)',ymin=0.0,ymax=4.5, setLine=1,saveFile=fname)

    fname=dirout+ "\\" + now+"_"+cou+ "_ReprNumberRt.jpg"

    #print (np.shape(xr0),np.shape(r0))
    plotFig(dates1=xr0,q1=r0,lab1='Rt KS',dates2=xr0J,q2=r0J,lab2='rt JRC',
            dates3=xr0C,q3=r0C,lab3='Rt CRAN',dates4=xr0RKI,q4=r0RKI,lab4='Rt RKI',
            title=cou+': Reproduction Number',ymin=0.0,ymax=3., setLine=1,saveFile=fname, areafill=True)

