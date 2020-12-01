#!/usr/bin/env python3

import os
from policiesScenario import mainScenario
#from Epidemic_Model_v0710.main_wrapper import mainOptimize
from datetime import datetime
from analisi import analisiGlob,getAllRt, getPopRegion, getVariabilityRegions

import socket
machine=socket.gethostname()
print(machine)
if machine=='HV-Covid_local':
    system='LINUX_LOCAL'
elif machine=='HV-Covid':
    system='LINUX_JRC'
elif machine=='ECMLAA-NB01':
    system='WINDOWS'
else:
    print('system not found')
    end
   
print ('----  system: '+system)
if system=='WINDOWS':
    #Windows
    dire0='E:/CV/Modelling_Activity/python/forePolicies/'
    dirExecEpi='E:/CV/Modelling_Activity/python/Epidemic Modelling v0710/'
    direPub='E:/CV/Modelling_Activity/FACTSHEETS/'
elif system=='LINUX_JRC':
    dirExecEpi='/home/critechuser/Epidemic_Modelling/python/Epidemic Modelling v0710/'
    dire0='/home/critechuser/Epidemic_Modelling/python/forePolicies/'
    direPub='/mnt/output/CV/FACTSHEETS_forecast/'
elif system=='LINUX_LOCAL':
    dirExecEpi='/mnt/diske/CV/Modelling_Activity/python/Epidemic Modelling v0710/'
    dire0='/mnt/diske/CV/Modelling_Activity/python/forePolicies/'
    direPub='/home/critechuser/Epidemic_Modelling/FACTSHEETS/'

print ('*******************************************************************************')
print ('**                Start fitting ' + datetime.utcnow().strftime('%Y-%m-%d %H:%M'))
print ('*******************************************************************************')
DAYFORE=datetime.utcnow().strftime('%Y-%m-%d')
SpecialCases=True

if SpecialCases:
    DAYFORE='2020-11-05'
    DAYFORE='2020-10-29'
    DAYFORE='2020-10-05'
    DAYFORE='2020-10-19'
    #DAYFORE='2020-09-29_final_0.95'
    DAYFORE='20200929_final_0.95_R'
    DAYFORE='20200929_final_0.95_RR'
    force=True
else:
    force=False

direIn=dirExecEpi+"_report/"+DAYFORE

c0=103;c1=103;analisiGlob(dire0,DAYFORE,c0,c1)
#getPopRegion(direIn, 'Hungary')
#getVariabilityRegions(direIn)

#getAllRt(dire0)

calibrationPeriod=30   # 30

fname1=direIn+ '/FITTING_ALLCOUNTRIES/data/FITTING_ALLCOUNTRIES_SIR_(-'+format(calibrationPeriod)+' -1).json'
fname2=direIn+ '/FITTING_ALLREGIONS/data/FITTING_ALLREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).json'

if not (os.path.exists(fname1) and os.path.exists(fname2))  and not force:
   print('calling optimization')
   # mainOptimize()
   # execute Epidemic model v0710
   cmd='cd "'+dirExecEpi+'";python3  main_wrapper.py -c '+format(calibrationPeriod)+' -s 0;cd ../forePolicies'
   print(cmd)
   os.system(cmd)


print('**                Start scenario ' + datetime.utcnow().strftime('%Y-%m-%d %H:%M'))
print('*******************************************************************************')

if (os.path.exists(fname1) and os.path.exists(fname2)) or force:
    if SpecialCases:
        #  FILTER   Country: [lockdown date],[target R0],[lockmax],[lockmin]
        #prefix='case1';    FILTER_COUNTRY={'Italy':['2020-11-06',0.7,0,400], 'France':['2020-10-30',0.7,0,400]}
        #prefix='case2';    FILTER_COUNTRY={'Italy':['2020-11-13',0.7,0,400], 'France':['2020-10-30',0.95,0,400]}
        #prefix='case3';    FILTER_COUNTRY={'Italy':['2020-11-06',0.95,0,400]}
        #prefix='case4';    FILTER_COUNTRY={'Italy':['2020-11-03',0.7,0,400]}
        prefix='';FILTER_COUNTRY={}
        #c0=-1;c1=-1
        #c0=105;c1=107
        #c0=103;c1=103
        #c0=-1;c1=-1;prefix='lock_50perc_ICU';FILTER_COUNTRY={'Italy':['2020-11-06',0.7,-1,-1,50,50]}
        #c0=-1;c1=-1;prefix='lock_20_20perc_wait14day';FILTER_COUNTRY={'Italy':['2020-10-19',0.7,-1,-1,20,20,14]}
        #c0=-1;c1=-1;prefix='lock_50_20perc_wait7day';FILTER_COUNTRY={'Italy':['2020-10-19',0.7,-1,-1,50,20,7]}
        #c0=-1;c1=-1;prefix='lock_50_20perc_wait14day';FILTER_COUNTRY={'Italy':['2020-10-19',0.7,-1,-1,50,20,14]}
        #c0=-1;c1=-1;prefix='imposed_lock_25_25_wait_60_0.7';FILTER_COUNTRY={'Italy':['2020-10-19',0.7,-1,-1,25,25,60]}
        #mainScenario(dire0,direIn,DAYFORE,c0,c1,FILTER_COUNTRY,prefix,calibrationPeriod)
        #for lockPerc in [25,50,75]:
        #    for unlockPerc in [lockPerc,lockPerc-25,lockPerc-50]:
        #        if unlockPerc<0:
        #                continue
        #        for dayswait in [2,7,14,30,60]:
        #           for targetRt in [0.7, 0.95, 1.2]:
        #            dire1=dire0+'/'+DAYFORE+'/imposed_lock_'+format(lockPerc)+'_'+format(unlockPerc)+'_wait_'+format(dayswait)+'_'+format(targetRt)+'_BYREGION'
        #            if not os.path.exists(dire1+'/000/Italy_newInfe.jpg'):
        #                c0=-1;c1=-1;prefix='lock_'+format(lockPerc)+'_'+format(unlockPerc)+'_wait_'+format(dayswait)+'_'+format(targetRt);FILTER_COUNTRY={'Italy':['2020-10-19',targetRt,-1,-1,lockPerc,unlockPerc,dayswait]}
        #                mainScenario(dire0,direIn,DAYFORE,c0,c1,FILTER_COUNTRY,prefix,calibrationPeriod)
        c0=103;c1=103;FILTER_COUNTRY='';prefix=''
        mainScenario(dire0,direIn,DAYFORE,c0,c1,FILTER_COUNTRY,prefix,calibrationPeriod)
    else:
        c0=100;c1=100;FILTER_COUNTRY='';prefix=''
        mainScenario(dire0,direIn,DAYFORE,c0,c1,FILTER_COUNTRY,prefix,calibrationPeriod)
    
    
    
    # publication
    direPubDay=direPub+datetime.utcnow().strftime('%Y%m%d')+'/forecast'
    if not os.path.exists(direPubDay):
      os.makedirs(direPubDay)
    readDir=dire0+DAYFORE+'/1000000_1000000_control_at_0_0.95_BYREGION/000'
    os.system('cp '+readDir+'/* '+direPubDay)
    
    
else:
    #print (fname1,fname2, 'one of those files not existing')
    print ('run optimization first')

print ('*******************************************************************************')
print ('**                End job  ' + datetime.utcnow().strftime('%Y-%m-%d %H:%M'))
print ('*******************************************************************************')
