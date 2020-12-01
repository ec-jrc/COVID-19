import random
import pandas as pd
import numpy as np
import datetime
import json
import numpy as np
import matplotlib.pyplot as plt
import os,random
from shutil import copyfile
from datetime import datetime
from datetime import timedelta

def prepare_cases(cases, cutoff=25,step=14):
    new_cases = cases.diff()
    smoothed = new_cases.rolling(step,
        win_type='gaussian',
        min_periods=1,
        center=True).mean(std=2).round()
    
    idx_start = np.searchsorted(smoothed, cutoff)
    smoothed = smoothed.iloc[idx_start:]
    original = new_cases.loc[smoothed.index]
    #if len(smoothed)==0:
    #    if step==7:
    #        original,smoothed=prepare_cases(cases,cutoff,1)
    #    else:
    #        original=cases
    #        smoothed=cases
    if len(smoothed)==0:
     cutoff -=1
     #print (cutoff)
     if cutoff>0:
         original,smoothed=prepare_cases(cases,cutoff,step)
     else:
         original=cases
         smoothed=cases                
    return original, smoothed

def getDirCase(dire0,DAYFORE,STUDY,case):
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
  #================================================  
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
#================================================
        elif case==110:
            lockMax=-1
            lockMin=-1
            StartTimeControl=31
            R0Target=0.7

        #dire1=dire0+'/'+DAYFORE+'/CASE '+format(case-100)+'  '+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+"_"+STUDY
        dire1=dire0+'/'+DAYFORE+'/'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+"_"+STUDY
        return dire1

def fvals(rowsCum,rowsNew,rowsLock,cou0):
		for i in range (len(rowsCum[0].split(','))):
			cou=rowsCum[0].split(',')[i]
			if cou=='days' or cou=='TOT\n' or cou=='Date':
					continue
			if '@' in cou:
				cou=cou.split("@")[1]
			if cou==cou0:
				vCum=float(rowsCum[1800].split(',')[i])-float(rowsCum[302].split(',')[i])
				maxNew=0
				totLock=0
				
				for j in range(1,len(rowsNew)):
					
					if maxNew<float(rowsNew[j].split(',')[i]):
						maxNew=float(rowsNew[j].split(',')[i])
				for j in range(len(rowsLock)):
					cou=rowsLock[j].split(',')[0]
					if cou==cou0:
						daysLock=float(rowsLock[j].split(',')[1])
						pop=float(rowsLock[j].split(',')[3])
						totLock=round(daysLock/pop,0)
						break
				return vCum,maxNew,totLock
		return 0,0,0
def getMaxEU27(rowsNew,maxNewFrance,maxNewDenmark,maxNewGreece,rowsCouNew):
		maxEU27=0
		timeline=[]
		ecList = "Austria, Belgium, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg, Malta, Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia, Spain, Sweden".split(', ')
		for j in range(1,len(rowsCouNew)):
			maxValue=0

			for i in range (len(rowsNew[0].split(','))):
				cou=rowsNew[0].split(',')[i]
				if cou=='days' or cou=='TOT\n' or cou=='Date':
					continue
				if not cou in ecList:
					continue
				if cou=='France' or cou=='Denmark' or cou=='Greece':
					#vNew=maxNewFrance
					ccs=rowsCouNew[0].split(',')
					for g in range(len(ccs)):
						if ccs[g]==cou:
							vNew=rowsCouNew[j].split(',')[g]
							break
				#elif \cou=='Denmark':
				#	vNew=maxNewDenmark
				#	vNew=rowsCouNew[j].split(',')[i]
				#elif cou=='Greece':
				#	vNew=maxNewGreece			
				#	vNew=rowsCouNew[j].split(',')[i]
				else:
					vNew=rowsNew[j].split(',')[i]
				#print(cou,vNew)
				maxValue+=float(vNew)
			timeline.append(maxValue)
			if maxValue>maxEU27:
				maxEU27=maxValue
			#print (maxValue,maxEU27)
		return maxEU27,timeline

def openFiles(dirOut):
	f=open(dirOut+'\\000\\cumPos_byCountry.csv','r')
	rowsCum=f.readlines()
	f.close
	f=open(dirOut+'\\000\\newPos_byCountry.csv','r')
	rowsNew=f.readlines()
	f.close
	f=open(dirOut+'\\000\\lockTime2.csv','r')
	rowsLock=f.readlines()
	f.close
	return rowsCum, rowsNew, rowsLock


def getParamGlobal(dire0,DAYFORE,case):
	STUDY='BYCOUNTRY';dirOutCou=getDirCase(dire0,DAYFORE,STUDY,case)
	rowsCouCum,rowsCouNew,rowsCouLock=openFiles(dirOutCou)	
	vCumFrance,maxNewFrance,totLockFrance=fvals(rowsCouCum,rowsCouNew,rowsCouLock,'France')
	vCumDenmark,maxNewDenmark,totLockDenmark=fvals(rowsCouCum,rowsCouNew,rowsCouLock,'Denmark')
	vCumGreece,maxNewGreece,totLockGreece=fvals(rowsCouCum,rowsCouNew,rowsCouLock,'Greece')

	STUDY='BYCOUNTRYNAT';dirOutCou=getDirCase(dire0,DAYFORE,STUDY,case)
	rowsCouCum,rowsCouNew,rowsCouLock=openFiles(dirOutCou)	

	#STUDY='BYCOUNTRYNAT_ANYREG';dirOutCou=getDirCase(dire0,DAYFORE,STUDY,case)
	STUDY='BYREGION';dirOutReg=getDirCase(dire0,DAYFORE,STUDY,case)
	rowsRegCum,rowsRegNew,rowsRegLock=openFiles(dirOutReg)	

	avgLock=0
	for i in range (len(rowsCouCum[0].split(','))):
		cou=rowsCouCum[0].split(',')[i]
		
		if cou=='days' or cou=='TOT\n' or cou=='Date':
			continue
		vCumNat,vNewNat,totLockNat=fvals(rowsCouCum,rowsCouNew,rowsCouLock,cou)
		if cou=='France':
			vCum=vCumFrance
			vNew=maxNewFrance
			totLock=totLockFrance
		elif cou=='Denmark':
			vCum=vCumDenmark
			vNew=maxNewDenmark
			totLock=totLockDenmark
		elif cou=='Greece':
			vCum=vCumGreece
			vNew=maxNewGreece
			totLock=totLockGreece
		else:
			vCum,vNew,totLock=fvals(rowsRegCum,rowsRegNew,rowsRegLock,cou)

		print(cou,',',vCum,',',vNew,',',vCumNat,',',vNewNat,',',totLock,',',totLockNat)
		avgLock +=totLock
	avgLock /=27
	avgLock=avgLock /150*100
	maxReg,timelineReg=getMaxEU27(rowsRegNew,maxNewFrance,maxNewDenmark,maxNewGreece,rowsCouNew)
	maxCou,timelineCou=getMaxEU27(rowsCouNew,maxNewFrance,maxNewDenmark,maxNewGreece,rowsCouNew)

	print('maxEU27,,', maxReg,',,',maxCou,',',round(avgLock,0))
	return timelineReg,timelineCou

def analisiGlob(dire0,DAYFORE,c0,c1):
	casi=[]
	for case in range(c0,c1+1):
		print('=====================')
		print('case ',case)
		print('=====================')
		timelineReg,timelineCou=getParamGlobal(dire0,DAYFORE,case)
		casi.append(timelineReg)

	f=open(dire0+'/'+DAYFORE+'/outcasi.txt','w')
	for j in range(len(casi[0])):
		riga=''
		for c in range(c0,c1+1):
			riga += format(casi[c-c0][j])+','
		f.write(riga+'\n')
	f.close()
		

def getMinMax(dire):
		EUtot=[]
		for mc in range(0,50):
			print(mc)
			fname=dire+"BYCOUNTRYNAT/"+"{:03d}".format(mc)+"/newPos_byCountry.csv"
			f=open(fname,'r')
			r=f.read().split('\n')
			f.close()
			fra=[]
			for riga in r:
				if riga=='':
					continue
				fra.append(riga.split(',')[11])

			fname=dire+"BYREGION/"+"{:03d}".format(mc)+"/newPos_byCountry.csv"
			f=open(fname,'r')
			r=f.read().split('\n')
			f.close()
			if mc==0:
				dat=[None]*(len(r)-1)
				max=[-1e10]*(len(r)-1)
				min=[1e10] *(len(r)-1)
				ref=[0]*(len(r)-1)
				value=[0]*(len(r)-1)
			
			for j in range(1,len(r)):
				
				riga=r[j]
				if riga=='':
					continue
				dat[j]=riga.split(',')[1]
				tot=0
				for k in range(2,27):
					tot +=int(riga.split(',')[k])
				tot +=int(fra[j])
				value[j]=tot
				if tot>max[j]:
					max[j]=tot
				if tot<min[j]:
					min[j]=tot
				if mc==0:
					ref[j]=tot
			EUtot.append(value.copy())
		return dat,ref,min,max,EUtot

def analisi(c0,c1,dire0,DAYFORE):
	#dire0='E:/CV/Programmi/epidemic_fitting_vers20200819/forePolicies/'
	#DAYFORE='20200826'
	#STUDY='BYREGION'
	#STUDY='BYCOUNTRY'
	random.seed()
	#rows=getrows(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS)
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
	
		#dire=dire0+'/'+DAYFORE+'/'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"
		dire=dire0+'/'+DAYFORE+'/'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+"_"
		#dire=dire0+'/'+DAYFORE+'/case'+format(int(lockMax)) +"_"+format(int(lockMin))+"_control_at_"+format(StartTimeControl)+"_"+format(R0Target)+"_"
		dire=dire.replace('//','/')
		dat,ref,min,max, EU=getMinMax(dire)
	
		fname=dire+"BYREGION/MC_result.csv"
		f=open(fname,'w')
		for j in range(len(dat)):
			f.write(format(dat[j])+',')
			f.write(format(ref[j])+',')
			f.write(format(min[j])+',')
			f.write(format(max[j])+',,')
			for k in range(len(EU)):
				
				f.write(format(EU[k][j])+',')
			
			f.write('\n')
		f.close()

def extractDemo():
	import json
	ff=open('region_comb.json','r',encoding='UTF-8')
	data=json.load(ff)
	ff.close()

	for feat in data['features']:
	  pro=feat['properties']
	  cou=pro['Cnt_Name']
	  reg=pro['Region']
	  if reg=='': reg=cou
	  ghsl_sum=pro['ghsl_sum']
	  print(cou+','+reg+','+format(ghsl_sum))
	
def calcR0_RKI(cou,reg='',ndays=4,states=[]):
	if cou=='Czech_Republic': cou=cou.replace("_"," ") #!!! PB
	if reg=='' and len(states)==0:
		url = 'https://github.com/ec-jrc/COVID-19/raw/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
		states = pd.read_csv(url, usecols=['Date', 'CountryName', 'CumulativePositive'], parse_dates=['Date'], index_col=['CountryName', 'Date'], squeeze=True).sort_index()  
	if reg!='' and len(states)==0:
		url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
		states = pd.read_csv(url,usecols=['Date', 'CountryName','Region','CumulativePositive'],	parse_dates=['Date'], index_col=['CountryName','Region','Date'],squeeze=True).sort_index()
	x,r0=calcR0_RKIs(cou,states,reg,ndays)
	return x,r0,states

def calcR0_RKIs(state_name,states,reg='',ndays=4):
    if reg=='':
        cases0 = states.xs(state_name).rename(f"{state_name} cases")
    else:
        cases0 = states.xs(state_name).xs(reg).rename(f"{state_name+' '+reg} cases")
    cases = cases0[0:len(cases0)]

    original, smoothed = prepare_cases(cases,10,5)
    sdif=pd.array(smoothed)
    sdif=pd.array(original)
    dd=pd.array(smoothed.index)
    dd=pd.array(original.index)
    
    dates=[];r0=[];rmin=[];rmax=[]
    for k in range(2*ndays,sdif.shape[0]):
        week1=0.0
        week2=0.0
        for j in range(0,ndays):
            week1 += sdif[k-2*ndays+j]
            week2 += sdif[k-ndays+j]
        if week1 !=0:
            rrki=week2/week1
        else:
            rrki=0.0
        if rrki<0:
            rrki=0.0
        dates.append(datetime.datetime.strptime(format(dd[k],'%Y-%m-%d'),'%Y-%m-%d'))
        r0.append(rrki)
    
    return dates,r0

def getAllRt(dirOut):
	ecList = "Austria, Belgium, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg, Malta, Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia, Spain, Sweden".split(', ')
	states=[]
	EU=[]
	for cou in ecList:
		print('calculating R0 for '+cou)
		x,r0,states=calcR0_RKI(cou,'',7,states)
		EU.append([cou,x,r0])

	f=open (dirOut+'allR0.csv','w')
	riga=''
	for eu in EU:
		cou=eu[0]
		riga +='time,'+cou+' r07,,'
	f.write(riga+'\n')
	for k in range(len(EU[0][1])+1):
		riga=''
		for eu in EU:
			try:
				riga +=format(eu[1][k],'%Y-%m-%d')+','+format(eu[2][k].round(2))+',,'
			except:
				riga +=',,,'
		f.write(riga+'\n')
	f.close()

	f=open (dirOut+'allR0_3.csv','w')
	riga='date,'
	for eu in EU:
		cou=eu[0]
		riga +=cou+' r07,'
	f.write(riga+'\n')

	startDate=datetime.datetime(2020,2,1)
	endDate=datetime.datetime.utcnow()
	totDays=(endDate-startDate).days

	for k in range(totDays+1):
		day=startDate+datetime.timedelta(days=k)
		riga=format(day,'%Y-%m-%d')+','
		for eu in EU:
			found=False
			for jj in range(len(eu[1])):
				if eu[1][jj]==day:
					try:
						riga +=format(eu[2][jj].round(2))+','
					except:
						riga +=format(eu[2][jj])+','
					found=True
					break
			if not found:
				riga +=','
		f.write(riga+'\n')
	f.close()

def getlist(fname):
    f=open(fname,'r',encoding="utf-8")
    rows=f.read().split('\n')
    f.close()
    return rows

def getrows(direIn):
	FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS=direIn+'/jrc-covid-19-all-days-by-regions.csv'

	if not os.path.exists(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS):
        #os.remove(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS)
		url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
		cmd='curl -o "'+ FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS+ '" '+url
		print (cmd)
		os.system(cmd)

	rows=getlist(FILE_EPIDEMIOLOGY_WEBCRITECH_REGIONS)
	regions={}
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

	return regions

def getPopRegion(direIn,cou):
    regions=getrows(direIn)

    calibrationPeriod=30
    fname=direIn+ '/FITTING_ALLREGIONS/data/FITTING_ALLREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).json'
    if not os.path.exists(fname):
        fname=direIn+ '/FITTING_SELECTEDREGIONS/data/FITTING_SELECTEDREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).json'
    data=json.load(open(fname))
    for (reg,v) in data.items():
        try:
            couReg=regions[reg][0]
        except:
            couReg=''

        if cou==couReg:
            r0=data[reg]['r0']
            Trecov=data[reg]['Trecov']
            N=data[reg]['population']
            print(cou,',',reg,',',N)

def getVariabilityRegion(direIn,cou,regions=''):
    if regions=='':
        regions=getrows(direIn)

    calibrationPeriod=30
    fname=direIn+ '/FITTING_ALLREGIONS/data/FITTING_ALLREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).json'
    if not os.path.exists(fname):
        fname=direIn+ '/FITTING_SELECTEDREGIONS/data/FITTING_SELECTEDREGIONS_SIR_(-'+format(calibrationPeriod)+' -1).json'
    data=json.load(open(fname))
    values=[[],[]]
    for (reg,v) in data.items():
        try:
            couReg=regions[reg][0]
        except:
            couReg=''

        if cou==couReg:
            r0=data[reg]['r0']
            Trecov=data[reg]['Trecov']
            N=data[reg]['population']
            R=data[reg]['_ivp']['R']
            I=data[reg]['_ivp']['I']
            S=N-I-R
            
            dIdt=r0/Trecov*S*I/N

            SpecIncid=dIdt/N*1e5*14
            CurrRisk=SpecIncid*r0
            values[0].append(SpecIncid)
            values[1].append(CurrRisk)
            #print(cou,',',reg,',',N,',',SpecIncid,',',CurrRisk)

    print(cou,np.mean(values[0]),np.std(values[0]),np.std(values[0])/np.mean(values[0]))
    return regions

def getVariabilityRegions(direIn):
	ecList = "Austria, Belgium, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg, Malta, Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia, Spain, Sweden".split(', ')
	reg=''
	for cou in ecList:
		reg=getVariabilityRegion(direIn,cou,reg)
