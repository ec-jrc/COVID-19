import sys,os
import datetime
import numpy as np

def readMob(nuts,op, typerec):
	print('Read mobility data')
	fname=r'E:/CV/Modelling_Activity/python/Mobility_Indicators/mobility-indicators-NUTS0123.csv'
	f=open(fname,'r')
	rows=f.read().split('\n')
	f.close()
	#   0         1   2       3  4  5   6                    7
	# 2020-02-01,AT,NUTS0_ID,AT,1,total,2.1828063578875407,normalised
	mobiAll={}
	print('Analysing mobility data')
	for r in rows:
		if r=='' : continue
		try:
#			if not r.split(',')[1]=='IT':
#				continue
			#print(r)
			#if 'ITF,' in r:
				#dd=''
				#print(r)
			p=r.split(',')
			nutsid=p[2]
			
			oper=p[4]
			tip= p[5]
			
			if nutsid==nuts and oper==op and tip==typerec:
				nutscode=p[3]
				if not nutscode in mobiAll:
					#print(nutscode,nutsid)
					mobiAll[nutscode]=[nutscode,[],[]]
				
				if nutscode=='PT30':
					a=1
				if not p[6]=='':
					dd=datetime.datetime.strptime(p[0], '%Y-%m-%d')
					mobVal=float(p[6])
					mobiAll[nutscode][1].append(dd)
					mobiAll[nutscode][2].append(mobVal)
				else:
					print ('** nullo',r)

		except:
			print("Unexpected error:", sys.exc_info()[0])
			raise

	print('End read mobility data')
	return mobiAll

def getRatio(mobiAll,nutscode, date1,period1,period2):
	times=mobiAll[nutscode][1]
	values=mobiAll[nutscode][2]
	i0=len(times)
	for i in range(len(times)):
		if times[i]>date1:
			i0=i
			break
	av1=np.average(values[i0-period1:i0])
	av2=np.average(values[i0-period2:i0-period1])
	if av2!=0:
		ratio=av1/av2
	else:
		ratio=1
	return ratio

def retrieveData(dirEpi,cou):
	for firstDay in range(-260,-30):
		ANALYSIS_HORIZON = (firstDay, firstDay+29)
		lastDay=datetime.datetime.utcnow()+datetime.timedelta(days=ANALYSIS_HORIZON[1])
		#print("Start "+format(ANALYSIS_HORIZON[0]),lastDay)
		DIR_REPORTING = dirEpi+ "/_report/"+lastDay.strftime('%Y-%m-%d')
		if not os.path.isdir(DIR_REPORTING):
			continue
		for STUDY_TYPE in ["FITTING_SELECTEDCOUNTRIES"]:
			csv=DIR_REPORTING+'/'+STUDY_TYPE+'/data/FITTING_SELECTEDCOUNTRIES_SIR_('+format(ANALYSIS_HORIZON[0])+' '+format(ANALYSIS_HORIZON[1])+').csv'
			json=DIR_REPORTING+'/'+STUDY_TYPE+'/data/FITTING_SELECTEDCOUNTRIES_SIR_('+format(ANALYSIS_HORIZON[0])+' '+format(ANALYSIS_HORIZON[1])+').json'
			if os.path.exists(csv):
				f=open(csv,'r')
				rows=f.read().split('\n')
				f.close()
				for r in rows:
					# 'Italy,3.499999999999...1133646668', '''
					if r.split(',')[0]==cou:
						r0=r.split(',')[1]
						trecov=r.split(',')[2]
						print(lastDay,cou,r0,trecov)