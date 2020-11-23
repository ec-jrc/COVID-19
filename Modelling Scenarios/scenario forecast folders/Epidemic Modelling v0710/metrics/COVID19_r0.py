import os, sys, getopt
import pandas as pd
import numpy as np
import datetime
from scipy import stats as sps
from scipy.interpolate import interp1d

# Gamma is 1/serial interval
# https://wwwnc.cdc.gov/eid/article/26/7/20-0282_article
# https://www.nejm.org/doi/full/10.1056/NEJMoa2001316
GAMMA = 1/7


def highest_density_interval(pmf, p=.9, debug=False):
    # If we pass a DataFrame, just call this recursively on the columns
    if(isinstance(pmf, pd.DataFrame)):
        return pd.DataFrame([highest_density_interval(pmf[col], p=p) for col in pmf],
                            index=pmf.columns)
    
    cumsum = np.cumsum(pmf.values)
    
    # N x N matrix of total probability mass for each low, high
    total_p = cumsum - cumsum[:, None]
    
    # Return all indices with total_p > p
    lows, highs = (total_p > p).nonzero()
    
    # Find the smallest range (highest density)
    try:
        best = (highs - lows).argmin()
    except:
        best=highs
        
    low = pmf.index[lows[best]]
    high = pmf.index[highs[best]]
    
    return pd.Series([low, high],
                         index=[f'Low_{p*100:.0f}',
                                f'High_{p*100:.0f}'])
    
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

def get_posteriors(sr, sigma=0.15):
    
    # We create an array for every possible value of Rt
    R_T_MAX = 12
    r_t_range = np.linspace(0, R_T_MAX, R_T_MAX*100+1)

    # (1) Calculate Lambda
    lam = sr[:-1].values * np.exp(GAMMA * (r_t_range[:, None] - 1))

    
    # (2) Calculate each day's likelihood
    likelihoods = pd.DataFrame(
        data = sps.poisson.pmf(sr[1:].values, lam),
        index = r_t_range,
        columns = sr.index[1:])
    
    # (3) Create the Gaussian Matrix
    process_matrix = sps.norm(loc=r_t_range,
                              scale=sigma
                             ).pdf(r_t_range[:, None]) 

    # (3a) Normalize all rows to sum to 1
    process_matrix /= process_matrix.sum(axis=0)
    
    # (4) Calculate the initial prior
    #prior0 = sps.gamma(a=4).pdf(r_t_range)
    prior0 = np.ones_like(r_t_range)/len(r_t_range)
    prior0 /= prior0.sum()

    # Create a DataFrame that will hold our posteriors for each day
    # Insert our prior as the first posterior.
    #try:
    posteriors = pd.DataFrame(
            index=r_t_range,
            columns=sr.index,
            data={sr.index[0]: prior0}
        )
    #except:
    #    print ('error in data')
    # We said we'd keep track of the sum of the log of the probability
    # of the data for maximum likelihood calculation.
    log_likelihood = 0.0

    # (5) Iteratively apply Bayes' rule
    for previous_day, current_day in zip(sr.index[:-1], sr.index[1:]):

        #(5a) Calculate the new prior
        current_prior = process_matrix @ posteriors[previous_day]
        
        #(5b) Calculate the numerator of Bayes' Rule: P(k|R_t)P(R_t)
        numerator = likelihoods[current_day] * current_prior
        
        #(5c) Calcluate the denominator of Bayes' Rule P(k)
        denominator = np.sum(numerator)
        
        # Execute full Bayes' Rule
        posteriors[current_day] = numerator/denominator
        
        # Add to the running sum of log likelihoods
        log_likelihood += np.log(denominator)
    
    return posteriors, log_likelihood

def calcR0s(state_name,states,reg=''):
    # if states=='':
    #     url = 'https://github.com/ec-jrc/COVID-19/raw/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
    #     states = pd.read_csv(url,
    #                          usecols=['Date', 'CountryName', 'CumulativePositive'],
    #                          parse_dates=['Date'],
    #                          index_col=['CountryName', 'Date'],
    #                      squeeze=True).sort_index()
    
#    state_name = 'Belgium'
    if reg=='':
        cases0 = states.xs(state_name).rename(f"{state_name} cases")
    else:
        cases0 = states.xs(state_name).xs(reg).rename(f"{state_name+' '+reg} cases")
    cases = cases0[0:len(cases0)]
    #print (cases)
    original, smoothed = prepare_cases(cases)
    

    
    # Note that we're fixing sigma to a value just for the example
    posteriors, log_likelihood = get_posteriors(smoothed, sigma=.25)
    
    # Note that this takes a while to execute - it's not the most efficient algorithm
    hdis = highest_density_interval(posteriors, p=.9)
    
    most_likely = posteriors.idxmax().rename('ML')
    
    # Look into why you shift -1
    result = pd.concat([most_likely, hdis], axis=1)
    csv=result.to_csv()
    rows=csv.replace('\r','').split('\n')
    x=[];r0=[];rl=[];rh=[]
    for r in rows:
        
        if r.startswith('Date'):continue
        if r=='':continue
        p=r.split(',')
        d=datetime.datetime.strptime(p[0],'%Y-%m-%d')
        
        try:
            r0f=float(p[1]);rlf=float(p[2]);rhf=float(p[3])
            r0.append(r0f)
            rl.append(rlf)
            rh.append(rhf)
            x.append(d)
        except:
            #print ('error in data for ',d,r)
            r=r
    
    #print ('AA',np.shape(x),np.shape(r0),p)
    return x,r0,rl,rh, csv
    #result.to_csv(state_name+'rt.csv')

def calcR0(cou):
	if cou=='Czech_Republic': cou=cou.replace("_"," ") #!!! PB
	url = 'https://github.com/ec-jrc/COVID-19/raw/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
	states = pd.read_csv(url, usecols=['Date', 'CountryName', 'CumulativePositive'], parse_dates=['Date'], index_col=['CountryName', 'Date'], squeeze=True).sort_index()  
	x,r0,rl,rh, csv=calcR0s(cou,states)
	return x,r0,rl,rh,csv

def gd(a1,a0):
    return (a1-a0).days+(a1-a0).seconds/3600./24.

def calcR0_JRC(cou,reg=''):
	if cou=='Czech_Republic': cou=cou.replace("_"," ") #!!! PB
	if reg=='':
		url = 'https://github.com/ec-jrc/COVID-19/raw/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
		states = pd.read_csv(url, usecols=['Date', 'CountryName', 'CumulativePositive'], parse_dates=['Date'], index_col=['CountryName', 'Date'], squeeze=True).sort_index()  
	else:
		url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
		states = pd.read_csv(url,usecols=['Date', 'CountryName','Region','CumulativePositive'],	parse_dates=['Date'], index_col=['CountryName','Region','Date'],squeeze=True).sort_index()
	x,r0,rl,rh=calcR0_JRCs(cou,states,reg)
	return x,r0,rl,rh  

def calcR0_7d(cou,reg='',ndays=7):
	if cou=='Czech_Republic': cou=cou.replace("_"," ") #!!! PB
	if reg=='':
		url = 'https://github.com/ec-jrc/COVID-19/raw/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
		states = pd.read_csv(url, usecols=['Date', 'CountryName', 'CumulativePositive'], parse_dates=['Date'], index_col=['CountryName', 'Date'], squeeze=True).sort_index()  
	else:
		url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
		states = pd.read_csv(url,usecols=['Date', 'CountryName','Region','CumulativePositive'],	parse_dates=['Date'], index_col=['CountryName','Region','Date'],squeeze=True).sort_index()
	x,r0=calcR0_7ds(cou,states,reg,ndays)
	return x,r0

 
def calcR0_RKI(cou,reg='',ndays=4):
	if cou=='Czech_Republic': cou=cou.replace("_"," ") #!!! PB
	if reg=='':
		url = 'https://github.com/ec-jrc/COVID-19/raw/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
		states = pd.read_csv(url, usecols=['Date', 'CountryName', 'CumulativePositive'], parse_dates=['Date'], index_col=['CountryName', 'Date'], squeeze=True).sort_index()  
	else:
		url='https://raw.githubusercontent.com/ec-jrc/COVID-19/master/data-by-region/jrc-covid-19-all-days-by-regions.csv'
		states = pd.read_csv(url,usecols=['Date', 'CountryName','Region','CumulativePositive'],	parse_dates=['Date'], index_col=['CountryName','Region','Date'],squeeze=True).sort_index()
	x,r0=calcR0_RKIs(cou,states,reg,ndays)
	return x,r0

def calcR0_RKIs(state_name,states,reg='',ndays=4):
    if reg=='':
        cases0 = states.xs(state_name).rename(f"{state_name} cases")
    else:
        cases0 = states.xs(state_name).xs(reg).rename(f"{state_name+' '+reg} cases")
    cases = cases0[0:len(cases0)]

    original, smoothed = prepare_cases(cases,10,5)
    sdif=pd.array(smoothed)
    dd=pd.array(smoothed.index)
    
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
            
        dates.append(datetime.datetime.strptime(format(dd[k],'%Y-%m-%d'),'%Y-%m-%d'))
        r0.append(rrki)
    
    return dates,r0

def calcR0_7ds(state_name,states,reg='',ndays=7):
    if reg=='':
        cases0 = states.xs(state_name).rename(f"{state_name} cases")
    else:
        cases0 = states.xs(state_name).xs(reg).rename(f"{state_name+' '+reg} cases")
    cases = cases0[0:len(cases0)]

    original, smoothed = prepare_cases(cases,10,5)
    sdif=pd.array(smoothed)
    dd=pd.array(smoothed.index)
    
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
            
        dates.append(datetime.datetime.strptime(format(dd[k],'%Y-%m-%d'),'%Y-%m-%d'))
        r0.append(rrki)
    
    return dates,r0

def calcR0_JRCs(state_name,states,reg=''):
    if reg=='':
        cases0 = states.xs(state_name).rename(f"{state_name} cases")
    else:
        cases0 = states.xs(state_name).xs(reg).rename(f"{state_name+' '+reg} cases")
    cases = cases0[0:len(cases0)]

    original, smoothed = prepare_cases(cases,10,5)
    sdif=pd.array(smoothed)
    dd=pd.array(smoothed.index)
    genTime=7.0
    dates=[];r0=[];rmin=[];rmax=[]
    for k in range(8,sdif.shape[0]):
        if sdif[k]>0 and sdif[k-7]>0:
            vr0=(np.log(sdif[k])-np.log(sdif[k-7]))*genTime/gd(dd[k],dd[k-7])+1
            vrmax=(np.log(sdif[k])-np.log(sdif[k-7]))*(genTime+1)/gd(dd[k],dd[k-7])+1
            vrmin=(np.log(sdif[k])-np.log(sdif[k-7]))*(genTime-1)/gd(dd[k],dd[k-7])+1
            #print(r0,sdif[k],sdif[k-7],gd(dd[k],dd[k-7]),np.log(sdif[k]),np.log(sdif[k-7]))
        else:
            vr0=0.0
            vrmin=0.0
            vrmax=0.0
        dates.append( datetime.datetime.strptime(format(dd[k],'%Y-%m-%d'),'%Y-%m-%d'))
        r0.append(vr0)
        rmin.append(vrmin)
        rmax.append(vrmax)
    
    return dates,r0,rmin,rmax

def calcR0_CRAN(cou,readFolder,reg='',reg0=''):
    print(' calcR0_cran' ,cou,reg,readFolder)
    if not os.path.exists(readFolder):
        os.makedirs(readFolder)
    create=False
    sreg=''
    if reg!='':
        sreg=' '+reg
    print('checking '+readFolder+'\\'+cou+sreg+' _R0.csv')
    if not (os.path.exists(readFolder+'\\'+cou+sreg+' _R0.csv') and os.path.exists(readFolder+'\\'+cou+sreg+' _R0_COU.csv') and os.path.exists(readFolder+'\\'+cou+sreg+' _R0_confint.csv')):
        print('file not existing, to be created '+readFolder+'\\'+cou+reg+' _R0.csv')
        create=True
    else:
        st=os.stat(readFolder+'\\'+cou+sreg+' _R0.csv')
        fdate=datetime.datetime.fromtimestamp(st.st_mtime)
        age=datetime.datetime.now()-fdate
        #print(age.days+age.seconds/(3600.0*24.0))
        if age.days+age.seconds/(3600.0*24.0)>0.5:
            print('file existing but old, to be created')
            create=True
    print('>'+cou+'<>'+reg+'<')
    if create:
        if (reg==''):
            cmd=r'"C:\Program Files\R\R-4.0.0\bin\Rscript.exe"'+' calcR0.R -mode NAT -o '+readFolder+' -c '+cou.replace(" ","_")
        else:
            cmd=r'"C:\Program Files\R\R-4.0.0\bin\Rscript.exe"'+' calcR0.R -mode REG -o '+readFolder+' -c '+cou.replace(" ","_")   +' -r '+reg.replace(" ","_")         
        print('calcR0_CRAN command= ' + cmd)
        os.system(cmd)
        
    dates=[];r0=[];rlow=[];rhigh=[]
    fname=readFolder+'\\'+cou+sreg+' _R0.csv'
    if not os.path.exists(fname):
            reg=reg0
    fname=readFolder+'\\'+cou+sreg+' _R0.csv'            
    if not os.path.exists(fname):
        return  dates,r0,rlow,rhigh 

    fname=readFolder+'\\'+cou+sreg+' _R0.csv';f=open(fname,'r');textR0=f.read();f.close()
    fname=readFolder+'\\'+cou+sreg+' _R0_COU.csv';f=open(fname,'r',encoding='UTF-8');textCOU=f.read();f.close()
    fname=readFolder+'\\'+cou+sreg+' _R0_confint.csv';f=open(fname,'r');textCI=f.read();f.close()
    lines=textCOU.split('\n')
    for k in range(2,len(lines)-1):
        p=lines[k].split(',')
        #print(p)
        #print(p[len(p)-2])
        d=datetime.datetime.strptime(p[len(p)-2],'%Y-%m-%d')
        dates.append(d)

    lines1=textR0.split('\n')
    lines2=textCI.split('\n')
    r0.append(0);rlow.append(0);rhigh.append(0)
    #r0.append(0);rlow.append(0);rhigh.append(0)
    #r0.append(0);rlow.append(0);rhigh.append(0)
    
    for k in range(2,len(lines1)-2):
        if lines1[k]=='':continue
        p=lines2[k].split(',')
        r0.append(float(lines1[k]))
        rlow.append(float(p[0]))
        rhigh.append(float(p[1]))
    
   # print (len(r0),len(dates))
   # for k in range(len(dates)):
   #      print(dates[k],r0[k])
    return dates,r0,rlow,rhigh

    
if __name__ == "__main__":


    opts, args = getopt.getopt(sys.argv[1:],"c:o:i:h",["cou=","odir=","idir=","help="])

    if len(opts)<1:
        print ('INPUT ERROR')
        print ('COVID19_r0.py  -c <country|ALL> -o <outputdir>  ')
        print ('  example:')
        print ('python D:\mnt\output\COVID-19\scripts_py3\COVID19_r0.py -c Italy -o D:\mnt\output\COVID-19\pdfOut\ ')
        sys.exit()
        

    for opt, arg in opts:
        
        if opt in ("-h", "--help"):
            print ('python D:\mnt\output\COVID-19\scripts\COVID19_r0.py -c Italy -o D:\mnt\output\COVID-19\pdfOut\ ')
            sys.exit()
        elif opt in ("-o", "--odir"):
            dirout = arg			
        elif opt in ("-i", "--idir"):
            readFolder = arg			
        elif opt in ("-c", "--cou"):
            cou = arg
        else:
            print ('error parameter not available: '+opt)
            sys.exit()



    url = 'https://github.com/ec-jrc/COVID-19/raw/master/data-by-country/jrc-covid-19-all-days-by-country.csv'
    states = pd.read_csv(url,usecols=['Date', 'CountryName', 'CumulativePositive'],parse_dates=['Date'],index_col=['CountryName', 'Date'],squeeze=True).sort_index()  
    
    now=datetime.datetime.now().strftime("%Y%m%d")        
    rowLast="country,date,r0_ML,Low_90',High_90\n"
    print (cou)
    if cou=="ALL":
        ecList = "Austria, Belgium, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Latvia, Lithuania, Luxembourg, Malta, Netherlands, Poland, Portugal, Romania, Slovakia, Slovenia, Spain, Sweden".split(', ')
        eucmParticip = "Iceland,Montenegro,North Macedonia,Norway,Turkey,United Kingdom,Switzerland".split(',')
        for cou in ecList:
            print(cou)
            try:
                x,r0,rl,rh,csv=calcR0s(cou,states)
                f=open(dirout+"\\"+now+"_"+cou+"_r0.csv","w")
                f.write(csv.replace("\n",""))
                f.close()
                rowLast +=cou+','+format(x[len(x)-1])+','+"%.2f" % r0[len(r0)-1]+','+"%.2f" % rl[len(rl)-1]+','+"%.2f" % rh[len(rh)-1]+'\n'
            except:
                rowLast +=cou+",,,,\n"

        for cou in eucmParticip:
            print(cou)
            try:
                x,r0,rl,rh,csv=calcR0s(cou,states)
                f=open(dirout+"\\"+now+"_"+cou+"_r0.csv","w")
                f.write(csv.replace("\n",""))
                f.close()
                rowLast +=cou+','+format(x[len(x)-1])+','+"%.2f" % r0[len(r0)-1]+','+"%.2f" % rl[len(rl)-1]+','+"%.2f" % rh[len(rh)-1]+'\n'
            except:
                rowLast +=cou+",,,,\n"

    else:
        #try:
            x,r0,rl,rh,csv=calcR0s(cou,states)
            #print(r0)
            f=open(dirout+"\\"+now+"_"+cou+"_r0.csv","w")
            f.write(csv.replace("\n",""))
            f.close()
            rowLast +=cou+','+format(x[len(x)-1])+','+"%.2f" % r0[len(r0)-1]+','+"%.2f" % rl[len(rl)-1]+','+"%.2f" % rh[len(rh)-1]+'\n'
            if readFolder !='':
                x,r0,rl,rh=calcR0_CRAN(cou,cou,readFolder)
            
       # except:
       #    print('Error calling calcR0')
       #    rowLast +=cou+",,,,"
    
    f=open(dirout+"\\"+now+"_r0_last.csv","w")
    f.write(rowLast)
    f.close()
# #x,r0,rl,rh=calcR0('Italy')
# x,r0,rl,rh=calcR0('Austria')
# print(r0)