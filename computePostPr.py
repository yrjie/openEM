import sys,os
import numpy as np
import statsmodels.api as sm

if len(sys.argv)<4:
    print 'Usage: pos_peakTag neg_peakTag allPeakTag'
    exit(1)

posfile=sys.argv[1]
negfile=sys.argv[2]
testfile=sys.argv[3]

interval=[100, 180, 247, 315, 473, 558, 615, 1000000]
eps=1e-8
pi=0.5

def genEcdf(infile):
    nowId=''
    tagDistr=[[] for i in range(len(interval)+1)]
    tagOne=[0]*len(tagDistr)
    fi=open(infile)
    for line in fi:
    	line=line.strip()
    	if len(line)<1:
    	    continue
    	temp=line.split('\t')
    	fragLen=int(temp[8])
    	if fragLen<0:
    	    fragLen=-fragLen
    	id='\t'.join([temp[0],temp[1],temp[2]])
        if nowId=='':
            nowId=id
    	if id != nowId:
            if tagOne[-1]>0:
        	    for i in range(len(interval)):
        	    	tagOne[i]/=1.0*tagOne[-1]
    	    for i in range(len(tagOne)):
    	    	tagDistr[i].append(tagOne[i])
    	    tagOne=[0]*len(tagOne)
    	    nowId=id
        for i,x in enumerate(interval):
    	    if fragLen<x:
    	    	tagOne[i]+=1
    	    	break
        tagOne[-1]+=1
    fi.close()
    for i in range(len(interval)):
    	tagOne[i]/=1.0*tagOne[-1]
    for i in range(len(tagOne)):
    	tagDistr[i].append(tagOne[i])
    num=len(tagDistr[0])
    retEcdf=[]
    for x in tagDistr:
#         print np.mean(x)
    	retEcdf.append(sm.distributions.ECDF(x))
#     for i in range(len(tagDistr)):
#         x=tagDistr[i]
#         y=retEcdf[i]
#         temp=y([min(x)+i*(max(x)-min(x))/10 for i in range(11)])-y([min(x)+i*(max(x)-min(x))/10 for i in range(-1,10)])
#         print '\t'.join([str(m) for m in temp])
    return num, retEcdf

def genPostPr(infile, ecdfP, ecdfN):
    nowId=''
    tagOne=[0]*(len(interval)+1)
    a=pi
    b=1-pi
    fi=open(infile)
    for line in fi:
    	line=line.strip()
    	if len(line)<1:
    	    continue
    	temp=line.split('\t')
    	fragLen=int(temp[8])
    	if fragLen<0:
    	    fragLen=-fragLen
        id='\t'.join([temp[0],temp[1],temp[2]])
        if nowId=='':
            nowId=id
        if id != nowId and nowId != '':
            if tagOne[-1]>0:
                for i in range(len(interval)):
                    tagOne[i]/=1.0*tagOne[-1]
            for i in range(len(interval)):
                x=tagOne[i]
                a*=ecdfP[i](x)-ecdfP[i](x-0.1)
                b*=ecdfN[i](x)-ecdfN[i](x-0.1)
            x=tagOne[-1]
            a*=ecdfP[-1](x)-ecdfP[-1](x-10)
            b*=ecdfN[-1](x)-ecdfN[-1](x-10)
#             print ecdfP[-1](x)-ecdfP[-1](x-1)
#             print ecdfN[-1](x)-ecdfN[-1](x-1)
            if a+b<eps:
                b=eps
            print '\t'.join([nowId,str(a/(a+b))])
            tagOne=[0]*len(tagOne)
            a=pi
            b=1-pi
            nowId=id
        for i,x in enumerate(interval):
            if fragLen<x:
                tagOne[i]+=1
                break
        tagOne[-1]+=1
    fi.close()
    if tagOne[-1]>0:
        for i in range(len(interval)):
            tagOne[i]/=1.0*tagOne[-1]
    for i in range(len(interval)):
        x=tagOne[i]
        a*=ecdfP[i](x)-ecdfP[i](x-0.1)
        b*=ecdfN[i](x)-ecdfN[i](x-0.1)
    x=tagOne[-1]
    a*=ecdfP[-1](x)-ecdfP[-1](x-10)
    b*=ecdfN[-1](x)-ecdfN[-1](x-10)
    if a+b<eps:
        b=eps
    print '\t'.join([nowId,str(a/(a+b))])
    
	

num1, ecdfP=genEcdf(posfile)
num2, ecdfN=genEcdf(negfile)
pi=0.3
# pi=1.0*num1/(num1+num2)
# print pi
postPr=genPostPr(testfile, ecdfP, ecdfN)
