import sys,os
import numpy as np
import statsmodels.api as sm
import scipy.optimize as op

if len(sys.argv)<4:
    print 'Usage: pos_peakTag neg_peakTag allPeakTag'
    print 'The 9th column should be the fragment length'
    exit(1)

posfile=sys.argv[1]
negfile=sys.argv[2]
testfile=sys.argv[3]

# interval=[100, 180, 247, 315, 473, 558, 615, 1000000]
interval=[85, 164, 1000000]
eps=1e-16
pi=0.5

def zeroTest(x):
    if x<eps:
        x=eps
    return x

def zeroLog(x):
    if x<eps:
        x=eps
    return np.log(x)

def inv_logit(p):
    t1=np.exp(p)
    return t1/(1+t1)

def q_sum(x, data, Z, lam):
    mu1=x[0]
    sig1=x[1]
    mu0=x[2]
    sig0=x[3]
    a=-sum(Z*scipy.log(ss.norm.pdf(data, loc=mu1, scale=sig1)*lam))
    b=-sum((1-Z)*scipy.log(ss.norm.pdf(data, loc=mu0, scale=sig0)*(1-lam)))
    ret=a+b
    #ret=-sum(Z*scipy.log(ss.norm.pdf(data, loc=mu, scale=sig)*lam))
    #print sum1,sum0
    return ret

def E_step(x1,x0,lam,data):
    a1=scipy.log(ss.norm.pdf(data, loc=x1[0], scale=x1[1])*lam)
    a0=scipy.log(ss.norm.pdf(data, loc=x0[0], scale=x0[1])*(1-lam))
    lratio=a1-a0
    return inv_logit(lratio)

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
    	    if fragLen<=x:
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
    lratio=np.log(pi)-np.log(1-pi)
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
                lratio+=zeroLog(ecdfP[i](x)-ecdfP[i](x-0.1))
                lratio-=zeroLog(ecdfN[i](x)-ecdfN[i](x-0.1))
            x=tagOne[-1]
            a*=ecdfP[-1](x)-ecdfP[-1](x-10)
            b*=ecdfN[-1](x)-ecdfN[-1](x-10)
            lratio+=zeroLog(ecdfP[-1](x)-ecdfP[-1](x-10))
            lratio-=zeroLog(ecdfN[-1](x)-ecdfN[-1](x-10))
#             print ecdfP[-1](x)-ecdfP[-1](x-1)
#             print ecdfN[-1](x)-ecdfN[-1](x-1)
#             if a+b<eps:
#                 b=eps
#             print '\t'.join([nowId,str(a/(a+b))])
            print '\t'.join([nowId, str(inv_logit(lratio))])
            tagOne=[0]*len(tagOne)
            a=pi
            b=1-pi
            lratio=np.log(pi)-np.log(1-pi)
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
        lratio+=zeroLog(ecdfP[i](x)-ecdfP[i](x-0.1))
        lratio-=zeroLog(ecdfN[i](x)-ecdfN[i](x-0.1))
    x=tagOne[-1]
    a*=ecdfP[-1](x)-ecdfP[-1](x-10)
    b*=ecdfN[-1](x)-ecdfN[-1](x-10)
    lratio+=zeroLog(ecdfP[-1](x)-ecdfP[-1](x-10))
    lratio-=zeroLog(ecdfN[-1](x)-ecdfN[-1](x-10))
#     if a+b<eps:
#         b=eps
#     print '\t'.join([nowId,str(a/(a+b))])
    print '\t'.join([nowId, str(inv_logit(lratio))])
    





def genEcdf2(infile):
    nowId=''
    tagDistr=[[] for i in range(2)]
    tagNum1=0
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
            tagDistr[1].append(tagNum1)
            tagNum1=0
            nowId=id
        tagDistr[0].append(fragLen)
        tagNum1+=1
    fi.close()
    tagDistr[1].append(tagNum1)
    num=len(tagDistr[1])
#     print num
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

def genPostPr2(infile, ecdfP, ecdfN):
    nowId=''
    tagNum1=0
    a=pi
    b=1-pi
    lratio=np.log(pi)-np.log(1-pi)
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
#             print '\t'.join([nowId,str(a/(a+b))])
            lratio+=zeroLog(ecdfP[1](tagNum1)-ecdfP[1](tagNum1-10))
            lratio-=zeroLog(ecdfN[1](tagNum1)-ecdfN[1](tagNum1-10))
            print '\t'.join([nowId, str(inv_logit(lratio))])
            tagNum1=0
            lratio=np.log(pi)-np.log(1-pi)
            nowId=id
#         lratio+=zeroLog(ecdfP[0](fragLen)-ecdfP[0](fragLen-1))
#         lratio-=zeroLog(ecdfN[0](fragLen)-ecdfN[0](fragLen-1))
        tagNum1+=1
    fi.close()
    lratio+=zeroLog(ecdfP[1](tagNum1)-ecdfP[1](tagNum1-10))
    lratio-=zeroLog(ecdfN[1](tagNum1)-ecdfN[1](tagNum1-10))
#     print '\t'.join([nowId,str(a/(a+b))])
    print '\t'.join([nowId, str(inv_logit(lratio))])


pi=0.3
# num1, ecdfP=genEcdf(posfile)
# num2, ecdfN=genEcdf(negfile)
# postPr=genPostPr(testfile, ecdfP, ecdfN)

num1, ecdfP=genEcdf2(posfile)
num2, ecdfN=genEcdf2(negfile)
postPr=genPostPr2(testfile, ecdfP, ecdfN)
