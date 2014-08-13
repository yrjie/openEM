import sys,os
import numpy as np
import statsmodels.api as sm

if len(sys.argv)<4:
    print 'Usage: pos_peakTag neg_peakTag allPeakTag [binN]'
    print 'The 9th column should be the fragment length'
    exit(1)

posfile=sys.argv[1]
negfile=sys.argv[2]
testfile=sys.argv[3]

binN=3
if len(sys.argv)>4:
    binN=int(sys.argv[4])

if binN==3:
    interval=[85, 164, 1000000]
else:
    interval=[100, 180, 247, 315, 473, 558, 615, 1000000]
    
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

def genMultiNPr(infile):
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
#             print tagOne
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
    retPr=[]
    for i in range(len(interval)):
        retPr.append(np.mean(tagDistr[i]))
    retPr.append(sm.distributions.ECDF(tagDistr[-1]))
#         print np.mean(x)
#         retEcdf.append(sm.distributions.ECDF(x))
#     for i in range(len(tagDistr)):
#         x=tagDistr[i]
#         y=retEcdf[i]
#         temp=y([min(x)+i*(max(x)-min(x))/10 for i in range(11)])-y([min(x)+i*(max(x)-min(x))/10 for i in range(-1,10)])
#         print '\t'.join([str(m) for m in temp])
    return retPr

def genPostPr(infile, prP, prN):
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
#             if tagOne[-1]>0:
#                 for i in range(len(interval)):
#                     tagOne[i]/=1.0*tagOne[-1]
            for i in range(len(interval)):
                lratio+=tagOne[i]*(zeroLog(prP[i])-zeroLog(prN[i]))
            x=tagOne[-1]
#             lratio+=zeroLog(prP[-1](x)-prP[-1](x-10))
#             lratio-=zeroLog(prN[-1](x)-prN[-1](x-10))
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
#     if tagOne[-1]>0:
#         for i in range(len(interval)):
#             tagOne[i]/=1.0*tagOne[-1]
    for i in range(len(interval)):
        lratio+=tagOne[i]*(zeroLog(prP[i])-zeroLog(prN[i]))
    x=tagOne[-1]
#     lratio+=zeroLog(prP[-1](x)-prP[-1](x-10))
#     lratio-=zeroLog(prN[-1](x)-prN[-1](x-10))
    print '\t'.join([nowId, str(inv_logit(lratio))])
    

pi=0.3
prP=genMultiNPr(posfile)
prN=genMultiNPr(negfile)
print prP
print prN
# genPostPr(testfile, prP, prN)
# pi=1.0*num1/(num1+num2)
# print pi

