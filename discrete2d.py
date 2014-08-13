import sys,os
import numpy as np
import statsmodels.api as sm

if len(sys.argv)<4:
    print 'Usage: pos_peakTag neg_peakTag allPeakTag [binN]'
    print 'The 9th column should be the fragment length'
    exit(1)

interval=[85, 164, 1000000]

posfile=sys.argv[1]
negfile=sys.argv[2]
testfile=sys.argv[3]

allFeat=[[] for i in range(len(interval)+1)]

binN=3
if len(sys.argv)>4:
    binN=int(sys.argv[4])

eps=1e-16

def zeroLog(x):
    if x<eps:
        x=eps
    return np.log(x)

def inv_logit(p):
    t1=np.exp(p)
    return t1/(1+t1)

def genDiscretePr(infile, retTh=[]):
    # pr: [lowB:[tag1, tag2, tag3], medianB:[...], hiB:[...]]
    # threshold: [bRatio:[th1, th2, max], tag:[th1, th2, max]]
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
    
    retPr=[[0]*binN for i in range(binN)]
    num=0
    
    if not retTh:
        # b ratio
        tmp=[]
        for i in range(binN):
            tmp.append(np.percentile(tagDistr[1],(i+1)*100.0/binN))
        retTh.append(tmp)
        
        # tag count
        tmp=[]
        for i in range(binN):
            tmp.append(np.percentile(tagDistr[-1],(i+1)*100.0/binN))
        retTh.append(tmp)
        
        for i in range(len(allFeat)):
            allFeat[i]+=tagDistr[i]
    
    num=len(tagDistr[0])
    for i in range(num):
        for idx0, y in enumerate(retTh[0]):
            if tagDistr[1][i]<=y:
                break
        for idx1, y in enumerate(retTh[1]):
            if tagDistr[-1][i]<=y:
                break
        retPr[idx0][idx1]+=1
    
    for i in range(binN):
        for j in range(binN):
            retPr[i][j]/=1.0*num
        
    return retPr, retTh

def genPostPr(infile, prP, thP, prN, thN):
    nowId=''
    tagOne=[0]*(len(interval)+1)
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
            for idx0, y in enumerate(thP[0]):
                if tagOne[1]<=y:
                    break
            for idx1, y in enumerate(thP[1]):
                if tagOne[-1]<=y:
                    break
            lratio+=zeroLog(prP[idx0][idx1])
            
            for idx0, y in enumerate(thN[0]):
                if tagOne[1]<=y:
                    break
            for idx1, y in enumerate(thN[1]):
                if tagOne[-1]<=y:
                    break
            lratio-=zeroLog(prN[idx0][idx1])
            
            print '\t'.join([nowId, str(inv_logit(lratio))])
            tagOne=[0]*len(tagOne)
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
    for idx0, y in enumerate(thP[0]):
        if tagOne[1]<=y:
            break
    for idx1, y in enumerate(thP[1]):
        if tagOne[-1]<=y:
            break
    lratio+=zeroLog(prP[idx0][idx1])
    
    for idx0, y in enumerate(thN[0]):
        if tagOne[1]<=y:
            break
    for idx1, y in enumerate(thN[1]):
        if tagOne[-1]<=y:
            break
    lratio-=zeroLog(prN[idx0][idx1])
            
    print '\t'.join([nowId, str(inv_logit(lratio))])
    
def printOut(pr, th):
    for x in pr:
        print '\t'.join([str(y) for y in x])
    print '\n'
    for x in th:
        print '\t'.join([str(y) for y in x])
    print '\n'

def getThAll():
    retTh=[]
    num=0
    # b ratio
#     print len(allFeat[1])
    tmp=[]
    for i in range(binN):
        tmp.append(np.percentile(allFeat[1],(i+1)*100.0/binN))
    retTh.append(tmp)
    
    # tag count
    tmp=[]
    for i in range(binN):
        tmp.append(np.percentile(allFeat[-1],(i+1)*100.0/binN))
    retTh.append(tmp)
    
    return retTh

pi=0.3
prP, thP=genDiscretePr(posfile, [])
prN, thN=genDiscretePr(negfile, [])

thAll=getThAll()
prP, thP=genDiscretePr(posfile, thAll)
prN, thN=genDiscretePr(negfile, thAll)

# printOut(prP, thP)
# printOut(prN, thN)
genPostPr(testfile, prP, thP, prN, thN)


