import sys,os
import numpy as np
import statsmodels.api as sm
import operator

if len(sys.argv)<3:
    print 'Usage: allPeakTag binN [debug]'
    print 'The 9th column should be the fragment length'
    exit(1)

minB=50
maxB=200

interval=[85, 164, 1000000]

infile=sys.argv[1]
allPeak=[]


class Peak:
    id=''
    lenNum=[]
    tagNum=0
    bRatio=0
    sig=0
    z=0.5
    def __init__(self, id, lenNum, tagNum, sig):
        self.id=id
        self.lenNum=lenNum
        self.tagNum=tagNum
        self.sig=sig

binN=int(sys.argv[2])

debug=0
if len(sys.argv)>3:
    debug=int(sys.argv[3])

eps=1e-16
epsEM=1e-3
maxEM=1000

def zeroLog(x):
    if x<eps:
        x=eps
    return np.log(x)

def inv_logit(p):
    t1=np.exp(p)
    return t1/(1+t1)

def readFile(infile, allPeak):
    # pr: [lowB:[tag1, tag2, tag3], medianB:[...], hiB:[...]]
    # threshold: [bRatio:[th1, th2, max], tag:[th1, th2, max]]
    nowId=''
    tagDistr=[[] for i in range(maxB-minB+1)]
    tagOne=[0]*len(tagDistr)
    tagNum=0
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
        sig=float(temp[3])
        if nowId=='':
            nowId=id
        if id != nowId:
#             print tagOne
#             if tagNum>0:
#                 for i in range(len(tagOne)):
#                     tagOne[i]/=1.0*tagNum
#                 tagOne=list(np.cumsum(tagOne))
            for i in range(len(tagOne)):
                tagDistr[i].append(tagOne[i])
            allPeak.append(Peak(nowId, tagOne, tagNum, sig))
            nowId=id
            sig=float(temp[3])
            tagOne=[0]*len(tagOne)
            tagNum=0
        if fragLen>=minB and fragLen<=maxB:
            tagOne[fragLen-minB]+=1
        tagNum+=1
    fi.close()
    if tagNum>0:
        for i in range(len(tagOne)):
            tagOne[i]/=1.0*tagNum
    for i in range(len(tagOne)):
        tagDistr[i].append(tagOne[i])
    allPeak.append(Peak(nowId, tagOne, tagNum, sig))
    if debug:
        print len(allPeak)

def genDistr(peakSet, th, intA, intB):
    num=len(peakSet)
    allBratio=[x.bRatio for x in peakSet]
    allTagnum=[x.tagNum for x in peakSet]
    retPr=[[0]*binN for i in range(binN)]
    for i in range(num):
        for idx0, y in enumerate(th[0]):
            if allBratio[i]<=y:
                break
        for idx1, y in enumerate(th[1]):
            if allTagnum[i]<=y:
                break
        retPr[idx0][idx1]+=1
    
    for i in range(binN):
        for j in range(binN):
            retPr[i][j]/=1.0*num
    
    return retPr

def E_step(peakSet, intA, intB, edistrP, edistrN, th, pi):
    num=len(peakSet)
    allBratio=[x.bRatio for x in peakSet]
    allTagnum=[x.tagNum for x in peakSet]
    for i in range(num):
        for idx0, y in enumerate(th[0]):
            if allBratio[i]<=y:
                break
        for idx1, y in enumerate(th[1]):
            if allTagnum[i]<=y:
                break
        lratio=0
        lratio+=zeroLog(edistrP[idx0][idx1])+zeroLog(pi)
        lratio-=zeroLog(edistrN[idx0][idx1])+zeroLog(1-pi)
        peakSet[i].z=inv_logit(lratio)
#         print pi, edistrP[idx0][idx1], edistrN[idx0][idx1], lratio
#         print peakSet[i].z

def M_step(allPeak):
    allZ=[x.z for x in allPeak]
    pi=sum(allZ)/len(allZ)
#     print "maxZ: ", max(allZ)
#     loZ=np.percentile(allZ, 40*pi)
#     loZ=np.percentile(allZ, 100*(1-pi))
#     hiZ=np.percentile(allZ, 100-80*pi)
#     hiZ=0.5
#     loZ=0.5
#     print len(allPeak), hiZ, loZ
#     posPk=[x for x in allPeak if x.z>hiZ-eps]
#     negPk=[x for x in allPeak if x.z<loZ+eps]
#     smpSize=10000
#     if len(posPk)>smpSize*pi:
#         posPk=posPk[0:int(smpSize*pi)]
#     if len(negPk)>smpSize*(1-pi):
#         negPk=negPk[0:int(smpSize*(1-pi))]
    
    num=len(allPeak)
    allPeak.sort(key=operator.attrgetter('z'))
    posPk=[allPeak[i] for i in range(int((1-pi)*num), num)]
    negPk=[allPeak[i] for i in range(0,int((1-pi)*num))]
    
    if debug:
        print len(posPk), len(negPk)
    posTag=sum([x.tagNum for x in posPk])
    negTag=sum([x.tagNum for x in negPk])
    posLen=[]
    negLen=[]
    diffB=0   
    intA=minB
    intB=maxB
    
    diffNP=0
    intA1=minB
    intB1=maxB
    th=[[0]*binN for i in range(2)]
    for i in range(minB, maxB):
        posLen.append(1.0*sum([x.lenNum[i-minB] for x in posPk])/posTag)
        negLen.append(1.0*sum([x.lenNum[i-minB] for x in negPk])/negTag)
    tmpDiff=0
    j=0
    for i in range(len(posLen)):
        t1=posLen[i]-negLen[i]
        if t1<0:
            tmpDiff=0
            j=i
        else:
            tmpDiff+=t1
        if diffB<tmpDiff:
            diffB=tmpDiff
            intA=j+minB
            intB=i+minB
#     for i in range(minB, maxB):
#         for j in range(i+1, maxB):
#             posB=sum([x.lenCum[j-minB]-x.lenCum[i-minB] for x in posPk])
#             negB=sum([x.lenCum[j-minB]-x.lenCum[i-minB] for x in negPk])
#             tmpDiff=(1.0*posB/posTag-1.0*negB/negTag)/(j-i)
#             if tmpDiff>diffB:
#                 diffB=tmpDiff
#                 intA=i
#                 intB=j
#             if tmpDiff<diffNP:
#                 diffNP=tmpDiff
#                 intA1=i
#                 intB1=j
    if debug:
        print "diffPN", intA, intB, diffB
#         print "diffNP", intA1, intB1, diffNP
#     intA=85
#     intB=164
    for x in allPeak:
        x.bRatio=1.0*sum(x.lenNum[intA-minB+1:intB-minB+1])/x.tagNum
    allBratio=[x.bRatio for x in allPeak].sort()
    for i in range(binN):
        th[0][i]=allBratio[(i+1)*num/binN]
    
    allTagnum=[x.tagNum for x in allPeak].sort()
    for i in range(binN):
        th[1][i]=allTagnum[(i+1)*num/binN]
        
    edistrP=genDistr(posPk, th, intA, intB)
    edistrN=genDistr(negPk, th, intA, intB)
    
    return intA, intB, edistrP, edistrN, th, pi

def runEM():
    oldA=intA=minB
    oldB=intB=maxB
#     allTagNum=[x.tagNum for x in allPeak]
#     loTag=np.percentile(allTagNum,10)
#     hiTag=np.percentile(allTagNum,90)
#     
#     allSig=[x.sig for x in allPeak]
#     loSig=np.percentile(allSig,10)+eps
#     hiSig=np.percentile(allSig,90)-eps
#     tmpPk=[x for x in allPeak if x.tagNum<=loTag]
#     tmpPk=[x for x in allPeak if x.sig<=loSig]
#     for x in tmpPk:
#         x.z=0
#     tmpPk=[x for x in allPeak if x.tagNum>=hiTag]
#     tmpPk=[x for x in allPeak if x.sig>=hiSig]
#     for x in tmpPk:
#         x.z=1
    allPeak.sort(key=operator.attrgetter('sig'))
    num=len(allPeak)
    for i in range(0,num/10):
        allPeak.z=0
    
    for i in range(num*9/10, num):
        allPeak.z=1
    
    oldPi=pi=0.5
    oldZ=[x.z for x in allPeak]
    iter=0
    while 1:
        intA, intB, edistrP, edistrN, th, pi = M_step(allPeak)
        # B region: (intA, intB]
#         print intA, intB, pi
        E_step(allPeak, intA, intB, edistrP, edistrN, th, pi)
        nowZ=[x.z for x in allPeak]
        diffZ=max([abs(nowZ[i]-oldZ[i]) for i in range(len(nowZ))])
        iter+=1
        if debug:
            print diffZ, pi
        if (abs(oldPi-pi)<epsEM and diffZ<epsEM) or iter>maxEM:
            break
        oldA=intA
        oldB=intB
        oldPi=pi
        oldZ=nowZ
    
    if debug:
        printOut(edistrP, th)
        printOut(edistrN, th)
    
    if not debug:
        for x in allPeak:
            print '\t'.join([x.id, str(x.z)])
            
def printOut(pr, th):
    for x in pr:
        print '\t'.join([str(y) for y in x])
    print '\n'
    for x in th:
        print '\t'.join([str(y) for y in x])
    print '\n'

# pi=0.3
readFile(infile, allPeak)
runEM()

# printOut(prP, thP)
# printOut(prN, thN)
# genPostPr(testfile, prP, thP, prN, thN)


