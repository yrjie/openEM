import sys,os

if len(sys.argv)<2:
    print 'Usage: peakTag'
    exit(1)

file=sys.argv[1]

interval=[85, 164, 1000000]

def getRatio(infile):
    nowId=''
    tagOne=[0]*(len(interval)+1)
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
            print '\t'.join([nowId] + [str(x) for x in tagOne])
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
    print '\t'.join([nowId] + [str(x) for x in tagOne])

getRatio(file)