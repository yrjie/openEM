import sys, os

if len(sys.argv)<2:
    print 'Usage: roc.tab'
    exit(1)

xmax=0.06

allRoc=[]
allRoc.append([0, 0])

print sys.argv[1]

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=[float(x) for x in line.split('\t')]
    if temp[0]>xmax:
    	break
    allRoc.append([temp[0], temp[1]])
fi.close()

x1=allRoc[-1][0]
y1=allRoc[-1][1]
x2=allRoc[-2][0]
y2=allRoc[-2][1]

ymax=y2+(y1-y2)/(x1-x2)*(xmax-x2)
if ymax>1:
    ymax=1
allRoc.append([xmax, ymax])
print allRoc[-1]

auc=0
for i, x in enumerate(allRoc):
    if i==0:
    	continue
    dx=x[0]-allRoc[i-1][0]
    ymean=(x[1]+allRoc[i-1][1])/2
    auc+=dx*ymean
print auc
