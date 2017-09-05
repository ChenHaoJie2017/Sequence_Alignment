'''Scoring Matrix:
+---+---+---+---+---+---+
|   | - | A | T | G | C |
+---+---+---+---+---+---+
| - | 8 | -6| -6| -6| -6|
+---+---+---+---+---+---+
| A | -6| 8 | -4| -2| -4|
+---+---+---+---+---+---+
| T | -6| -4| 8 | -4| -2|
+---+---+---+---+---+---+
| G | -6| -2| -4| 8 | -4|
+---+---+---+---+---+---+
| C | -6|-4 | -2| -4| 8 |
+---+---+---+---+---+---+
'''


#I develop this program for research study.free and open for this codes.
#I develop this program for research study.free and open for this codes.
import math
class Needleman_Wunsch:

#initializes the algorithm
    def __init__(self,s,t):
        self.s=s
        self.t=t
#to get the result
    def Start(self):
        print '**********************The Needleman Wunsch Algorithm**********************'
        print ''
        print '***************************Global Alignment*******************************'
        print ''
        print '=========================================================================='
        h,d=HD(self.s,self.t)
        i=1
        print 'result '+str(i)+'~'+str(i+1)+':'
        while getBack(self.s,self.t,d)!=getBack(self.s,self.t,d):
            i+=2
            print 'result '+str(i)+'~'+str(i+1)+':'
        print '**********************Author_info**********************'
        print 'Author:Chen_Haojie,If your have some problem or advices '+'\n'+'you can connect me with QQ(2631838120) or WeChat(CHJ2631838120) '
        print '-----------------Free and Open-----------------'



#calculate the maximum value and return a direction for tracing back.
def judge(number,max_value):
    if number==max_value:
        return '1'
    else:
        return '0'
#help to construct two matrixs
def maximum(one,two,three):
    max1=max(one,two,three)
    pos=0
    dic={'100':1,'010':2,'001':3,'110':4,'011':5,'101':6,'111':7}
    check=''
    array=[one,two,three]
    for i in array:
        check+=judge(i,max1)
    return max1,dic[check]
    
#get the H score matrix
#get the D direction matrix
def HD(s,t):
    import numpy as np
    h=np.zeros((len(s)+1,len(t)+1))
    d=np.zeros((len(s)+1,len(t)+1))
    for i in range(len(s)+1):
        h[i][0]=-2*i
        if i!=0:
            d[i][0]=1
    for i in range(len(t)+1):
        h[0][i]=-2*i
        if i!=0:
            d[0][i]=3
    for i in range(1,len(s)+1):
        for j in range(1,len(t)+1):
            score=-1
            if s[i-1]==t[j-1]:
                score=1
            h[i][j],d[i][j]=maximum(h[i-1][j]+(-2),h[i-1][j-1]+score,h[i][j-1]+(-2))
    return h,d

#trace back
def getBack(s,t,d):
    s1=''
    t1=''
    i=len(s)
    j=len(t)
    last=[]
    while i>0 and j>0:
        k=d[i][j]
        if k==1:
            s1=s[i-1]+s1
            t1='-'+t1
            i-=1
            s1,t1=edge(i,j,s,t,s1,t1)
        if k==2:
            s1=s[i-1]+s1
            t1=t[j-1]+t1
            i-=1
            j-=1
            s1,t1=edge(i,j,s,t,s1,t1)
        if k==3:
            s1='-'+s1
            t1=t[j-1]+t1
            j-=1
            s1,t1=edge(i,j,s,t,s1,t1)
        if k==4:
            s1=s[i-1]+s1
            t1='-'+t1
            last.append([i,j])
            i-=1
            s1,t1=edge(i,j,s,t,s1,t1)
        if k==5:
            s1='-'+s1
            t1=t[j-1]+t1
            last.append([i,j])
            j-=1
            s1,t1=edge(i,j,s,t,s1,t1)
        if k==6:
            s1=s[i-1]+s1
            t1='-'+t1
            last.append([i,j])
            i-=1
            s1,t1=edge(i,j,s,t,s1,t1)
        if k==7:
            s1=s[i-1]+s1
            t1='-'+t1
            last.append([i,j])
            i-=1
            s1,t1=edge(i,j,s,t,s1,t1)
    d=modify(d,last)
    between=process(s1,t1)
    show(s1,between,t1)
    return s1,between,t1
    
#process the result
def process(s1,t1):
    between=''
    for i in range(len(s1)):
        if s1[i]==t1[i]:
            between+='|'
        elif s1[i]!='-'and t1[i]!='-':
            between+=':'
        else:
            between+=' '
    return between

#display the result
def show(s1,between,t1):
    count=0
    for i in range(int(math.floor(len(s1)/float(60)))):
        print s1[i*60:(60+i*60)]+'\n'+between[i*60:(60+i*60)]+'\n'+t1[i*60:(60+i*60)]
        count=60+i*60
    print s1[count:]+'\n'+between[count:]+'\n'+t1[count:]    
#to recognize the edge
def edge(i,j,s,t,s1,t1):
    if j==0 and i!=0:
        s1=s[:i]+s1
        t1='-'*i+t1
    if j!=0 and i==0:
        s1='-'*j+s1
        t1=t[:j]+t1
    return s1,t1      
#to modify the direction matrix
def modify(d,last):
    if len(last)>0:
        x=last[-1][0]
        y=last[-1][1]
        if d[x][y]==4:
            d[x][y]=2
        elif d[x][y]==5:
            d[x][y]=2
        elif d[x][y]==6:
            d[x][y]=3
        elif d[x][y]==7:
            d[x][y]=5
    return d

s=s_fabH_Ecoli_k12='ATGTATACGAAGATTATTGGTACTGGCAGCTATCTGCCCGAACAAGTGCGGACAAACGCCGATTTGGAAAAAATGGTGGACACCTCTGACGAGTGGATTGTCACTCGTACCGGTATCCGCGAACGCCACATTGCCGCGCCAAACGAAACCGTTTCAACCATGGGCTTTGAAGCGGCGACACGCGCAATTGAGATGGCGGGCATTGAGAAAGACCAGATTGGCCTGATCGTTGTGGCAACGACTTCTGCTACGCACGCTTTCCCGAGCGCAGCTTGTCAGATTCAAAGCATGTTGGGCATTAAAGGTTGCCCGGCATTTGACGTTGCAGCAGCCTGCGCAGGTTTCACCTATGCATTAAGCGTAGCCGATCAATACGTGAAATCTGGGGCGGTGAAGTATGCTCTGGTCGTCGGTTCCGATGTACTGGCGCGCACCTGCGATCCAACCGATCGTGGGACTATTATTATTTTTGGCGATGGCGCGGGCGCTGCGGTGCTGGCTGCCTCTGAAGAGCCGGGAATCATTTCCACCCATCTGCATGCCGACGGTAGTTATGGTGAATTGCTGACGCTGCCAAACGCCGACCGCGTGAATCCAGAGAATTCAATTCATCTGACGATGGCGGGCAACGAAGTCTTCAAGGTTGCGGTAACGGAACTGGCGCACATCGTTGATGAGACGCTGGCGGCGAATAATCTTGACCGTTCTCAACTGGACTGGCTGGTTCCGCATCAGGCTAACCTGCGTATTATCAGTGCAACGGCGAAAAAACTCGGTATGTCTATGGATAATGTCGTGGTGACGCTGGATCGCCACGGTAATACCTCTGCGGCCTCTGTCCCGTGCGCGCTGGATGAAGCTGTACGCGACGGGCGCATTAAGCCGGGGCAGTTGGTTCTGCTTGAAGCCTTTGGCGGTGGATTCACCTGGGGCTCCGCGCTGGTTCGTTTCTAG'
t=t_fabH_Xoo='ATGAGCAAGCGGATCTATTCCAGAATCGCGGGCACGGGTAGCTATTTGCCCGAAAAGGTGTTGACCAACGACGACATGTCGAAGATCGTCGACACCAGCGATGAGTGGATCTTCTCGCGCACTGGTATTCGTGAACGTCACATCGCTGCCGACGACCAGGCCACCAGCGACCTGGCGTATTTCGCCTCGCTGAAGGCGATGGAAGCGGCCGGTGTGACTGCCGACGAGATCGATCTGATCGTCATCGGCACCACCACCCCGGACCTGATATTTCCGTCGACCGCCTGTCTGCTGCAGGCGCGGCTGGGCAACGTCGGTTGCGGTGCAATGGACGTCAACGCTGCCTGTTCGGGGTTTGTCTACGCGCTTAGCGTGGCCGACAAGTTCGTGCGCAGCGGCGACGCCAAGACCGTGCTGGTGGTTGGCGCAGAAACGCTGACCCGCATCCTCGACTGGACCGACCGCACCACCTGCGTGCTATTCGGCGATGGTGCCGGCGCGGTGATCCTCAAGGCCGACGAAGAGACCGGTATTCTCAGCACCCACTTGCGTGCCGATGGCAGCAAGAAAGAGCTGCTGTGGGATCCGGTTGGCGTGTCGGTTGGCTTCGGTGAAGGCAAGAATGGCGGCGGCGCGTTGCTGATGAAAGGCAACGACGTCTTCAAATACGCGGTCAAGGCGCTGAACTCGGTGGTCGATGAGACTCTGGCGGCCAACGGCTACGACGCGCACGACCTGGATTGGCTGATTCCACATCAGGCCAATCTGCGCATCATCGAAGCCACCGCCAAGCGGCTGGACCTGCCGATGGAGCAAGTGGTGGTCACCGTGGATCGCCACGGCAATACTTCCTCGGCGTCGGTGCCGCTGGCATTGGACGAAGCGGTGCGCTCGGGTCGCGTGCAGCGCGGCCAGTTGCTGCTGCTGGAAGCCTTCGGCGGCGGGTTTGCCTGGGGCTCGGCGCTGTTGCGCTACTGA'

test= Needleman_Wunsch(s,t)
test.Start()


