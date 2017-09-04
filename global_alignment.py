import numpy as np
class Needleman_Wuncsh(object):
#match=1,mismatch=0,gap=-1
    def __init__(self,seq1,seq2):
        self.seq1=seq1
        self.seq2=seq2
        self.scores_matrix=np.empty((len(seq1)+1,len(seq2)+1))
    def __repr__(self):
        for i in range(len(self.scores_matrix)):
            self.scores_matrix[i][0]=-i
        for i in range(len(self.scores_matrix[0])):
            self.scores_matrix[0][i]=-i
        for i in range(1,len(self.scores_matrix)):
            for j in range(1,len(self.scores_matrix[0])):
                score=0
                if self.seq1[i-1]==self.seq2[j-1]:
                    score=1
                self.scores_matrix[i][j]=np.max([self.scores_matrix[i-1][j]-1,self.scores_matrix[i][j-1]-1,self.scores_matrix[i-1][j-1]+score])
        print self.scores_matrix
        return 'Result show above'
print Needleman_Wuncsh('ACGC','GACTAC')    
