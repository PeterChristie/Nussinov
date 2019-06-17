'''
 Peter Christie
 Created for CS 315, Computational Biology @ Williams College 
 Program for calculating the folding pattern of RNA using Nussinov's algorithm for conserved engergy 
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC350273/pdf/pnas00498-0062.pdf

 USAGE 
 python3 nussinov.py <sequence> <minimum loop length>

 keep sequence length below 40 for best results 
'''

import sys
import numpy as np
np.set_printoptions(threshold=np.nan)

class Nussinov:

    def __init__(self, sequence, m):
        self.sequence = sequence
        self.m = m

    def makeMatrix(self, n): 
        """
        Creates a matrix for population. Zeros added to the middle m rows
        Matrix is copied on top and bottom, so calculations only occur on one side
        """     
        scoreMatrix = np.empty((n,n))
        scoreMatrix[:] = np.NAN
        for k in range(0, m):
            for i in range(n-k):
                j = i + k
                scoreMatrix[i][j] = 0
        return scoreMatrix

    def isPair(self, pairing):
        """
        Determines whether or not two nucleotides are a pair 
        """
        matches = [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')] 
        if pairing in matches:
            return True
        else: 
            return False

    def findOptimal(self, i, j, sequence):
        """
        Function to populate a matrix with the optimal score at each position
        """
        # cases are too close to pair 
        if i >= j-m:  
            return 0

        # either i and j are paired or not paired. If not paired, we shrink the window 
        # with findOptimal(i, j-1, sequence)
        else:
            unpaired = self.findOptimal(i, j-1, sequence)
            possiblePairs = []
            for t in range(i, j-4):
                if self.isPair((sequence[t], sequence[j])):
                    possiblePairs.append(1 + self.findOptimal(i, t-1, sequence) + self.findOptimal(t+1, j-1, sequence))
            if not possiblePairs: 
                possiblePairs = [0]
            # take maximum value from list of pairs 
            # compare that to unpaired case
            paired = max(possiblePairs)

            return max(unpaired, paired)

    def traceback(self, i, j, pairList, scoreMatrix, sequence):
        """
        Funtion to scan back through the matrix and determine optimal pairing positions 
        """
        # made it to the end of the matrix 
        if j <= i:     
            return

        # j will be unpaired in this case, so score we not change if we exlude it
        elif scoreMatrix[i][j] == scoreMatrix[i][j-1]:  
            self.traceback(i, j-1, pairList, scoreMatrix, sequence)

        # j is paired in these cases
        else: 
            # check from which sub-structure the pair score is made
            tbList = []
            for x in range(i, j-m):
                if self.isPair((sequence[x], sequence[j])):
                    tbList.append(x)
            for k in tbList:
                if k-1 < 0:
                    if scoreMatrix[i][j] == scoreMatrix[k+1][j-1] + 1:
                        pairList.append((k,j))
                        self.traceback(k+1, j-1, pairList, scoreMatrix, sequence)
                        break

                # add (k,j) to the list of pairs, traceback on the two sub-structures formed by the pair 
                elif scoreMatrix[i][j] == scoreMatrix[i][k-1] + scoreMatrix[k+1][j-1] + 1:
                    pairList.append((k,j)) 
                    self.traceback(i, k-1, pairList, scoreMatrix, sequence) 
                    self.traceback(k+1, j-1, pairList, scoreMatrix, sequence)
                    break

    def compute(self, sequence):
        """
        Driver code for both the matrix population and traceback
        """
        # create matrix the same size as sequence 
        n = len(sequence)
        scoreMatrix = self.makeMatrix(n)
        pairList = []

        # populate matrix 
        for k in range(m, n): 
            for i in range(n-k):
                j = i + k
                scoreMatrix[i][j] = self.findOptimal(i,j, sequence)

        # mirror the matrix to the opposite side to avoid using zeros
        for k in range(n):
            for i in range(0,k):
                scoreMatrix[k][i] = scoreMatrix[i][k]

        self.traceback(0, n-1, pairList, scoreMatrix, sequence)
        
        # create visual representation of the folding
        output = ['.' for base in range(len(sequence))]
        for pair in pairList:
            output[pair[0]] = '('
            output[pair[1]] = ')'

        # output
        print(scoreMatrix)
        print(pairList)
        print(''.join(output))

if __name__ == "__main__":
    # run the program for a given sequence and min loop length
    sequence, m = sys.argv[1], int(sys.argv[2])
    N = Nussinov(sequence, m)
    N.compute(sequence)
    

