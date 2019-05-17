#!/usr/bin/python
import sys
from copy import deepcopy
from Bio import SeqIO
from Bio.Alphabet.IUPAC import extended_protein
d = dict.fromkeys(list(extended_protein.letters), 0)
d["."]=0

uniqueAa=0
aaCount=0
overAverageCount=0
overAas=""
oAve = {}
aAstring=""
aAstringSub=""

del d["."]

#Define function to return subsequences
def all_subs(word):
  sub = []
  for i in range(len(word)):
    for j in range(i, len(word)):
#      if len(word[i:j+1]) == overAverageCount:
    	sub.append(word[i:j+1])
    	print word[i:j+1]
  return sub


#Define function to compare 2 dictionaries
def dicts_equal(d1,d2):
    """ return True if all keys and values are the same """
    return all(k in d2 and d1[k] == d2[k]
               for k in d1) \
        and all(k in d1 and d1[k] == d2[k]
               for k in d2)


for record in SeqIO.parse(sys.argv[1], "fasta"):
	#The sequence is in the `seq` attribute of the record
# print record.seq
	for aa in record.seq:
		d[aa] += 1
		aAstring = aAstring + aa
	if (aa != '.'):
		aaCount += 1


def ncsub(seq, s=0):
    if seq:
        x = seq[:1]
        xs = seq[1:]
        p2 = s % 2
        p1 = not p2
        return [x + ys for ys in ncsub(xs, s + p1)] + ncsub(xs, s + p2)
    else:
        return [[]] if s >= 3 else []

testString="ABCDEFG"
for i in ncsub(range (1,len(testString))):
	subString=""
	if len(i)==3:
		for j in range(0,len(i)):
			subString= subString + testString[i[j]-1]
		print  subString
