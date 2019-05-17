#!/usr/bin/python
import sys
import itertools
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

#Define function to return continuous subsequences
def all_subs(word):
	sub = []
	for i in range(len(word)):
		for j in range(i, len(word)):
			if len(word[i:j+1]) == overAverageCount:
				sub.append(word[i:j+1])
				print word[i:j+1]
	return (sub)

#Define function to return non continuous subsequences  string indicies
def ncsub(seq, s=0):
    if seq:
        x = seq[:1]
        xs = seq[1:]
        p2 = s % 2
        p1 = not p2
        return [x + ys for ys in ncsub(xs, s + p1)] + ncsub(xs, s + p2)
    else:
        return [[]] if s >= 3 else []

#Function to convert the string indicies into substring if correct length
def return_nc_sub(str):
	for i in ncsub(range (1,len(str))):
		subString=""
		sub = []
		if len(i)==oAve:
			for j in range(0,len(i)):
				subString= subString + str[i[j]-1]
			print subString
			sub.append(subString)
	return(sub)



#Define function to compare 2 dictionaries
def dicts_equal(d1,d2):
    """ return True if all keys and values are the same """
    return all(k in d2 and d1[k] == d2[k]
               for k in d1) \
        and all(k in d1 and d1[k] == d2[k]
               for k in d2)

#Define function to print all sequential substrings in a string
def find_sequential_substrings(string,substring):
	substringIndex = 0
	substrings = ""
	for i in string:
		if (substringIndex == len(substring)):
			substringIndex=0
		elif (i == substring[substringIndex]):
			substrings = substrings + i
			substringIndex+=1
	return substrings


#Read in the input
for record in SeqIO.parse(sys.argv[1], "fasta"):
        #The sequence is in the `seq` attribute of the record
#       print record.seq
        for aa in record.seq:
                d[aa] += 1
		aAstring = aAstring + aa
                if (aa != '.'):
                        aaCount += 1

#for key,val in d.items():
#        print key, val, str((float(100)/float(aaCount)* val)) +"%"

#loop through each key/val pair and store some info for later
for key,val in d.items():	
	#Count the number of Aas with a value>0
	if (val > 0):
		uniqueAa+=1

#find the average Aa number
averageAa= aaCount/uniqueAa

#Loop through each to find those overrepresented retaining the total and the number over for each aa
for key,val in d.items():
	if (val > averageAa):
		overAverageCount+= val - averageAa
		oAve[key]= val - averageAa
		overAas= overAas+ key

subD=deepcopy(oAve)
subD=dict.fromkeys(subD,0)

#Loop over the original Aa string and if the Aa is overrepresented append it to a string
for i in aAstring:
	if i in oAve.keys():
		aAstringSub= aAstringSub + i

'''
print aAstringSub
print d
print oAve
'''

subsequence = sys.argv[2]
'''
potentialSubString= find_sequential_substrings(aAstring,subsequence)
print potentialSubString
print len(potentialSubString)
'''



longestLength=0
longestSeq=""

for i in list(itertools.permutations(list(subsequence))):
	sub=find_sequential_substrings(aAstring, ''.join(i))
	if len(sub) > longestLength:
		longestLength= len(sub)
		longestSeq= sub
		longestSub= ''.join(i)

print  longestSeq
print longestLength
print longestSub
'''
########### CODE BELOW HERE IS THE NON DYNAMIC VERSION



#Loop through each substring  and check if the replaced nucleotides are equal to those which need replacing



print (return_nc_sub(aAstringSub))

#get the continuous
substrings=  all_subs(aAstringSub)

#get the non continuous
ncsub(aAstringSub)


#Loop over all substrings and check for the right aas to replace
for i in substrings:
	subDInLoop= deepcopy(subD)
	for aa in i:
		subDInLoop[aa] += 1
	if dicts_equal(subDInLoop,oAve):
		print i


#Get subsequences the length of the number overrepresented from the overrepresented subsequence
for subAa in (aAstringSub[i:i+overAverageCount] for i in range (len(aAstringSub)-overAverageCount+1)):
	subDict= dict.fromkeys(list(extended_protein.letters), 0)
	for aA in subAa:
		subDict[aA] += 1
	print subDict



#ErrorChecks
print overAverageCount
print averageAa
print oAve
print overAas
print aAstringSub



#print all subsequences of length which needs to be replaced

# check each subsequence and see if any are composed of the number of nucleotides. Can make into python dicts and just compare them# 
# see if the subsequences that do can be collapsed into a shorter subsequence, e.g. GEEKSEEKGEEK into GEEKGEEK
  

'''
