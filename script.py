import re

class AminoAcid:
  def __init__(self, namecode, secstr, actualstructure, jpred, porter):
    self.namecode = namecode
    self.secstr = secstr
    self.actualstructure = actualstructure
    self.jpred = jpred
    self.porter = porter



# calculate Qa
def calcqa(w,x,y,z):
    per_a = 100*w/(w + y)
    per_na = 100*x/(x + z)
    qa = (per_a + per_na)/2
    return qa


# check if an aminoacid from secstr is a given structure
def findstructureonsecstr(lineaa,structure):
    structures = {
        'T': 1,
        'B': 2,
        'H': 3  
    }
    aaarray = lineaa.split("|")
    if structure in aaarray[structures[structure]]:
        return True
    else:
        return False


# returns a list with the first value all the possible helix aas and the second value all the non possible helixaa
def findpossiblehelixaa(filename, aminoacids):
    f = open(filename, "r")
    lines = f.readlines()
    i = 0
    for line in lines:
        if re.search('^\s*[0-9]',line):
            aminoacids.append(AminoAcid(line[9], findstructureonsecstr(line, 'H'), False, False, False))
        i += 1
    return aminoacids


# on existed aminoacid list we add the structure according to dssp
def findhelicesaadssp(filename, aminoacids):
    f = open(filename, "r")
    lines = f.readlines()
    foundheader = False
    i = 0 
    for line in lines: 
        if foundheader:
            if (re.search("[HGI]", line[16])):
                aminoacids[i].actualstructure = True
            i +=1
        if not foundheader and line.startswith('  #'):
            foundheader = not foundheader
    return aminoacids


def findhelicesaajpred(filename, aminoacids):
    f = open(filename, "r")
    lines = f.readlines()
    i = 0 
    for a in lines[0]: 
        if (re.search("[A-Z]", a) and re.search("[HGI]", lines[1][i])):
            aminoacids[i].jpred = True
        i +=1
    return aminoacids


# ex.1 script: calculate PBTI Qa from Secstr
aminoacids = []
filename_secstr = "secstr_pbti.txt"
aminoacids = findpossiblehelixaa(filename_secstr, aminoacids)
filename_dssp = "pbti.dssp"
aminoacids = findhelicesaadssp(filename_dssp, aminoacids)
w = 0
x = 0
y = 0
z = 0
for aa in aminoacids:
    if aa.secstr and aa.actualstructure:
        w +=1
    elif not aa.secstr and not aa.actualstructure:
        x +=1
    elif not aa.secstr and aa.actualstructure:
        y +=1
    else:
        z +=1


qa = calcqa(w,x,y,z)
print ("result {:f}".format(qa))


# ex2 compare SECSTR, JPRED and PORTER on 1l2y protein
# read files
# aminoacids = []
# filename_secstr = "secstr_1l2y.txt"
# aminoacids = findpossiblehelixaa(filename_secstr, aminoacids)
# filename_dssp = "1l2y.dssp"
# aminoacids = findhelicesaadssp(filename_dssp, aminoacids)
# filename_jpred = "jpred_1l2y.txt"
# aminoacids = findhelicesaajpred(filename_jpred, aminoacids)


# # calculate Qa Sectr
# w = 0
# x = 0
# y = 0
# z = 0
# for aa in aminoacids:
#     if aa.secstr and aa.actualstructure:
#         w +=1
#     elif not aa.secstr and not aa.actualstructure:
#         x +=1
#     elif not aa.secstr and aa.actualstructure:
#         y +=1
#     else:
#         z +=1

# qasecstr = calcqa(w,x,y,z)
# print ("Qa Secstr: {:f}".format(qasecstr))


# # calculate Qa JPRED
# w = 0
# x = 0
# y = 0
# z = 0
# for aa in aminoacids:
#     if aa.jpred and aa.actualstructure:
#         w +=1
#     elif not aa.jpred and not aa.actualstructure:
#         x +=1
#     elif not aa.jpred and aa.actualstructure:
#         y +=1
#     else:
#         z +=1


# qajpred = calcqa(w,x,y,z)
# print ("Qa JPRED: {:f}".format(qajpred))
