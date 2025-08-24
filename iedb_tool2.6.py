#!/usr/local/bin/python2.7
# IEDB Tool 2
#Ver. 2.6
#Pre.L

import argparse,  re,  os, sys, csv
import epitopefinder.epitopes as ep
import CommonFunctions as CF
from Bio import pairwise2
pairwise2.MAX_ALIGNMENTS = 15 # limiting pairwise2 max alignments to 15

#========================"Global Flags and Variables"=======================#

HomologyFloat = 0.8
FileSuffix_1 = '_DataBaseMatches.txt'
FileSuffix_2 = '_KnownHLA_Epitopes.txt'
FileSuffix_3 = '_UnknownHLA_Epitopes.txt'
CurrentRefID = 'Result'
WriteDIR = ''
Write = 'w'
Append = 'a'

#==========================="Parsing IEDB CSV"===========================#
# Not used in this version.

def parse_iedb_file(file):
    
    return ep.ReadCompactIEDB(file)


#==========================="What epitopes are found on the genome"===========================#

def epitopeSearch(predFile, seqHash):
    tempLineNum = 1
    predCount = 0
    predHandler = open(predFile)
    writeFile = open(WriteDIR + CurrentRefID + FileSuffix_1, 'w')
    writeFile.write("Number" + '\t' + "Epitope" + '\t' + "HLA-Type" + '\n')
    keyRecord = set([])
    
    for predLine in predHandler:
        if(predCount < 1):
            predCount += 1
            continue
        predInfo = predLine.rstrip().split('\t')
        
        for key in sorted(seqHash.keys()):
            alignment = pairwise2.align.localxx(predInfo[0], key)
            if(len(alignment) < 1):
                continue
            tempList = []
            for a in alignment:
                tempList.append(float(list(a)[2]))                        
            homologyScore = float(max(tempList) / min(len(predInfo[0]),len(key)))
            if(homologyScore >= HomologyFloat):
                if(key not in keyRecord):                 
                    keyRecord.add(key)
                    hla_list = re.sub("'", '', str(list(seqHash[key])))
                    writeFile.write(str(tempLineNum) + '\t' + key + '\t'+ hla_list +'\n')
                    tempLineNum += 1
            
    writeFile.close()


#==========================="Collect Unique Epitope"===========================#
# Not used in this version

def collectUniqueEpitope(iedbList):
        
    tempHash = {}
     
    for item in iedbList:
        if(len(item.sequence) >= 12):
            continue
        if(item.sequence not in tempHash):
            tempHash[item.sequence] = set([item.mhcallele])
        else:
            tempHash[item.sequence].add(item.mhcallele)    
    
    return tempHash    

#==========================="Read IEDB file and collect Unique Epitope"===========================#
# Replaces parse_iedb_file and collectUniqueEpitope function
# because the package epitopefinder is outdated.
# Currently under development to fit in with previous code.

def readIEDBcsv(iedbFile):

    headerCount = 2

    uniq_epiHash = {}

    with open(iedbFile, 'rb') as csvfile:
        csv_rows = csv.reader(csvfile)
        for row in csv_rows:
            if(headerCount > 0):
                headerCount -= 1
                continue
            
            allele = row[101] # The column that says Allele Name. HongAn (Column CX 101 [102-1 because python list start from 0])
            atgn_epi_rel = row[104] # The column that says Antigen Epitope Relation. Only take in if it is "Epitope" HongAn(Column DA 104)
            atgn_obj_type = row[105] # The column that says ANtigen Object Type. Only take in if it is "Linear peptide" HongAn(Column DB 105)
                   
            # Formats are (X is char, N is a int):
            # HLA-XXXN*NN:NN / HLA-XX / HLA-XXN / HLA-XXXN
            # HLA-XXXNN 
            # HLA class II
            
            atgn_desc = row[106] #The column that says Antigen Description. HongAn(Column DC 106)
            # print atgn_obj_type
            
            
            #if(count >= 0):
            if (re.search("Epitope", atgn_epi_rel) and  
                re.search("Linear peptide", atgn_obj_type)):
                
                #print allele, atgn_desc    
                if(atgn_desc in uniq_epiHash):
                    uniq_epiHash[atgn_desc].add(str(allele))
                else:
                    #print allele                    
                    uniq_epiHash[atgn_desc] = set([allele])
            #    count -= 1
            #else:
            #    break
            
    return uniq_epiHash


#==========================="Extract reference information"===========================#

def extractRefInfo(record):
    
    global CurrentRefID
    refSeq = ''
    
    for item in record:
        refSeq = item.seq
        #CurrentRefID = str(item.id)
    return refSeq


#==========================="Extract HLA information"===========================#

def extractHLA(file, *rfile):
    
    handle = open(file)
    HLA_Set = set([])
    
    if(rfile):
        rHandle = open(rfile[0])
        tempSet = set(rHandle)
        for line in handle:
            if line in tempSet:
                HLA_Set.add(line)
    else:
        HLA_Set = set(handle)
    handle.close()
    return HLA_Set


#==========================="create alternate HLA names"===========================#

def createAltTable(set_name):
    
    altHash = {}
    
    for hla in list(set_name):        
        
        if(len(hla.rstrip().split('*')) == 2):
            cleanHLA = hla.rstrip().split('*')        
            prefix = cleanHLA[0]
            suffix = cleanHLA[1].split(':')        
            altHash[hla.rstrip()] = prefix + str(int(suffix[0]))
        elif(len(hla.rstrip().split('-')) == 3):
            cleanHLA = hla.rstrip().split('-')        
            prefix = cleanHLA[0]
            mid = cleanHLA[1]
            suffix = cleanHLA[2]
            altHash[hla.rstrip()] = prefix + str(mid) +"-" + suffix
    #print altHash
    return altHash    


#==========================="Epitopes found on the genome PLUS matching HLA to subject"===========================#

def subHLA_vs_IEDB_known(alt_hash):
    
    subject_IEDBKnown = {}
    undeterminedEpi = {}
    possibleSpecial = {}
    
    for key in alt_hash:        
        lineCount = 0
        IEDB_Known = open(WriteDIR + CurrentRefID + FileSuffix_1)
        
        for line in IEDB_Known:            
            if(lineCount < 1):
                lineCount += 1
                continue
            lineList = line.rstrip().split('\t')            
            keySubbed = re.sub('\*', "\*", key) #Check the HLA to escape * character            
            
            if(re.search(keySubbed, lineList[2])): #Matches to HLA-X*XX:XX convention (e.g. HLA-A*02:01)
                if(lineList[1] not in subject_IEDBKnown):
                    subject_IEDBKnown[lineList[1]] = set([key])
                else:
                    subject_IEDBKnown[lineList[1]].add(key)                    
            
            if(re.search(alt_hash[key], lineList[2])):#Matches to HLA-XX convention (e.g. HLA-A25)                
                if(lineList[1] not in subject_IEDBKnown):
                    subject_IEDBKnown[lineList[1]] = set([alt_hash[key]])
                else:
                    subject_IEDBKnown[lineList[1]].add(alt_hash[key])
            
            if(not re.search(keySubbed, lineList[2]) and not re.search(alt_hash[key], lineList[2])):                
                if(re.search("allele undetermined", lineList[2])):                    
                    if(lineList[1] not in undeterminedEpi):                        
                        undeterminedEpi[lineList[1]] = lineList[2]
                else:                    
                    possibleSpecial[lineList[1]] = lineList[2]
            
        IEDB_Known.close()
    
    writeFile = open(WriteDIR + CurrentRefID + FileSuffix_2, 'w')
    writeUnknown = open(WriteDIR + CurrentRefID + FileSuffix_3, 'w')
    writeFile.write('Number'+ '\t' + 'Epitope' + '\t' + 'Matched HLA Type' + '\n')
    writeUnknown.write('Number'+ '\t' + 'Epitope' + '\t' + 'Matched HLA Type' + '\n')
    lineCount = 1
    unknownCount = 1
    
    for key in subject_IEDBKnown:        
        hla_list = re.sub("'", '', str(list(subject_IEDBKnown[key])))
        writeFile.write(str(lineCount) + '\t' + key + '\t' + hla_list + '\n')
        lineCount += 1
    
    for key in undeterminedEpi:        
        if(key not in subject_IEDBKnown):
            writeUnknown.write(str(unknownCount) + '\t' + key + '\t' + undeterminedEpi[key] + '\n')
            unknownCount += 1
    
    keys = undeterminedEpi.keys()
    for item in keys:
        if(item in subject_IEDBKnown):
            del undeterminedEpi[item]
            
    writeFile.close()
    writeUnknown.close()
    return (subject_IEDBKnown, undeterminedEpi, possibleSpecial)


#==========================="Predicted Epitopes vs Mutated Epitopes Prediction"===========================#
# Look for Epitope that have mutations in it.
# Look for Epitopes that didn't mutate too.

def mutPred_vs_Pred(pred_file, mpred_file, hlaTable, mutSeq):
    
    predVSmpred = {}
    noChange = {}
    predHandler = open(pred_file)    
    predCount = 0
    
    for predLine in predHandler:
        if(predCount < 1):
            predCount += 1
            continue
        predInfo = predLine.rstrip().split('\t')
        if(predInfo[1] not in hlaTable):
            continue
        mpredCount = 0
        mpredHandler = open(mpred_file)
        tempHash = {}
        tempHash2 = {}
        noPartner = True
        
        for mpredLine in mpredHandler:
            if(mpredCount < 1):
                mpredCount += 1
                continue                
            mpredInfo = mpredLine.rstrip().split('\t')
            if(mpredInfo[0] == predInfo[0]): # Dealing with epitopes that don't change
                if(predInfo[1] == mpredInfo[1]): # And their HLA matches, then we store info
                    tempHash2[predInfo[1]] = [mpredInfo[0],mpredInfo[1],predInfo[2],mpredInfo[2],
                                              predInfo[3],predInfo[4],4]
                    noPartner = False
                    break
                
            if(predInfo[3] == mpredInfo[3] and predInfo[4] == mpredInfo[4]):# Dealing with epitopes that do change
                if(predInfo[1] == mpredInfo[1]):# And their HLA matches, then we store info
                    tempHash[predInfo[1]] = [mpredInfo[0],mpredInfo[1],predInfo[2],mpredInfo[2],
                                             predInfo[3],predInfo[4],5]
                    noPartner = False
                    break
                                             
        if(noPartner):#Dealing with epitopes are predicted Immune Esc. because no prediction results
            tempHash[predInfo[1]] = [mutSeq[int(predInfo[3])-1:int(predInfo[4])],"NULL",predInfo[2],"NULL",
                                     predInfo[3],predInfo[4], 6]
        
        if(len(tempHash) > 0):
            if(predInfo[0] in predVSmpred):
                innerHash = predVSmpred[predInfo[0]]
                for key in tempHash.keys():
                    innerHash[key] = tempHash[key] #Adds new DIFFERENT HLA entry into existing hash
                predVSmpred[predInfo[0]] = innerHash
            else:
                predVSmpred[predInfo[0]] = tempHash #Initialises the hash with the first tempHash
        
        if(len(tempHash2) > 0):
            if(predInfo[0] in noChange):
                innerHash = noChange[predInfo[0]]
                for key in tempHash2.keys():
                    innerHash[key] = tempHash2[key] #Adds new DIFFERENT HLA entry into existing hash
                noChange[predInfo[0]] = innerHash
            else:
                noChange[predInfo[0]] = tempHash2
        mpredHandler.close()
    predHandler.close()    
    return predVSmpred, noChange


#==========================="Stage 2a Categorising"===========================#

def stage2a(speciHash, predChangeHash, noChange, alt_Hash):
    
    predVSknown = open(WriteDIR + CurrentRefID + '_EpitopeList.txt', 'w')        
    predVSknown.write('No.\tEpitope(K)\tEpitope(S)\tVariant\tHLA-DB\tHLA-Sub\tStart Position\tEnd Position\tFounder Rank\tPRC\tCategory\n')
    lineCount = 1
    
    #Do the epitopes in NoChange first
    for key1 in sorted(speciHash.keys()):
        for key2 in noChange:            
            homologyScore = getHomologyScore(key1, key2)
            if(homologyScore >= HomologyFloat):
                noChangeInner = noChange[key2]
                for innerKey in noChangeInner.keys():                     
                    if(innerKey in speciHash[key1] or alt_Hash[innerKey] in speciHash[key1]):
                        
                        scoreChange = float(noChangeInner[innerKey][2]) - float(noChangeInner[innerKey][3])
                        predVSknown.write(str(lineCount)+'\t'+key1+'\t'+key2+'\tNo-Change\t'+str(list(speciHash[key1])))
                        predVSknown.write('\t'+innerKey+'\t'+str(noChangeInner[innerKey][4])+'\t'+str(noChangeInner[innerKey][5]))
                        predVSknown.write('\t'+str(noChangeInner[innerKey][2])+'\t'+str(scoreChange))
                        predVSknown.write('\t'+str(noChangeInner[innerKey][6])+'\n')                        
                        lineCount+=1
            
    #Then do the epitopes that do change
    for key1 in sorted(speciHash.keys()):
        for key2 in sorted(predChangeHash.keys()):
            homologyScore = getHomologyScore(key1, key2)
            if(homologyScore >= HomologyFloat):
                predInner = predChangeHash[key2]
                for innerKey in predInner.keys():
                    if(innerKey in speciHash[key1] or alt_Hash[innerKey] in speciHash[key1]):
                        predVSknown.write(str(lineCount)+'\t'+key1+'\t'+key2+'\t'+str(predInner[innerKey][0])+'\t'+str(list(speciHash[key1])))
                        predVSknown.write('\t'+innerKey+'\t'+str(predInner[innerKey][4])+'\t'+str(predInner[innerKey][5]))
                        if(predInner[innerKey][3] != 'NULL'):
                            scoreChange = float(predInner[innerKey][2]) - float(predInner[innerKey][3])
                            predVSknown.write('\t'+str(predInner[innerKey][2])+'\t'+str(scoreChange)+'\t'+str(predInner[innerKey][6])+'\n')
                        else:
                            predVSknown.write('\t'+str(predInner[innerKey][2])+'\t'+'NULL'+'\t'+str(predInner[innerKey][6])+'\n')
                        lineCount+=1
    
    predVSknown.close()
    return lineCount


#==========================="Stage 2b Categorising"===========================#

def stage2b(funkyHash, pred_file, hlaTable, lineCount, instruction):
    
    predHandler = open(pred_file)
    writeFile = open(WriteDIR + CurrentRefID + "_EpitopeList-SingleTP.txt", instruction)
    if(instruction == 'w'):
        writeFile.write('No.\tEpitope(K)\tEpitope(S)\tVariant\tHLA-DB\tHLA-Sub\tStart Position\tEnd Position\tFounder Rank\tPRC\tCategory\n')    
    predCount = 0
    
    for predLine in predHandler:
        if(predCount < 1):
            predCount += 1
            continue
        predInfo = predLine.rstrip().split('\t')
        if(predInfo[1] not in hlaTable):
            continue
        for key in funkyHash:
            homologyScore = getHomologyScore(predInfo[0], key)
            if(type(funkyHash[key]) == set):
                hashList = list(funkyHash[key])
            else:
                hashList = funkyHash[key]
            if(homologyScore >= HomologyFloat):
                if(instruction == 'w'):
                    if(predInfo[1] in funkyHash[key] or hlaTable[predInfo[1]] in funkyHash[key]):
                        writeFile.write(str(lineCount)+'\t'+key+'\t'+predInfo[0]+'\tN/A\t'+str(hashList)+'\t'+predInfo[1]+'\t')                    
                        writeFile.write(str(predInfo[3])+'\t'+str(predInfo[4]))
                        writeFile.write('\t'+str(predInfo[2])+'\tN/A\tN/A\n')
                        lineCount += 1
                else:                    
                    writeFile.write(str(lineCount)+'\t'+key+'\t'+predInfo[0]+'\tN/A\t'+str(hashList)+'\t'+predInfo[1]+'\t')                    
                    writeFile.write(str(predInfo[3])+'\t'+str(predInfo[4]))
                    writeFile.write('\t'+str(predInfo[2])+'\tN/A\tN/A\n')
                    lineCount += 1
                
    writeFile.close()
    return lineCount


#==========================="Stage 3a Categorising"===========================#

def stage3a(unspeciHash, predChangeHash, noChange, alt_Hash, lineCount):
    
    predVSknown = open(WriteDIR + CurrentRefID + '_EpitopeList.txt', 'a')    
    
    
    #Do the epitopes in NoChange first
    for key1 in sorted(unspeciHash.keys()):
        for key2 in noChange:            
            homologyScore = getHomologyScore(key1, key2)
            if(homologyScore >= HomologyFloat):
                noChangeInner = noChange[key2]
                for innerKey in noChangeInner.keys():                     
                    if(innerKey not in unspeciHash[key1] and alt_Hash[innerKey] not in unspeciHash[key1]):
                        scoreChange = float(noChangeInner[innerKey][2]) - float(noChangeInner[innerKey][3])
                        predVSknown.write(str(lineCount)+'\t'+key1+'\t'+key2+'\tNo-Change\t'+unspeciHash[key1]+'\t')
                        #predVSknown.write(str(lineCount)+'\t'+key1+'\t'+key2+'\t'+str(noChangeInner[innerKey][0])+'\t'+innerKey+'\t')
                        predVSknown.write(innerKey+'\t'+str(noChangeInner[innerKey][4])+'\t'+str(noChangeInner[innerKey][5]))
                        predVSknown.write('\t'+str(noChangeInner[innerKey][2])+'\t'+str(scoreChange)+'\t'+str(noChangeInner[innerKey][6]-3)+'\n')                        
                        lineCount+=1
            
    #Then do the epitopes that do change
    for key1 in sorted(unspeciHash.keys()):
        for key2 in sorted(predChangeHash.keys()):
            homologyScore = getHomologyScore(key1, key2)
            if(homologyScore >= HomologyFloat):                
                predInner = predChangeHash[key2]
                for innerKey in predInner.keys():
                    if(innerKey not in unspeciHash[key1] and alt_Hash[innerKey] not in unspeciHash[key1]):                        
                        predVSknown.write(str(lineCount)+'\t'+key1+'\t'+key2+'\t'+str(predInner[innerKey][0])+'\t'+unspeciHash[key1]+'\t')
                        predVSknown.write(innerKey+'\t'+str(predInner[innerKey][4])+'\t'+str(predInner[innerKey][5]))
                        if(predInner[innerKey][3] != 'NULL'):
                            scoreChange = float(predInner[innerKey][2]) - float(predInner[innerKey][3])
                            predVSknown.write('\t'+str(predInner[innerKey][2])+'\t'+str(scoreChange)+'\t'+str(predInner[innerKey][6]-3)+'\n')
                        else:
                            predVSknown.write('\t'+str(predInner[innerKey][2])+'\t'+'NULL'+'\t'+str(predInner[innerKey][6]-3)+'\n')                            
                        lineCount+=1
    
    predVSknown.close()
    return lineCount


#==========================="Special Stage A"===========================#

def specialStageA(specialHash, predChangeHash, noChange, alt_Hash):
    
    specialFile = open(WriteDIR + CurrentRefID + '_SpecialList.txt', 'w')        
    specialFile.write('No.\tEpitope(K)\tEpitope(S)\tVariant\tHLA-DB\tHLA-Sub\tStart Position\tEnd Position\tFounder Rank\tPRC\tCategory\n')
    lineCount = 1
    
    #Do the epitopes in NoChange first
    for key1 in sorted(specialHash.keys()):
        if key1 in noChange:
            noChangeInner = noChange[key1]
            for innerKey in noChangeInner:
                if(innerKey in specialHash[key1] or alt_Hash[innerKey] in specialHash[key1]):
                    continue
                if(float(noChangeInner[innerKey][2]) <= 5.0 or float(noChangeInner[innerKey][3]) <= 5.0):
                    scoreChange = float(noChangeInner[innerKey][2]) - float(noChangeInner[innerKey][3])
                    specialFile.write(str(lineCount)+'\t'+key1+'\t'+key1+'\tNo-Change\t'+specialHash[key1]+'\t'+innerKey+'\t')                        
                    specialFile.write(str(noChangeInner[innerKey][4])+'\t'+str(noChangeInner[innerKey][5]))
                    specialFile.write('\t'+str(noChangeInner[innerKey][2])+'\t'+str(scoreChange)+'\t'+'7\n')                        
                    lineCount+=1
            
    #Then do the epitopes that do change
    for key1 in sorted(specialHash.keys()):
        if key1 in predChangeHash:
            predInner = predChangeHash[key1]
            for innerKey in predInner.keys():
                if(innerKey in specialHash[key1] or alt_Hash[innerKey] in specialHash[key1]):
                    continue
                if(predInner[innerKey][3] != 'NULL'):
                    if(float(predInner[innerKey][2]) <= 5.0 or float(predInner[innerKey][3]) <= 5.0):
                        specialFile.write(str(lineCount)+'\t'+key1+'\t'+key1+'\t'+str(predInner[innerKey][0])+'\t'+specialHash[key1]+'\t'+innerKey+'\t') 
                        specialFile.write(str(predInner[innerKey][4])+'\t'+str(predInner[innerKey][5]))
                        scoreChange = float(predInner[innerKey][2]) - float(predInner[innerKey][3])
                        specialFile.write('\t'+str(predInner[innerKey][2])+'\t'+str(scoreChange)+'\t'+'7\n')
                else:
                    if(float(predInner[innerKey][2]) <= 5.0):
                        specialFile.write(str(lineCount)+'\t'+key1+'\t'+key1+'\t'+str(predInner[innerKey][0])+'\t'+specialHash[key1]+'\t'+innerKey+'\t') 
                        specialFile.write(str(predInner[innerKey][4])+'\t'+str(predInner[innerKey][5]))
                        specialFile.write('\t'+str(predInner[innerKey][2])+'\t'+'NULL'+'\t'+'7\n')
                lineCount+=1
    
    specialFile.close()


#==========================="Special Stage B"===========================#

def specialStageB(specialHash, predFile, alt_Hash):
    
    predHandler = open(predFile)
    specialFile = open(WriteDIR + CurrentRefID + '_SpecialSingleTP.txt', 'w')
    specialFile.write('No.\tEpitope(K)\tEpitope(S)\tVariant\tHLA-DB\tHLA-Sub\tStart Position\tEnd Position\tFounder Rank\tPRC\tCategory\n')
    lineCount = 1
    predCount = 0
    
    for predLine in predHandler:
        if(predCount < 1):
            predCount += 1
            continue
        predInfo = predLine.rstrip().split('\t')    
    
        #Check the specialHash
        if(predInfo[0] in specialHash):
            if(predInfo[1] not in alt_Hash):
                continue
            if(predInfo[1] in specialHash[predInfo[0]] or alt_Hash[predInfo[1]] in specialHash[predInfo[0]]):
                continue
            hashList = specialHash[predInfo[0]]
            if(float(predInfo[2]) <= 5.0):
                specialFile.write(str(lineCount)+'\t'+predInfo[0]+'\t'+predInfo[0]+'\tN/A\t'+str(hashList)+'\t'+predInfo[1]+'\t')                    
                specialFile.write(str(predInfo[3])+'\t'+str(predInfo[4]))
                specialFile.write('\t'+str(predInfo[2])+'\tN/A\tN/A\n')                      
                lineCount+=1
    specialFile.close()


#==========================="Retrieve alignment and homology score"===========================#

def getHomologyScore(string1, string2):
    
    alignment = pairwise2.align.localxx(string1, string2)
    homologyScore = 0
        
    tempList = []
    if(len(alignment) >= 1):
        for a in alignment:
            tempList.append(float(list(a)[2]))                        
        homologyScore = float(max(tempList) / len(string1)) 
    return homologyScore


#==========================="main function"===========================#

if __name__ == '__main__':
    
    refHash = {}
    
    program_function = """
    *****
    
    IEDB Tool 2
    
    Processes IEDB data in "compact" format (option available in IEDB) and IEDB epitope prediction
    algorithm outputs. THis tool has two modes.
    
    ***Mode 1 - Single Time Point Data Processing***
        By providing only:
            1) The .csv file from IEDB.
            2) Prediction file parsed using iedbPredictionParser tool.
            3) HLA allele list relevant to subject.
        
        This tool will produce a list of epitopes matching with what is found in IEDB with a user
        defined homology (defaults at 0.8).
    
    ***Mode 2 - Double Time Point Data Processing***
        By providing:
            1) The .csv file from IEDB.
            2) Prediction file parsed using iedbPredictionParser tool.
            3) HLA allele list relevant to subject.
            4) Fasta file in amino acids with mutations subbed in (relative to prior time point).
            5) Prediction file of the above fasta file parsed using iedbPredictedParser tool.
            
        This tool will produce a list of epitopes that are sorted into 6 categories:
            
            Category 1-3: IEDB lists the epitope to have a HLA with an undetermined allele.
            Cat. 1 - Epitopes that do not undergo mutation (in 2nd provided time point)
            Cat. 2 - Epitopes that mutate (in 2nd provided time point) but still predicted
                     to be an epitope.
            Cat. 3 - Epitopes that mutate (in 2nd provided time point) but NOT predicted to
                     be an epitope.
            
            Category 4-7: IEDB clearly defines HLA alleles for the epitope.
            Cat. 4 - Epitopes that do not undergo mutation (in 2nd provided time point)
            Cat. 5 - Epitopes that mutate (in 2nd provided time point) but still predicted
                     to be an epitope.
            Cat. 6 - Epitopes that mutate (in 2nd provided time point) but NOT predicted to
                     be an epitope.
            Cat. 7 - Epitope that have been predicted to be good targets (Percentile Rank < 5.0)
                     however IEDB reports a non-matching HLA type for that specific epitope.
                
    Ver. 2.6
    
    Written by PresDawgz (~w~)v
    
    *****
    """
    
    parser = argparse.ArgumentParser()
    
    #Required Arguments:
    parser.add_argument('-csv', '--IEDB_CSV', help ='Filename of compact file from IEDB.', required = True)
    parser.add_argument('-pfile', '--IEDB_PRED', help = 'Parsed prediction file from iedbPredictionParser.', required = True )
    parser.add_argument('-hfile', '--HLA_LIST', help = 'A file containing subject HLA information. See template => SAMPLE-HLA_Set.txt', required = True)
    
    #Optional Arguments:
    parser.add_argument('-res', '--R_LIST', help = 'File containing a list of HLA that can be synethesized into Dextramers')
    parser.add_argument('-hv', '--H_VALUE', help = 'Percentage of homology between IEDB Known epitope and IEDB Prediction. Default = 0.8', type = float)
    parser.add_argument('-o', '--OUT_DIR', help = 'Output directory. Otherwise files are written to current directory.')
    parser.add_argument('-extra', '--EXTRA_FILES', help = 'Prevents removal of potentially useful output files that iedb_tool.py makes.', action = 'store_true' )
    parser.add_argument('-mfile', '--IEDB_MPRED', help = 'Parsed prediction file of mutated epitopes from iedbPredictionParser.')
    parser.add_argument('-mref', '--MUT_SEQ', help = 'A fasta sequence file in amino acids with mutations subbed in. Output of fixationSubber.py')
    parser.add_argument('-n', '--OUT_NAME', help = 'Prefix to add to output file name. Otherwise names files with default prefix "Result"')
    
    
    if(len(sys.argv) <= 1):
        print(program_function)
        parser.print_help()
        sys.exit()
    
    args = parser.parse_args()
    
    #Checking Required Arguments:
    if(not CF.okFile(args.IEDB_CSV)):
        raise CF.InputError(args.IEDB_CSV,  "Invalid input file: ")
    if(not CF.okFile(args.IEDB_PRED)):
        raise CF.InputError(args.IEDB_PRED,  "Invalid input file: ")
    if(not CF.okFile(args.HLA_LIST)):
        raise CF.InputError(args.HLA_LIST,  "Invalid input file: ")    
    
    #Checking Optional Arguments:
    if(args.R_LIST):
        if(not CF.okFile(args.R_LIST)):
            raise CF.InputError(args.R_LIST, "Invalid input file: ")
    if(args.H_VALUE):
        HomologyFloat = args.H_VALUE
    if(args.OUT_DIR):
        if(os.path.isdir(args.OUT_DIR)):
            WriteDIR = args.OUT_DIR + '/'
        else:
            raise CF.InputError(args.OUT_DIR, "Invalid directory: ")
    if(args.IEDB_MPRED):
        if(not CF.okFile(args.IEDB_MPRED)):
            raise CF.InputError(args.IEDB_MPRED,  "Invalid input file: ")        
    if(args.MUT_SEQ):
        if(not CF.okFile(args.MUT_SEQ)):
            raise CF.InputError(args.MUT_SEQ,  "Invalid input file: ")
        else:
            mut_record = CF.getFile(args.MUT_SEQ)
            CF.recordCheck(mut_record, args.MUT_SEQ)
    if(args.OUT_NAME):
        CurrentRefID = args.OUT_NAME
    
    ###########
    # Stage 1 # => Parsing input files and splitting specified and unspecified HLA info from IEDB database
    ###########
    
    #ref_record = CF.getFile(args.REF_SEQ) #reparse the reference file, because iterator does not reset.
    #refSeq = extractRefInfo(ref_record)    
    #iedbList = parse_iedb_file(args.IEDB_CSV) #iedibList contains info from IEDB known epitops.
    #uniqueEpitopeHash = collectUniqueEpitope(iedbList)#uniqueEpitopeHash contains a set of non-duplicated epitopes and their list of HLA types.    
    
    uniqueEpitopeHash = readIEDBcsv(args.IEDB_CSV)
    
    if(args.R_LIST):
        HLA_Set = extractHLA(args.HLA_LIST, args.R_LIST)
    else:
        HLA_Set = extractHLA(args.HLA_LIST) #Store the HLA information of the subject from input.    
    
    HLA_alt_table = createAltTable(HLA_Set) #Make alternative table because some people store HLA-A*02:01 and HLA-A*02:06 as HLA-A2
    epitopeSearch(args.IEDB_PRED, uniqueEpitopeHash)
    (speciHash, unspeciHash, specialHash) = subHLA_vs_IEDB_known(HLA_alt_table) #returns hashtables
    epitopeNumStart = 1
    
    ############
    # Stage 2a # => Find the top 50 or so epitopes that are under specified HLA. Those that don't change and change
    ############
    if(args.MUT_SEQ and args.IEDB_MPRED):    
        mut_record = CF.getFile(args.MUT_SEQ) 
        mutSeq = extractRefInfo(mut_record)
    
        predictionChange, noChange = mutPred_vs_Pred(args.IEDB_PRED, args.IEDB_MPRED, HLA_alt_table, mutSeq)    
        epitopeNum = stage2a(speciHash, predictionChange, noChange, HLA_alt_table)
    
    ############
    # Stage 3a # => If we do not have 50. Then we start searching the unspecified ones.
    ############
        #if(epitopeNum < 50):
        #    pass
        epitopeNum2 = stage3a(unspeciHash, predictionChange, noChange, HLA_alt_table, epitopeNum)
    else:
    ############
    # Stage 2b # => Find the top 50 or so epitopes that have specified HLA using only one time point of data.
    ############    
        epitopeNum = stage2b(speciHash, args.IEDB_PRED, HLA_alt_table, epitopeNumStart, Write)
        
        #if(eptiopeNum < 50):
        #    pass
        epitopeNum2 = stage2b(unspeciHash, args.IEDB_PRED, HLA_alt_table, epitopeNum, Append)    
    
    
    #################
    # Special Stage # => Epitopes with really high ranking prediction. This stage disregards HLA matching
    #################    and only considers 100% match in sequence
    
    # Reminder to self: does not require unspecified because it would've turned up previously
    # as allele undetermined
        
    if(args.MUT_SEQ and args.IEDB_MPRED):
        specialStageA(specialHash, predictionChange, noChange, HLA_alt_table)
    else:
        #print specialHash
        specialStageB(specialHash, args.IEDB_PRED, HLA_alt_table)
    if(not args.EXTRA_FILES):
        os.remove(WriteDIR + CurrentRefID + FileSuffix_1)
        os.remove(WriteDIR + CurrentRefID + FileSuffix_2)
        os.remove(WriteDIR + CurrentRefID + FileSuffix_3)
    
