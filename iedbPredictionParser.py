#!/usr/local/bin/python2.7
# IEDB Prediction Parser
#Ver. 0.01
#Pre.L

import argparse,  re,  os, sys
import CommonFunctions as CF

#==========================="Global Variables"===========================#

RankingCutOffLower = 10.0 # 10%
RankingCutOffUpper = 0.0 # 0%
OutputPath = ''

#==========================="Parsing IEDB Recommended Prediction"===========================#
# record[0], record[5], record[7] are allele, pattern and rank.
# record[2] and record[3] are start and stop positions respectively.

def parse_iedb_rec(file, cutoffLOWER, cutoffUPPER):
    
    fileHandle = open(file)
    patternHash = {}
    lineCount = 0    
    HLA_Set = set([])
    
    for line in fileHandle:
        if(lineCount < 4):
            lineCount += 1
            continue
        record = line.rstrip().split('\t')
        if (len(record) > 1):            
            if(float(record[7]) > cutoffLOWER):
                continue
            if(float(record[7]) < cutoffUPPER):
                continue    
            if(record[5] in patternHash):
                patternHash[record[5]].append(CF.Epitope_Info(record[0], record[7], record[2], record[3]))
            else:
                patternHash[record[5]] = [CF.Epitope_Info(record[0], record[7], record[2], record[3])]
            if(record[0] not in HLA_Set):
                HLA_Set.add(record[0])
                      
    return (patternHash, HLA_Set)

#==========================="Output hash to File"===========================#

def outputHash(hash):
    
    writeHandle = open(OutputPath + '.txt', 'w')
    writeHandle.write('Epitope' + '\t' + 'HLA-Type' + '\t' + 'Percentile Rank' 
                                +'\t' + 'Start Position'+ '\t' + 'End Position' + '\n')
    for key in sorted(hash):
        for item in hash[key]:
            writeHandle.write(key + '\t' + item.HLA_type() + '\t' + item.percentile_rank() 
                                +'\t' + item.start_position()+ '\t' + item.end_position() + '\n')
                                
    writeHandle.close()
            
#print key, "\t", item.HLA_type(), "\t", item.percentile_rank(), item.start_position(), item.end_position()


#==========================="Output HLA information to File"===========================#

def outputHLA(HLA_Set):
    
    HLA_List = list(HLA_Set)
    writeHandle = open(OutputPath + 'HLA_Set.txt', 'w')
    
    for item in HLA_List:
        writeHandle.write(item + '\n')
        
    writeHandle.close()
    
#==========================="Main"===========================#

if __name__ == '__main__':
    
    program_function = """
    *****

    IEDB Prediction File Parser
    
    Processes IEDB prediction text files from the Recommended method.
    
    
    *****
    """
    
    parser = argparse.ArgumentParser()
    
    #Required Arguments
    parser.add_argument('-p', '--IEDB_prediction', help ='Filename the text file from IEDB MHC Class-I Prediction tool.',  required = True)
    parser.add_argument('-o', '--output_name', help = 'Output name for the text file.', required = True)    
    #Optional Arguments
    parser.add_argument('-lower',  '--ranking_threshold_lower',  help ='Cut off of epitopes at a defined percentile ranking (Range 0-100). Default = 10.0', type = float)
    parser.add_argument('-upper',  '--ranking_threshold_upper',  help ='Cut off of epitopes at a defined percentile ranking (Range 0-100). Default = 0.0', type = float)
    
    if(len(sys.argv) <= 1):
        print(program_function)
        parser.print_help()
        sys.exit()
    args=parser.parse_args()
        
    #Checking Required Arguments
    if(not CF.okFile(args.IEDB_prediction)):
        raise CF.InputError(args.IEDB_prediction,  "Invalid input file: ")    
    OutputPath = args.output_name
    
    
    #Checking Optional Arguments
    if(args.ranking_threshold_lower):
        RankingCutOffLower = args.ranking_threshold_lower
    if(args.ranking_threshold_upper):
        RankingCutOffUpper = args.ranking_threshold_upper
                
    (patternHash, HLA_Set) = parse_iedb_rec(args.IEDB_prediction, RankingCutOffLower, RankingCutOffUpper)
    outputHash(patternHash)
    outputHLA(HLA_Set)
    
    
        
    
    
