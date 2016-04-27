#!/usr/local/bin/python2.7
# FixationSubber
#Ver. 0.01
#Pre.L

import argparse,  re,  os, sys
import epitopefinder.epitopes as ep
import CommonFunctions as CF

CutOff = 0.75
#==========================="Subbin' some fixation"===========================#

def subFixation(rec, snpfile):
    
    refSeq = ''
    
    for line in rec:
        refSeq = line.seq
        refID = line.id 
    refList = list(refSeq)
    snpHandle = open(snpfile)
#print refList
    for line in snpHandle:        
        lineInfo = line.split()
        if(float(lineInfo[3]) >= CutOff):
            mt = lineInfo[2].split('>')[1]
            refList[int(lineInfo[1]) - 1] = mt
            
    print '>' + refID + "-FixSub"    
    print re.sub(r'[\'\[\], ]', '', str(refList))    

#==========================="main function"===========================#

if __name__ == '__main__':
    refHash = {}    
    
    program_function = """
    *****
    FixationSubber
    
    Ver. 0.01    
    
    
    WARNING:
    Uses the old lofreq, because the new ones are in a new output format that cannot be handled
    by this tool. 
    
    
    Requires update on algorithm for parsing. Currently have no plans of implementing the new
    algorithm to perform this function.
    
    Written by PresDawgz (~w~)v
    
    *****
    """
    
    parser = argparse.ArgumentParser()
    
    #Required Arguments:    
    parser.add_argument('-ref',  '--REF_SEQ',  help = 'A fasta sequence file in nucleotide.', required = True)
    parser.add_argument('-snp',  '--SNP_FILE',  help = 'Lofreq output.', required = True)
    
    #Optional Arguments:
    parser.add_argument('-t', '--FIX_THRESHOLD', help = 'User defined fixation cut-off. Default = 0.75', type = float)
    
    if(len(sys.argv) <= 1):
        print(program_function)
        parser.print_help()
        sys.exit()
    
    args=parser.parse_args()
    
    
    #Checking Required Arguments:    
    if(not CF.okFile(args.SNP_FILE)):
        raise CF.InputError(args.SNP_FILE,  "Invalid input file: ")
    if(not CF.okFile(args.REF_SEQ)):
        raise CF.InputError(args.REF_SEQ,  "Invalid input file: ")
    else:
        ref_record = CF.getFile(args.REF_SEQ)
        CF.recordCheck(ref_record, args.REF_SEQ)            
        ref_record = CF.getFile(args.REF_SEQ) #reparse the reference file, because iterator does not reset.
    
    #Checking Optional Arguments:
    if(args.FIX_THRESHOLD):
        CutOff = args.FIX_THRESHOLD
    
    subFixation(ref_record, args.SNP_FILE)        