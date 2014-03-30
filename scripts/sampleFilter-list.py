#!/usr/bin/env python
from vcf import SampleFilter
from common import return_file_basename
import sys
import os
import argparse

#import pdb


def main():
    """ Perform sample filtering based based on samples read from a file (each samplename on seperate line) """
    usage = "usage: %prog [options]"
    parser = argparse.ArgumentParser(description='Given a gzipped vcf file and pedigree file, generate a new vcf with samples either kept/removed from LIST ')
    
    parser = argparse.ArgumentParser(description='Program description')
    parser.add_argument('-list', dest='listfile', type=str, help="file with list of samples")
    parser.add_argument('vcfile',  type=str,help='*.vcf.gz file')
    parser.add_argument('-i',dest='invert',action='store_true', help="keep samples in LIST  rather than removing from the vcf")
    #parser.add_argument('--no-invert',dest='invert',action='store_false', help="remove samples in LIST from the vcf")
    args=parser.parse_args()
    


    vcfroot, ext = os.path.splitext(args.vcfile)
    if ext == '.gz':
        vcf_basename = return_file_basename(return_file_basename(args.vcfile))
    else:
        vcf_basename = return_file_basename(args.vcfile)
        
    listfh=open(args.listfile, 'r')
    sampleList=[]
    for line in listfh:
        sampleList.append( line.strip() )
    
    

    keepstring=",".join(sampleList)
    
   
    #pdb.set_trace()
    sf = SampleFilter(infile=args.vcfile, outfile=sys.stdout, filters=keepstring, invert=args.invert)
    

if __name__ == "__main__":
    main()

