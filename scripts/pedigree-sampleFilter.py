
from vcf import SampleFilter
from PedPy import Ped, Pedfile
from common import return_file_basename
import sys
import os
import argparse




def main():
    """ Perform sample filtering based based on samples read from a *.ped file and *.vcf.gz file """
    usage = "usage: %prog [options]"
    parser = argparse.ArgumentParser(description='Given a gzipped vcf file and pedigree file, generate a new vcf with only those samples present in the pedigree (ped file) ')
    
    parser = argparse.ArgumentParser(description='Program description')
    parser.add_argument('-ped', dest='pedfile', type=str, help="*.ped file")
    parser.add_argument('vcfile',  type=str,help='*.vcf.gz file')
    parser.add_argument('--invert',dest='invert',action='store_true', help="remove samples not in the pedigree from the vcf")
    parser.add_argument('--no-invert',dest='invert',action='store_false', help="remove samples in the pedigree from the vcf")
    args=parser.parse_args()
    
    pedigree_basename=return_file_basename(args.pedfile)
    
    vcfroot, ext = os.path.splitext(args.vcfile)
    if ext == '.gz':
        vcf_basename = return_file_basename(return_file_basename(args.vcfile))
    else:
        vcf_basename = return_file_basename(args.vcfile)
        
        
    
    
    outfile=".".join([vcf_basename, pedigree_basename,'vcf'])
    
    pedobj=Pedfile(args.pedfile)
    pedobj.parsePedfile()
    keepList=pedobj.returnIndivids()
    keepstring=",".join(keepList)


    sf = SampleFilter(infile=args.vcfile, outfile=sys.stdout, filters=keepstring, invert=True)
    

if __name__ == "__main__":
    main()

