#!/usr/bin/env python
import itertools, argparse
import vcf
import numpy as np
from VcfSampleEval import *
import sys
from common import grouper, melt_lol

""" generator function that yields vcf.model._Record objects """
def py_recordgen(reader_obj):
    for rec in reader_obj:
        yield rec

""" This program takes in two (bg)zipped vcf file. -goldVcf are the gold standard genotype records
    and -evalVcf are the genotypes you want to evaluate against the gold standard.
    
    The program assumes each VCF has the same Vcf records. It will compute per-sample genotype concordance metrics
    non-reference sensitivity and non-reference discrepancy as described in this paper:
    http://www.nature.com/ng/journal/v43/n5/full/ng.806.html%3FWT.ec_id%3DNG-201105
    
    The output is written to STDOUT with sample, NRS, and NRD
    
    It will also print out csv files for each sample that represents the genotype concordance matrix cell counts
    They can be visualized in R  with fluctuation plots using this 
    script: https://github.com/indapa/BioinfoPipeline/blob/master/Rscripts/make-fluctuation-plots.R 
    
    TODO: make the script more general so you don't have to assume each VCF has the same records"""

def main():
    usage = "usage: %prog [options]  "
    
    parser = argparse.ArgumentParser(description='Calculate non-reference sensitivity (NRS) and non-reference discrepancy (NRD) of VCF files with the same records')
    parser.add_argument("-goldvcf", dest='gold', help="VCF with gold standard genotypes you want to compare to")
    parser.add_argument("-evalvcf", dest='eval', help="VCF you want to evaluate against the gold standard")
   
    args=parser.parse_args()
                                            

    nrsfh=open('NRS.log', 'w')
    nrdfh=open('NRD.log', 'w')
    #matrixfh=open('overall.wes.array.genotype.matrix.csv', 'w')
    concordancetable= np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )
    calledtable = np.matrix ( [ [0 ,0] , [0,0] ] )


    vcf_readerOne = vcf.Reader(open(args.eval, 'r'),compressed=True)
    vcf_readerTwo = vcf.Reader(open(args.gold, 'r'),compressed=True)




    vcf_gen1=py_recordgen(vcf_readerOne)
    vcf_gen2=py_recordgen(vcf_readerTwo)


    FLAG=True
    vcf_sample_eval_objects=[]
    common_samples=[]


    sys.stderr.write("computing per-sample concordance ....\n")
    for vrec1, vrec2 in itertools.izip(vcf_gen1, vcf_gen2):

        vrec1_samples=[elem.sample for elem in vrec1.samples ]
        vrec2_samples=[elem.sample for elem in vrec2.samples ]

        if vrec1.CHROM != vrec2.CHROM:
            sys.stderr.write("chromosome number does not match!\n")
            sys.stderr.write(vrec1.CHROM + " " + vrec2.CHROM + "\n")
            sys.exit(1)

        if vrec1.POS != vrec2.POS:
            sys.stderr.write("chromosome POS  does not match!\n")
            sys.stderr.write(vrec1.POS + " " + vrec2.POS + "\n")
            sys.exit(1)
    
        if vrec1.ID == None: vrec1.ID='.'
        if vrec2.ID == None: vrec2.ID='.'
    
        common_samples= [x for x in vrec1_samples if x in vrec2_samples  ]
        #print len(common_samples)
    
        if FLAG == True:
            #vcf_sample_eval_objects = [ VcfSampleEval ('array', 'wes', x) for x in common_samples ]
            vcf_sample_eval_objects = [ VcfSampleEval ('gold', 'eval', x) for x in common_samples ]
       
            FLAG=False
    
        gold_eval_genotypes=[] # list of tuples where (sample_name, eval_gt, compare(gold).gt is the order
        for s in common_samples:
       
            gold_eval_genotypes.append( [s, vrec1.genotype(s).gt_type, vrec2.genotype(s).gt_type] )
        
        
    #print gold_eval_genotypes
    
        for eval_obj, eval_genotypes in itertools.izip(vcf_sample_eval_objects, gold_eval_genotypes):
            if eval_genotypes[1] == None:
                eval_genotypes[1]=3
        
            if eval_genotypes[2] == None:
                eval_genotypes[2] = 3
        
            eval_obj.incrementcellcount(eval_genotypes[1],eval_genotypes[2])
            concordancetable[eval_genotypes[1], eval_genotypes[2] ]+=1
        
            if eval_genotypes[1] != eval_genotypes[2]:
                if (eval_genotypes[1] == 0 or eval_genotypes[1] == 3) and (eval_genotypes[2] == 1 or eval_genotypes[2] == 2):
                    nrsout="\t".join( [str(vrec1.CHROM), str(vrec1.POS),eval_genotypes[0], vrec1.ID, eval_genotypes[0], str(eval_genotypes[1]), str(eval_genotypes[2]) ] )
                    nrsfh.write(nrsout+"\n")
                if eval_genotypes[1] != 3:
                    nrdout="\t".join([str(vrec1.CHROM), str(vrec1.POS),eval_genotypes[0], vrec1.ID, eval_genotypes[0], str(eval_genotypes[1]), str(eval_genotypes[2]) ])
                    nrdfh.write(nrdout+"\n")
        #print
    
    concordancefh=open("concordance.txt", 'w')
    matrixfh=open("genotype.matrix.csv", 'w')

    print "Sample\tNRS\tNRD"
    for (eval_obj, sample)  in itertools.izip(vcf_sample_eval_objects, common_samples):
        (NRS, NRD)=eval_obj.returnNRS_NRD()
        outstring="\t".join( [sample, str(NRS), str(NRD)])
        print outstring
        eval_obj.write_genotype_matrix()
    
    outstring=",".join( map(str,melt_lol(concordancetable.tolist())) )  
    #print outstring
    #matrixfh.write(outstring+"\n")

if __name__ == "__main__":
    main()
