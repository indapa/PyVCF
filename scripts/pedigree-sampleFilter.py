from vcf import SampleFilter
from PedPy import Ped, Pedfile
pedfile='/Users/indapa/Research/SequencingData/Pedfiles/Families/003.ped'
vcfile='/Users/indapa/Research/SequencingData/VCFs/C868.vcf.gz'
outfile='filter.vcf'
pedobj=Pedfile(pedfile)
pedobj.parsePedfile()
keepList=pedobj.returnIndivids()
keepstring=",".join(keepList)


sf = SampleFilter(infile=vcfile, outfile=outfile, filters=keepstring, invert=True)
