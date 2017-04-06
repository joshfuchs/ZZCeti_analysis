'''
Written by JT Fuchs, UNC Chapel Hill, 2016

Combine all fitting_solutions.txt files in individual directories into one file.

Combine all pseudogaussian fits to a single PDF file.
'''

import os
from PyPDF2 import PdfFileMerger
from glob import glob

#Get list of directories in current directory. Only keep if it is a date by looking for ./2
directories = [x[0] for x in os.walk('./') if x[0][0:3]=='./2']

outpdf = PdfFileMerger()

with open('all_solutions.txt','wb') as outfile:
    for xdir in directories:
        os.chdir(xdir)
        pdfs = glob('fit*master*.pdf')
        for f in pdfs:
            outpdf.append(open(f),'rb')
        with open('fitting_solutions.txt','rb') as infile:
            outfile.write(infile.read())
        os.chdir('../')

outpdf.write(open('Pseudo_fits_flux.pdf','wb'))

