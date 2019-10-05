# -*- coding: utf-8 -*- 
#!/usr/bin/env python 
import models.CPC2_beta.bin.CPC2 as CPC2 
import models.CPAT.bin.cpat as CPAT 
import models.PredLnc_GFStack.src.PredLnc_GFStack as GFStack 
#import models.longdist.bin.longdist as longdist
#import models.CPPred.bin.CPPred as CPPred
#import models.CNCI.CNCI as CNCI
import argparse
import os

def parse_args():
    
    parser = argparse.ArgumentParser(description = 'Lncrna Prediction Package')
#common argument
    parser.add_argument('-m', '--model', help='please choose a prediction model ')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name', type=str) 
    parser.add_argument('-i', '--input', dest='fasta', help='input file', default='fasta', type=str)
    parser.add_argument('-manual', action='store_true')
    parser.add_argument('-p', '--species', dest='species', help="choose a species") 

#CNCI group argument    
    parser.add_argument('--parallel', action='store_true', default=True,  help='assign the running CUP numbers')
    parser.add_argument('--gtf', dest='gtf', action='store_true',  help='please enter your gtf files')
    parser.add_argument('--directory', dest='directory',action='store',help='if your input file is gtf type please enter RefGenome directory')

#CPC2 group arguments
    parser.add_argument('-r', '--reverse', help="also check the reverse strand")

#CPAT group arguments
    parser.add_argument("-s","--start",action="store",dest="start_codons",default='ATG',help="Start codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). default=%default")
    parser.add_argument("-t","--stop",action="store",dest="stop_codons",default='TAG,TAA,TGA',help="Stop codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). Multiple stop codons should be separated by ','. default=%default")

#longdist group arguments
    parser.add_argument('-l', '--longdistspecies', dest='longdist_species', help="please choose a species type from Human, Mouse, Zebrafish")
    parser.add_argument('--purge', action="store_true", help='Purge all intermediate files. Intermediate files have .longdist in their names. Do not purge them if you want to run this method a second time with the same data')
    parser.add_argument('-z', '--size', nargs=1, metavar='<200>', default=[200], type=int, dest='size', help='Minimun sequence size to consider. Default is 200.')

    args = parser.parse_args()
    return args

modelDir = os.getcwd() + '/models/'

def getModel(model, saveDir):
    print("----------------------------------Start Loading %s models-------------------------"%model)
    downURL = 'https://github.com/LittleHannah/lncRNAPredModels/blob/master/'+model+'.tgz?raw=true'
    #print("downURL is ", downURL)
    modelAbsPath = saveDir + model 
    #print('modelAbsPath is ', modelAbsPath)
    downComm = 'wget -q -O ' + modelAbsPath + '.tgz ' + downURL
    if(os.path.exists(modelAbsPath)):
        print("----------------------------------Loading Completed-------------------------")
        return True
    #print('download command is ', downComm)
    os.system(downComm)
    extractComm = 'tar -xzvf ' + modelAbsPath + '.tgz'+' -C ' +saveDir
    #print('extract command is ' , extractComm)
    os.system(extractComm)
    rmComm = 'rm ' + modelAbsPath + '.tgz'
    #print('rm command is ' ,rmComm)
    os.system(rmComm)
    if(os.path.exists(modelAbsPath)):
        print("----------------------------------Loading Succeed-------------------------")
        return True
    print("----------------------------------Loading Failed-------------------------")
    return False
    
args = parse_args()
if args.manual and not args.model:
    os.system("cat README.md")
elif args.manual and args.model:
    os.system("cat docs/README_%s"%(args.model))
else:
    if args.model == "CPC2":
        CPC2.__main(args)
    elif args.model == "CPAT":
        CPAT.main(args)
    elif args.model == "GFStack":
        saveDir = modelDir + 'PredLnc_GFStack/src/'
        if(getModel('GFStack_mouse', saveDir)and getModel('GFStack_human', saveDir)):
            GFStack.main(args)
    elif args.model == "longdist":
        longdist.main(args)    
    elif args.model == "CPPred":
        CPPred.main(args)
    elif args.model == "CNCI":
        CNCI.main(args)
