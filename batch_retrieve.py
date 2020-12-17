import argparse
import os
import re

def retrieve_loc(uni):
    # get subcellular location
    temp_f = open('%s.txt'%(uni), 'r')
    location = ''
    for line in temp_f:
        s = re.search('SUBCELLULAR LOCATION:', line)
        if s is not None:
            locs = line[(s.span()[1]+1):len(line)].split('.')
            for i in range(len(locs)-1):
                location += locs[i].split(' {')[0]
                location += '~'
            location = location[:len(location)-1]
                
    
    return location

def retrieve_aa(uni):
    # get size of protein and amino acid composition


    aa = ''

    temp_f = open('%s.txt'%(uni), 'r')
    all_text = temp_f.readlines()
    for i in range(len(all_text)):
        all_text[i] = all_text[i].strip('\n')
        if all_text[i].startswith('SQ'):
            size = re.findall('([0-9]*) AA',all_text[i])[0]
            while True:
                i+=1
                if all_text[i].startswith(' '):
                    aa += all_text[i].replace(' ','')
                else:
                    break
            charged, polar, non_polar, hydropho = aa_com(aa)
            break

    return size, charged, polar, non_polar, hydropho


def aa_com(aa):
    # get the aa composition
    charged = 'RHKDE'
    polar = 'STNQ'
    non_polar = 'CGP'
    hydropho = 'AVILMFYW'

    charged_a = [aa.count(l) for l in charged]
    polar_a = [aa.count(l) for l in polar]
    non_polar_a = [aa.count(l) for l in non_polar]
    hydropho_a = [aa.count(l) for l in hydropho]

    return sum(charged_a)/len(aa), sum(polar_a)/len(aa), sum(non_polar_a)/len(aa), sum(hydropho_a)/len(aa)
    


parser = argparse.ArgumentParser()

parser.add_argument('--input', default=None, type=str, required=True, help='input file with uniprot ID')
parser.add_argument('--output', default=None, type=str, required=True, help='output file name')

args = parser.parse_args()

fw = open(args.output, 'w')
fw.write('\t'.join(['uniprot','location','size','charged','polar','non_polar','hydropho']))

with open(args.input, 'r') as f:
    for line in f:
        line = line.strip()
        os.system('wget https://www.uniprot.org/uniprot/%s.txt;'%(line))
        location = retrieve_loc(line)
        size, charged, polar, non_polar, hydropho = retrieve_aa(line)
        os.system('rm %s.txt;'%(line))
        fw.write('\t'.join([line, location, str(size), str(charged), str(polar), str(non_polar), str(hydropho)]))

fw.close()







