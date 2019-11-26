'''
    This file is part of Broccoli.

    Broccoli is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Broccoli is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Broccoli.  If not, see <https://www.gnu.org/licenses/>.
    
    contact: romain.derelle@gmail.com

'''


import argparse
import os
import sys
import shutil
from scripts import broccoli_step1
from scripts import broccoli_step2
from scripts import broccoli_step3
from scripts import broccoli_step4


######################### functions

def parse_args():
    # define and parse command-line arguments
    parser = argparse.ArgumentParser(description='            Broccoli v1.0', add_help=False, formatter_class=argparse.RawTextHelpFormatter, epilog=' \n')
    
    common = parser.add_argument_group(' general options')
    common.add_argument('-steps',         help='steps to be performed, comma separated (default = \'1,2,3,4\')', metavar='', type=str, default='1,2,3,4')    
    common.add_argument('-threads',       help='number of threads [default = 1]', metavar='', type=int, default=1)
    common.add_argument('-h','-help',     action="help", help="show this help message and exit")
    
    step1 = parser.add_argument_group(' STEP 1  kmer clustering')
    step1.add_argument('-dir',            help='name of the directory containing the proteome files [required]', metavar='')
    step1.add_argument('-ext',            help='extension of proteome files (default = \'.fasta\')', metavar='', type=str, default='.fasta')    
    step1.add_argument('-min_length',     help='minimum length of sequences [default = 10]', metavar='', type=int, default=10)
    step1.add_argument('-kmer_size',      help='length of kmers [default = 100]', metavar='', type=int, default=100)
    step1.add_argument('-kmer_min_aa',    help='minimum nb of different aa a kmer should have [default = 15]', metavar='', type=int, default=15)
    
    step2 = parser.add_argument_group(' STEP 2  phylomes')
    step2.add_argument('-path_diamond',   help='path of DIAMOND with filename [default = \'diamond\']', metavar='', type=str, default='diamond')    
    step2.add_argument('-path_fasttree',  help='path of FastTree with filename [default = \'fasttree\']', metavar='', type=str, default='fasttree')    
    step2.add_argument('-e_value',        help='e-value for similarity search [default = 0.001]', metavar='', type=float, default=0.001)
    step2.add_argument('-nb_hits',        help='maximum nb of hits per species [default = 6]', metavar='', type=int, default=6)
    step2.add_argument('-max_gap',        help='maximum fraction of gap per position [default = 0.7]', metavar='', type=float, default=0.7)
    
    step3 = parser.add_argument_group(' STEP 3  network analysis')
    step3.add_argument('-min_nb_hits',    help='minimum number of hits belonging to the OG [default = 2]', metavar='', type=int, default=2)
    step3.add_argument('-fusion_shared',  help='minimum fraction of connected nodes in each OG [default = 0.5]', metavar='', type=float, default=0.5)
    step3.add_argument('-fusion_nb_sp',   help='minimum nb of species in OGs involved in gene-fusions [default = 2]', metavar='', type=int, default=2)
    step3.add_argument('-fusion_overlap', help='maximum overlap between OGs in amino-acids [default = 10]', metavar='', type=int, default=10)

    step4 = parser.add_argument_group(' STEP 4  orthologous pairs')
    step4.add_argument('-ratio_ortho',    help='limit ratio ortho/total [default = 0.5]', metavar='', type=float, default=0.5)
    step4.add_argument('-not_same_sp',    help='ignore ortho relationships between proteins of the same species (QfO benchmark)', action="store_true")
    
    args = parser.parse_args()
    
    # clean directory name
    if args.dir:
        args.dir = clean_dir_name(args.dir)

    return args.steps, args.threads, \
    args.dir, args.ext, args.min_length, args.kmer_size, args.kmer_min_aa, \
    args.e_value, args.nb_hits, args.path_diamond, args.path_fasttree, args.max_gap, \
    args.min_nb_hits, args.fusion_shared, args.fusion_overlap, args.fusion_nb_sp, \
    args.ratio_ortho, args.not_same_sp



def clean_dir_name(d):
    d = d.replace('./','')
    d = d.replace('/','')
    d = './' + d + '/'
    return d.replace('//','/')


def check_python_version():
    if sys.version_info[0] != 3 and sys.version_info[1] < 6:
        sys.exit('\n            ERROR: your python is version '+ str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + ', please use version 3.6+\n\n')
  

def parse_steps(p):
    # split the steps and check them
    l = p.split(',')
    try:
        s = {int(x) for x in l}
    except:
        sys.exit('\n            ERROR: the list of steps should be composed of integers (-steps option)\n\n')
    # check for consecutiveness
    range_steps = max(s) - min(s)
    if range_steps != (len(s) - 1):
        sys.exit('\n            ERROR: the steps should be consecutive (-steps option)\n\n')
    return s


def pre_checking_pgms(p_diamond, p_fasttree):
    # check diamond
    if '/' in p_diamond:
        if not os.path.isfile(p_diamond):
            sys.exit("\n            ERROR: the path to DIAMOND is incorrect\n\n")
    elif not shutil.which(p_diamond):
        sys.exit("\n            ERROR: the path to DIAMOND is incorrect\n\n")

    # check FastTree
    if '/' in p_fasttree:
        if not os.path.isfile(p_fasttree):
            sys.exit("\n            ERROR: the path to FastTree is incorrect\n\n")
    elif not shutil.which(p_fasttree):
        sys.exit("\n            ERROR: the path to FastTree is incorrect\n\n")


######################### main part


if __name__ == "__main__":
    
    ## get all arguments
    pre_steps, nb_threads, \
    directory, extension, min_seq, length_kmer, min_aa, \
    evalue, max_per_species, path_diamond, path_fasttree, trim_thres, \
    min_nb_hits, limit_shared, limit_overlap, limit_nb_sp, \
    limit_ortho, not_same_sp = parse_args()


    print('\n            Broccoli v1.0\n')


    ## check python version
    check_python_version()
        
    ## parse steps
    steps = parse_steps(pre_steps)
    
    ## check if -dir option (cases of 1st step)
    if 1 in steps and directory is None:
        sys.exit('\n            ERROR: you need to specify an input directory (see -help)\n\n')
    
    ## check path of executables if step 2 (diamond, fasttree)
    if 2 in steps:
        pre_checking_pgms(path_diamond, path_fasttree)
        
    ## execute the steps
    if 1 in steps:
        broccoli_step1.step1_kmer_clustering(directory, extension, min_seq, length_kmer, min_aa, nb_threads)

    if 2 in steps:
        broccoli_step2.step2_phylomes(evalue, max_per_species, path_diamond, path_fasttree, trim_thres, nb_threads)
    
    if 3 in steps:
        broccoli_step3.step3_orthology_network(min_nb_hits, limit_shared, limit_overlap, limit_nb_sp, nb_threads)
    
    if 4 in steps:
        broccoli_step4.step4_orthologous_pairs(limit_ortho, not_same_sp, nb_threads)
    
    
    
    
    


