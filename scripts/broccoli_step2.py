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

import os
import sys
import collections
import csv
import subprocess
import gzip
import math
import re
from multiprocessing import Pool as ThreadPool 
from scripts import utils
try:
    from ete3 import PhyloTree
except:
    sys.exit("\n            ERROR: the ete3 library is not installed\n\n")




def step2_phylomes(eval, msp, pdia, pfas, tt, nt):

    # convert the parameters to global variables (horrible hack)
    global evalue, max_per_species, path_diamond, path_fasttree, trim_thres, nb_threads
    evalue, max_per_species, path_diamond, path_fasttree, trim_thres, nb_threads = eval, msp, pdia, pfas, tt, nt

    print('\n --- STEP 2: phylomes\n')
    print(' # parameters')
    print(' e_value  : ' + str(evalue))
    print(' nb_hit   : ' + str(max_per_species))
    print(' gap      : ' + str(trim_thres))
    print(' threads  : ' + str(nb_threads))
       	
    ## create output directory (or empty it if it already exists)
    utils.create_out_dir('./dir_step2')

    ## check directory input data
    global list_files
    print('\n # check input files')
    list_files = pre_checking_data('./dir_step1/')
    
    ## load all sequences
    global name_2_sp_phylip_seq, all_species
    name_2_sp_phylip_seq, all_species = create_dict_seq(list_files)

    ## create databases
    multithread_databases(list_files, nb_threads)
    
    ## process each proteome
    print('\n # build phylomes ... be patient')
    multithread_process_file(list_files, nb_threads)
    
    ## save prot 2 species dict (needed for steps 3 and 4)
    save_prot_2_sp(name_2_sp_phylip_seq)
    
    # clean directory from databases
    os.system('rm ./dir_step2/*.db.dmnd')
    
    print ('')

    
def pre_checking_data(directory):
       
    # check if directory exists
    if not os.path.isdir(directory):
        sys.exit("\n            ERROR STEP 2: the directory ./dir_step1 does not exist.\n\n")
     
    # list the input files inside that directory
    else: 
        list_files = list()
        for file in os.listdir(directory):
            if file.endswith('.fas'):
                list_files.append(file)
        # sort the list
        list_files.sort()
    
    # print message
    if len(list_files) == 0:
        sys.exit("\n            ERROR STEP 2: there is no input fasta file (*.fas) in the directory ./dir_step1/\n\n")
    else:
        print (' ' + str(len(list_files)) + ' input fasta files')
            
    return list_files


# --------------------- #

def create_dict_seq(l_files):
    d_seq = dict()
    d_sp  = dict()
    for file in l_files:
        # extract species index
        sp = int(file.split('.')[0])
        d_sp[sp] = 0
        # get sequences
        with open('./dir_step1/' + file) as fasta_content:
            for name_seq, seq in utils.read_fasta(fasta_content):
                # get name in phylip format
                phylip_name = name_seq + ' ' * (10 - len(name_seq))
                # save all within a tuple
                d_seq[name_seq] = (sp, phylip_name, seq)
    
    print (' ' + str(len(d_seq)) + ' sequences')
    return d_seq, d_sp


# --------------------- #

def multithread_databases(l_file, n_threads):
    # start multithreading
    pool = ThreadPool(n_threads) 
    tmp_res = pool.map_async(prepare_databases, l_file, chunksize=1)
    pool.close() 
    pool.join()    


def prepare_databases(file):
    subprocess.check_output(path_diamond + ' makedb --in ./dir_step1/' + file + ' --db ./dir_step2/' + file.replace('.fas','.db') + ' 2>&1', shell=True)


def multithread_process_file(l_file, n_threads):
    # start multithreading
    pool = ThreadPool(n_threads) 
    tmp_res = pool.map_async(process_file, l_file, chunksize=1)
    results_2 = tmp_res.get()
    pool.close() 
    pool.join()
    
    # load species dict
    dict_species = utils.get_pickle('./dir_step1/species_index.pic')
    
    # create log file
    log_file = open('./dir_step2/log_step2.txt', 'w+')
    log_file.write('#species_file	nb_phylo	nb_NO_phylo	nb_empty_ali_ali	nb_pbm_tree\n')
    
    # save log
    for l in results_2:
        log_file.write(dict_species[l[0]] + '	' + '	'.join(l[1:]) + '\n')
    log_file.close()


def analyse_species(dict_sp):
    present = 0
    dupli = 0
    for k in dict_sp.values():
        if k > 0:
            present += 1
        if k > 1:
            dupli += 1
    return present, dupli


''' new function that create aligned HSP sequences from cigar string '''
def extract_HSP(full_seq, start, cig):
    # split cigar 
    matches = re.findall(r'(\d+)([A-Z]{1})', cig)
    l_tup = [(m[1], int(m[0])) for m in matches]

    # reconstruct HSP
    hsp = ''
    position = start
    for t in l_tup:
        # case of aa matches
        if t[0] == 'M':
            hsp += full_seq[position:(position + t[1])]
            position += t[1]
        # case of deletion
        elif t[0] == 'D':
            position += t[1]
        # case of insertion
        elif t[0] == 'I':
            hsp += '-' * t[1]
    return hsp


def process_location(qu_start, qu_end, min_start, max_end):
    if qu_start < min_start:
        min_start = qu_start
    if qu_end > max_end:
        max_end = qu_end
    return min_start, max_end


def get_positions(ref_name, hits, t):
    # count number of '-' per position
    tmp_ = [0] * len(hits[ref_name])
    for seq in hits.values():
        for n in range(0, len(seq)):
            if seq[n] == '-':
                tmp_[n] += 1
    # select good positions
    nb_seq = len(hits)
    good = set()
    for i,v in enumerate(tmp_):
        if (v / nb_seq) < trim_thres:
            good.add(i)
    return good


def process_file(file):       
    ## extract index
    index = file.split('.')[0]
    
    ## create output directory (or empty it if it already exists)
    if not os.path.isdir('./dir_step2/' + index):
        os.mkdir('./dir_step2/' + index)
    else:
        os.system('rm ./dir_step2/' + index + '/*.gz')
    
    ## perform local search against each database
    for file_db in list_files:
        search_output = index + '_' + file_db.replace('.fas','.gz')
        subprocess.check_output(path_diamond + ' blastp --quiet --threads 1 --db ./dir_step2/' + file_db.replace('.fas','.db') + ' --max-target-seqs ' + str(max_per_species) + ' --query ./dir_step1/' + file + ' \
                 --compress 1 --more-sensitive -e ' + str(evalue) + ' -o ./dir_step2/' + index + '/' + search_output + ' --outfmt 6 qseqid sseqid qstart qend sstart cigar 2>&1', shell=True)
    
    ## get all hits in a dict of list
    all_output = collections.defaultdict(list)
    for f in os.listdir('./dir_step2/' + index):
        if f.endswith('.gz'):
            with gzip.open('./dir_step2/' + index + '/' + f, mode="rt") as f:
                file_content = csv.reader(f, delimiter='	')
                for line in file_content:
                    # save output
                    all_output[line[0]].append(line[1:])     
       
    ## analyse BLAST hits
    nb_phylo    = 0
    nb_NO_phylo = 0
    nb_empty_ali  = 0
    all_alis = dict()
    no_phylo = dict()
    
    for prot in all_output:
        min_start = math.inf
        max_end   = 0
        ref_species = dict(all_species)
        all_hits = dict()
        targets_done = set()
        
        # get all hits for this prot
        reduced = list()
        for ll in all_output[prot]:
            target = ll[0]
            species = name_2_sp_phylip_seq[target][0]
            qu_start = int(ll[1]) -1
            qu_end   = int(ll[2]) -1
            ta_start = int(ll[3]) -1
            cigar    = ll[4]
            
            if target in all_hits:
                # extract HSP and add to target seq
                HSP = extract_HSP(name_2_sp_phylip_seq[target][2], ta_start, cigar)
                all_hits[target] = all_hits[target][:qu_start] + HSP + all_hits[target][qu_end + 1:]
                min_start, max_end = process_location(qu_start, qu_end, min_start, max_end) 
                
            else:
                ref_species[species] += 1
                # create target seq
                all_hits[target] = '-' * len(name_2_sp_phylip_seq[prot][2])
                # extract HSP
                HSP = extract_HSP(name_2_sp_phylip_seq[target][2], ta_start, cigar)
                all_hits[target] = all_hits[target][:qu_start] + HSP + all_hits[target][qu_end + 1:]
                min_start, max_end = process_location(qu_start, qu_end, min_start, max_end)              
            
            # reduce output for pickle (convert all element to integers)
            reduced.append(tuple(int(x) for x in ll[:3]))    
        
        # save reduced output
        all_output[prot] = tuple(reduced)
        
        # add query to hits if not there (it happens sometimes when many similar sequences from the same species)
        if prot not in all_hits:
            all_hits[prot] = name_2_sp_phylip_seq[prot][2]
        
        # analyse species content
        nb_present, nb_dupli = analyse_species(ref_species)
        # case phylogenetic analysis
        if nb_present > 1 and nb_dupli > 0:
            nb_phylo += 1
            
            # find good positions
            good_positions = get_positions(prot, all_hits, trim_thres)
            # save alignment
            if len(good_positions) == 0:
                nb_empty_ali += 1
            elif len(good_positions) < 4988:  ## FastTtree2 limitation (5000 per line)
                new_ali = [str(len(all_hits)) + '	' + str(len(good_positions))]
                for name, seq in all_hits.items():
                    trimed_seq = [seq[n] for n in range(len(seq)) if n in good_positions]
                    new_ali.append(name_2_sp_phylip_seq[name][1] + ''.join(trimed_seq))
                all_alis[prot] = '\n'.join(new_ali)
            else:
                # take only the 4988 first positions (longer alignments are very rare anyway)
                new_ali = [str(len(all_hits)) + '	4988']
                for name, seq in all_hits.items():
                    trimed_seq = [seq[n] for n in range(len(seq)) if n in good_positions][:4988]
                    new_ali.append(name_2_sp_phylip_seq[name][1] + ''.join(trimed_seq)) 
                all_alis[prot] = '\n'.join(new_ali)
                       
        # case NO phylogenetic analysis
        else:
            nb_NO_phylo += 1
            # sort the list of names for further processing
            xx = list(all_hits)
            xx.sort()
            no_phylo[prot] = xx

    
    ## convert DIAMOND output (keys to integers) and save it to file
    all_output = {int(x):t for x,t in all_output.items()}
    utils.save_pickle('./dir_step2/' + index + '_output.pic', all_output)

    ## save set-aside groups to file
    utils.save_pickle('./dir_step2/' + index + '_blast_ortho.pic', no_phylo)
    
    ## save all alignments to file
    alis_file = open('./dir_step2/alis_' + index + '.phy', 'w+')
    all_ref_prot = list()
    for ref_prot, ali in all_alis.items():
        alis_file.write(ali + '\n')
        all_ref_prot.append(ref_prot)
    alis_file.close()
    
    # free memory
    nb_alis = len(all_alis)
    all_alis   = None
    all_output = None
    
    ## perform phylogenetic analyses and root trees
    all_trees  = dict()
    nb_pbm_tree = 0
    a = subprocess.check_output(path_fasttree + ' -quiet -nosupport -noml -nome -fastest -bionj -pseudo -n ' + str(nb_alis) + ' ./dir_step2/alis_' + index + '.phy 2>&1', shell=True)
    a2 = a.strip().decode("utf-8")
    a3 = a2.split('\n')
    c = -1
    for line in a3:
        # case the line is in the form 'Ignored unknown character ...'
        if line.startswith('Ign'):
            pass
        else:
            c += 1
            if not line.startswith('('):
                nb_pbm_tree += 1            
                # security
                if nb_pbm_tree > 100:
                    sys.exit("\n            ERROR STEP 2: too many errors in phylogenetic analyses -> stopped\n\n")
            else:
                # import tree in ete3 and root it
                ete_tree = PhyloTree(line)
                mid = ete_tree.get_midpoint_outgroup()
                try:
                    ete_tree.set_outgroup(mid)
                except:
                    pass
                # get reference protein name
                prot = all_ref_prot[c]
                # save rooted tree
                all_trees[prot] = ete_tree.write()
    
    ## save trees to file
    utils.save_pickle('./dir_step2/' + index + '_trees.pic', all_trees)
    
    # clean directory
    os.system('rm ./dir_step2/alis_' + index + '.phy')
    
    return [index, str(nb_phylo), str(nb_NO_phylo), str(nb_empty_ali), str(nb_pbm_tree)]

      
def save_prot_2_sp(d):
    # create new dict with only sp
    d_str = {k:t[0] for k,t in d.items()} 
    d_int = {int(k):t[0] for k,t in d.items()} 
    # save them on disk
    utils.save_pickle('./dir_step2/prot_str_2_species.pic', d_str)
    utils.save_pickle('./dir_step2/prot_int_2_species.pic', d_int)
   
