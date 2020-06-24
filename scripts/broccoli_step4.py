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
import itertools
from multiprocessing import Pool as ThreadPool 
from pathlib import Path
from scripts import utils
try:
    from ete3 import PhyloTree
except:
    sys.exit("\n            ERROR: the ete3 library is not installed\n\n")



def step4_orthologous_pairs(lo, nsp, nt):

    # convert the parameters to global variables (horrible hack)
    global limit_ortho, not_same_sp, nb_threads
    limit_ortho, not_same_sp, nb_threads = lo, nsp, nt

    print('\n --- STEP 4: orthologous pairs\n')
    print(' ## parameters')
    print(' ratio ortho  : ' + str(limit_ortho))
    print(' not same sp  : ' + str(not_same_sp))
    print(' threads      : ' + str(nb_threads))
    print('\n ## load data')

    ## create output directory (delete it first if already exists)
    global out_dir
    out_dir = utils.create_out_dir('dir_step4')
       	
    ## check directory
    files_blast_list, files_tree = pre_checking(Path('dir_step2'))
    
    ## get original and combined names
    original_name = utils.get_pickle(Path('dir_step1') / 'original_names.pic')
    combined_prot = utils.get_pickle(Path('dir_step1') / 'combined_names.pic')
    
    ## load all data
    global all_no_tree, all_trees
    all_no_tree, all_trees, all_OGs = load_all_data(files_blast_list, files_tree)

    global prot_2_sp
    ## load prot_name 2 species dict (string version)
    prot_2_sp = utils.get_pickle(Path('dir_step2') / 'prot_str_2_species.pic')
    
    ## analyse OGs 1 by 1 and save ortho relationships in file
    print('\n ## analyse ' + str(len(all_OGs)) + ' orthologous groups 1 by 1')
    multithread_process_OG(all_OGs, nb_threads, original_name, combined_prot, not_same_sp)

    print(' done\n')
    

def pre_checking(directory):
       
    # check if directory exists
    if not directory.exists():
        sys.exit("\n            ERROR STEP 4: the directory dir_step2 does not exist.\n\n")
   
    # list the input _similarity_ortho.pic and _trees.pic pickle files 
    p = Path(directory / 'dict_similarity_ortho').glob('*')
    list_1 = [str(x.parts[-1]) for x in p if x.is_file() and '_similarity_ortho.pic' in str(x.parts[-1])] 
    p = Path(directory / 'dict_trees').glob('*')
    list_2 = [str(x.parts[-1]) for x in p if x.is_file() and '_trees.pic' in str(x.parts[-1])] 
    list_1.sort()
    list_2.sort()

    # print error
    if len(list_1) != len(list_2):
        sys.exit("\n            ERROR STEP 4: the number of *_trees.pic and *_blast_ortho.pic are different\n\n")
    return list_1, list_2


def load_all_data(f_blast, f_tree):
    print(' load NO tree results')  
    no_tree = utils.get_multi_pickle(Path('dir_step2') / 'dict_similarity_ortho', '_similarity_ortho.pic')
           
    print(' load tree results') 
    trees = utils.get_multi_pickle(Path('dir_step2') / 'dict_trees', '_trees.pic')
    
    print(' load OGs')    
    OGs = utils.get_pickle(Path('dir_step3') / 'OGs_in_network.pic')
    l_OGs = [l for l in OGs.values()]
    
    return no_tree, trees, l_OGs


def multithread_process_OG(l_ogs, n_threads, original, combined, not_same_sp):
    
    # start multithreading
    pool = ThreadPool(n_threads) 
    tmp_res = pool.map_async(process_OG, l_ogs, chunksize=1)
    results_2 = tmp_res.get()
    pool.close() 
    pool.join() 
    
    # save ortho relationships in file
    outfile = open(out_dir / 'orthologous_pairs.txt','w+')
    for s in results_2:
        if s != '':
            l = s.split(' ')
            for k in l:
                # extract pair of proteins
                prot1, prot2 = k.split('-')
                # get back combined proteins
                if int(prot1) in combined:
                    prot1_l = combined[int(prot1)]
                else:
                    prot1_l = [prot1]
                if int(prot2) in combined:
                    prot2_l = combined[int(prot2)]
                else:
                    prot2_l = [prot2]
                    
                # do nothing if option 'not_same_sp' AND the 2 proteins are from the same species
                if not_same_sp and prot_2_sp[prot1] == prot_2_sp[prot2]:
                    pass
                # build ortho relationships between all of them
                else:
                    for k1 in prot1_l:
                        for k2 in prot2_l:
                            outfile.write(original[int(k1)] + '	' + original[int(k2)] + '\n')


def process_OG(l_OG):
    
    # convert to set and protein names to string format (to match with pickle dict)
    s_OG = set(str(x) for x in l_OG)
    # prepare dict of ortho and para
    tmp_ortho = dict()
    tmp_para = dict()
    for prot in s_OG:
        tmp_ortho[prot] = collections.defaultdict(int)
        tmp_para[prot]  = collections.defaultdict(int)
        
    # get over each prot
    for prot in s_OG:
        # case prot in blast list -> all ortho
        if prot in all_no_tree:
            # remove prot that are not in OG
            l2 = [x for x in all_no_tree[prot] if x in s_OG]   
            # get all ortho relationships             
            for s in itertools.combinations(l2, 2):
                tmp_ortho[s[0]][s[1]] += 1
                tmp_ortho[s[1]][s[0]] += 1
        
        # case prot in trees -> extract orthos and paras
        elif prot in all_trees:
            '''
            this part of the code is a modified version of the ETE function 'get_evol_events_from_leaf'
            https://github.com/etetoolkit/ete/tree/master/ete3/phylo 
            '''
            ref_leaf  = prot
            tree = PhyloTree(all_trees[prot])

            # remove prot that are not in OG
            for leaf in tree:
                if leaf.name not in s_OG:
                    node = tree.search_nodes(name = leaf.name)[0] 
                    node.delete()
                    
            # browse the tree from terminal node
            node = tree.search_nodes(name = ref_leaf)[0]
            current  = node
            sister_leaves  = set([])
            browsed_spcs   = set(prot_2_sp[leaf.name] for leaf in current)
            browsed_leaves = set([current])
            while current.up:
                # distances control (0.0 distance check)
                for s in current.get_sisters():
                    for leaf in s.get_leaves():
                        sister_leaves.add(leaf)
                # Process sister node only if there is any new sequence (previene dupliaciones por nombres repetidos)
                sister_leaves = sister_leaves.difference(browsed_leaves)
                if len(sister_leaves)==0:
                    current = current.up
                    continue
                # Gets species at both sides of event
                sister_spcs        = set(prot_2_sp[leaf.name] for leaf in sister_leaves)
                overlaped_spces    = len(browsed_spcs & sister_spcs)
                all_spcs           = len(browsed_spcs | sister_spcs)
                sp_only_in_sister  = len(sister_spcs - browsed_spcs)
                sp_only_in_browsed = len(browsed_spcs - sister_spcs)
                
                names_browsed = set(n.name for n in browsed_leaves)
                if '' in names_browsed:
                    names_browsed = set(filter(None, names_browsed)) 
                names_sister = set(n.name for n in sister_leaves)
                if '' in names_sister:
                    names_sister = set(filter(None, names_sister))
                
                if all_spcs == 1 or overlaped_spces == 0 or (overlaped_spces == 1 and sp_only_in_sister >= 2 and sp_only_in_browsed >= 2):
                    # save ortho 
                    for prot1 in names_browsed:
                        for prot2 in names_sister:
                            tmp_ortho[prot1][prot2] += 1
                            tmp_ortho[prot2][prot1] += 1
                else:
                    # save para
                    for prot1 in names_browsed:
                        for prot2 in names_sister:
                            tmp_para[prot1][prot2] += 1
                            tmp_para[prot2][prot1] += 1
                        
                # Updates browsed species
                browsed_spcs   |= sister_spcs
                browsed_leaves |= sister_leaves
                sister_leaves  = set([])
                    
                # And keep ascending
                current = current.up
                
    # save ortho         
    final_ortho = ''       
    for prot1 in tmp_ortho:
        for prot2, nb_ortho in tmp_ortho[prot1].items():
            if prot1 > prot2:   # only consider 1 out of 2
                ratio = nb_ortho / (nb_ortho + tmp_para[prot1][prot2])
                if ratio > limit_ortho:
                    final_ortho += ' ' + prot1 + '-' + prot2
    
    return final_ortho.strip()


