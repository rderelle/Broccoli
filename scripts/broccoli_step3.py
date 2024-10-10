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
import operator
import itertools
import pickle
import statistics
from multiprocessing import Pool 
from pathlib import Path
import shutil
import gc
from scripts import utils
try:
    from ete3 import PhyloTree
except:
    sys.exit("\n            ERROR: the ete3 library is not installed\n\n")



def step3_orthology_network(rov, mw, mnh, lm, nbsp, nt):

    # convert the parameters to global variables (horrible hack)
    global sp_overlap, min_weight, min_nb_hits, chimeric_edges, chimeric_species, nb_threads
    sp_overlap, min_weight, min_nb_hits, chimeric_edges, chimeric_species, nb_threads = rov, mw, mnh, lm, nbsp, nt
    
    print('\n --- STEP 3: network analysis\n')
    print(' ## parameters')
    print(' species overlap  : ' + str(sp_overlap))
    print(' min edge weight  : ' + str(min_weight))
    print(' min nb hits      : ' + str(min_nb_hits))
    print(' chimeric edges   : ' + str(chimeric_edges))
    print(' chimeric species : ' + str(chimeric_species))
    print(' threads          : ' + str(nb_threads))

    ## create output directory (or empty it if it already exists)
    global out_dir
    out_dir = utils.create_out_dir('dir_step3')
        
    ## get all species
    all_species = utils.get_pickle(Path('dir_step1') / 'species_index.pic')
    
    ## check directory
    files_trees = pre_checking(Path('dir_step2'))
    
    ## create log file
    global log_file
    log_file = open(out_dir / 'log_step3.txt', 'w+')
        
    global prot_2_sp
    ## load prot_name 2 species dict (string version)
    prot_2_sp = utils.get_pickle(Path('dir_step2') / 'prot_str_2_species.pic')
    
    print('\n ## get ortho and para')
    ## extract ortho from dir_step2
    extract_ortho(files_trees)
    
    ## extract para from dir_step2
    prot_2_sp = None
    extract_para(files_trees)
    
    ## load prot_name 2 species dict (integer version)
    prot_2_sp = utils.get_pickle(Path('dir_step2') / 'prot_int_2_species.pic')
    
    print('\n ## network analysis')
    ## build network
    global all_nodes, all_edges 
    all_nodes, all_edges = build_network()
    
    ## load all search outputs
    print(' load similarity search outputs')
    global other_hits
    other_hits = load_search_outputs(Path('dir_step2') / 'dict_output', '_output.pic')

    ## define maximum number of edges to consider for lcc (= 2 * nb_species or 10 if less species) -> improve speed
    global limit_degree, limit_nb_max
    limit_degree, limit_nb_max = get_limit_lcc(all_species)

    ## calculate the local clustering coefficient of each node
    global all_lcc
    all_lcc = multithread_lcc(all_nodes, nb_threads)

    ## get all connected_components
    list_cc = utils.get_connected_components(all_edges, all_nodes)
    
    ## analyse each connected_components
    communities, all_chimeric_prot = analyse_cc(list_cc)
    
    ## remove spurious hits
    cleaned_communities = remove_spurious_hits(communities)
    
    ## save OG lists, fusions and stats
    save_outputs(cleaned_communities, all_chimeric_prot, all_species) 
    
    ## clean dir
    Path.unlink(out_dir / 'd_ortho.pic')
    Path.unlink(out_dir / 'd_para.pic')
    Path.unlink(out_dir / 's_filter.pic')
    
    print('')   
    
   
def pre_checking(directory):
       
    # check if directory exists
    if not directory.exists():
        sys.exit("\n            ERROR STEP 3: the directory dir_step2 does not exist.\n\n")
   
    # list the input _similarity_ortho.pic and _trees.pic pickle files (only names of _trees.pic files are returned)
    p = Path(directory / 'dict_similarity_ortho').glob('*')
    list_1 = [str(x.parts[-1]) for x in p if x.is_file() and '_similarity_ortho.pic' in str(x.parts[-1])] 
    p = Path(directory / 'dict_trees').glob('*')
    list_2 = [str(x.parts[-1]) for x in p if x.is_file() and '_trees.pic' in str(x.parts[-1])] 
    list_2.sort()
    
    # print error
    if len(list_1) != len(list_2):
        sys.exit("\n            ERROR STEP 3: the number of *_trees.pic and *_blast_ortho.pic are different\n\n")
    return list_2

# --------------------- #

def extract_ortho(l_trees):
        
    print(' extract ortho from similarity')
    # load pickle files
    tmp_d = utils.get_multi_pickle(Path('dir_step2') / 'dict_similarity_ortho', '_similarity_ortho.pic')

    # extract ortho
    d_ortho = collections.defaultdict(int)
    for l in tmp_d.values():
        for sub in itertools.combinations(l, 2):
            pair_int = int(str(len(sub[0])) + sub[0] + sub[1])
            d_ortho[pair_int] += 1
    
    global path_tmp, path_tmp_ortho
    path_tmp = utils.create_out_dir('./dir_step3/tmp')
    path_tmp_ortho = utils.create_out_dir('./dir_step3/tmp_ortho')
    
    # extract ortho from tree files
    print(' extract ortho from trees')    
    files_start = zip(l_trees, itertools.repeat(path_tmp), itertools.repeat(path_tmp_ortho), itertools.repeat(prot_2_sp), itertools.repeat(sp_overlap))    
    pool = Pool(nb_threads) 
    tmp_res = pool.starmap_async(extract_ortho_from_trees, files_start, chunksize=1)
    pool.close() 
    pool.join()            
    
    # unpack ortho and save them
    for filename in l_trees:
        content_pickle = pickle.load(open(path_tmp_ortho / filename, 'rb'))
        for pair in content_pickle:
            d_ortho[pair] += 1 

    # free memory
    content_pickle = None
    shutil.rmtree(path_tmp_ortho)

    # remove ortho found only once
    print(' remove ortho found only once')
    d_ortho = {k:v for k,v in d_ortho.items() if v > 1}
    
    # save it to file
    utils.save_pickle(out_dir / 'd_ortho.pic', d_ortho)

    # save a simplified version (without 2 first digits) as set to file
    s_filter = set()
    for k in d_ortho:
        k2 = str(k)
        s_filter.add( int(k2[2:]) )
    utils.save_pickle(out_dir / 's_filter.pic', s_filter)
    

def extract_ortho_from_trees(filename, path_tmp, path_tmp_ortho, prot_2_sp, sp_overlap):
    # prepare output variables
    l_ortho = list()
    l_ortho_para = list()
    
    # load dict of trees
    tmp_d = utils.get_pickle(Path('dir_step2') / 'dict_trees' / filename)
    
    # analyse trees 1 by 1
    for ref_leaf, newick in tmp_d.items():
        
        # load tree and get all leaves
        tree = PhyloTree(newick)
        all_leaves = {leaf.name for leaf in tree}
                                   
        # get all leaves from last interesting nodes      
        ref_node = tree.search_nodes(name = ref_leaf)[0]
        ortho = custom_species_overlap(ref_node, prot_2_sp, sp_overlap)
        
        # add ref_leaf to ortho in case no good node selected
        if len(ortho) == 0:
            ortho.add(ref_leaf)
        
        # get para
        para = all_leaves - ortho
    
        # save ortho
        xx = list(ortho)    
        xx.sort()    
        for sub in itertools.combinations(xx, 2):
            pair_int = int(str(len(sub[0])) + sub[0] + sub[1])
            l_ortho.append(pair_int)
        
        # save ortho@para if there is a paralogous group
        if para:
            l_ortho_para.append(' '.join(ortho) + '@' + ' '.join(para))

    # save ortho @ para
    utils.save_pickle(path_tmp / filename, l_ortho_para)
    
    # save ortho
    utils.save_pickle(path_tmp_ortho / filename, l_ortho)
    
    return [0,0]


def custom_species_overlap(node, prot_2_sp, sp_overlap):
    '''
    this function is a modified version of the ETE function 'get_evol_events_from_leaf'
    https://github.com/etetoolkit/ete/tree/master/ete3/phylo 
    '''
    # Prepare to browse tree from leaf to root
    last_good = set()
    current  = node
    sister_leaves  = set([])
    #browsed_spcs   = set([current.species])
    browsed_spcs   = set(prot_2_sp[leaf.name] for leaf in current)
    browsed_leaves = set([current])

    while current.up:
        # distances control (0.0 distance check)
        for s in current.get_sisters():
            for leaf in s.get_leaves():
                sister_leaves.add(leaf)
        # Process sister node only if there is any new sequence (previene dupliaciones por nombres repetidos)
        sister_leaves = sister_leaves.difference(browsed_leaves)
        if len(sister_leaves) == 0:
            current = current.up
            continue
        # Gets species at both sides of event
        sister_spcs        = set(prot_2_sp[leaf.name] for leaf in sister_leaves)
        overlaped_spces    = len(browsed_spcs & sister_spcs)
        all_spcs           = len(browsed_spcs | sister_spcs)

        ratio_sister   = overlaped_spces / len(sister_spcs)
        ratio_browsed  = overlaped_spces / len(browsed_spcs)

        # Updates browsed species
        browsed_spcs   |= sister_spcs
        browsed_leaves |= sister_leaves
        sister_leaves  = set([])
        
        if all_spcs == 1 or (ratio_sister <= sp_overlap and ratio_browsed <= sp_overlap):
            last_good = set(n.name for n in browsed_leaves)
        
        # And keep ascending
        current = current.up
    return last_good


def extract_para(l_trees):
    global set_filter
    set_filter = pickle.load(open(out_dir / 's_filter.pic', 'rb'))
    
    global path_tmp_para
    path_tmp_para  = utils.create_out_dir('./dir_step3/tmp_para')
    
    # extract para from tree files
    print(' extract para from trees')   
    files_start = zip(l_trees, itertools.repeat(set_filter), itertools.repeat(path_tmp_para), itertools.repeat(path_tmp))    
    pool = Pool(nb_threads) 
    tmp_res = pool.starmap_async(extract_para_from_trees, files_start, chunksize=1)
    pool.close() 
    pool.join()

    # combine results list para
    d_para = utils.get_pickle(out_dir / 'd_ortho.pic')
    d_para = {x:0 for x in d_para}
    for filename in l_trees:
        # load pick file with list para
        content_pickle = pickle.load(open(path_tmp_para / filename, 'rb'))
        for pair in content_pickle:
            try:
                d_para[pair] += 1 
            except:
                pass
                
    # save it to file
    utils.save_pickle(out_dir / 'd_para.pic', d_para)
    
    # free memory
    set_filter = None
    shutil.rmtree(path_tmp_para)
    shutil.rmtree(path_tmp)


def extract_para_from_trees(filename, set_filter, path_tmp_para, path_tmp):
    
    # prepare list
    out_para = list()

    # load filename and save ortho@para
    tmp = open(path_tmp / filename, 'rb')
    content_pickle = pickle.load(tmp)

    # save para
    for st in content_pickle:
        ortho, para = st.split('@')
        l_ortho = ortho.split(' ')
        l_para = para.split(' ')
        for name1 in l_para:
            for name2 in l_ortho:
                if name1 < name2:
                    combined_name = str(len(name1)) + name1 + name2
                else:
                    combined_name = str(len(name2)) + name2 + name1
                
                if len(combined_name) > 2:
                # check if present in filter and save it
                    reduced = int(combined_name[2:])
                    if reduced in set_filter:
                        out_para.append(int(combined_name))
    
    # dump out_para 
    utils.save_pickle(path_tmp_para / filename, out_para)

    return [0,0]


# --------------------- #

def build_network():
     
    d_edges = collections.defaultdict(dict)
    d_nodes = dict()
    nb_edges = 0
    
    d_ortho_pairs = utils.get_pickle(out_dir / 'd_ortho.pic')
    d_para_pairs  = utils.get_pickle(out_dir / 'd_para.pic')
    
    print(' build network:')
    for k, nb_ortho in d_ortho_pairs.items():
        # build edge and nodes if nb ortho superior to nb para
        if nb_ortho > d_para_pairs[k]:
            # extract names and convert them to integers
            s_k = str(k)
            size_1st_int = int(s_k[0]) + 1
            name1, name2 = int(s_k[1:size_1st_int]), int(s_k[size_1st_int:])
            # build edges and nodes
            d_edges[name1][name2] = nb_ortho
            d_edges[name2][name1] = nb_ortho
            d_nodes[name1] = 0
            d_nodes[name2] = 0
            nb_edges += 1
    
    # print results
    print('      _ ' + str(len(d_nodes)) + ' nodes')
    print('      _ ' + str(nb_edges) + ' edges')
    
    # save to log file
    log_file.write('#network size:\n' + str(len(d_nodes)) + ' nodes\n' + str(nb_edges) + ' edges\n\n')
        
    # get maximum connected value for each node
    max_ortho = dict()
    for node in d_nodes:
        max_ortho[node] = max(d_edges[node].items(), key=operator.itemgetter(1))[1]
        
    # convert values to weight ratio nb_ortho / max_ortho
    log_file.write('\n#edge_weight	nb_edges\n')
    vector_weights = [0] * 21
    nb_removed = 0
    for node1 in d_edges:
        to_remove = set()
        for node2, nb_ortho in d_edges[node1].items():
            weight = nb_ortho / max_ortho[node1]
            # save weight in vector for log file
            r_weight = int(round(20 * weight))
            try:
                vector_weights[r_weight] += 1
            except:
                print(r_weight)
            # decide to remove edge or not based on the minimum edge weight
            if weight < min_weight:
                to_remove.add(node2)
            else:
                d_edges[node1][node2] = nb_ortho / max_ortho[node1]
        # remove weak edges
        for node2 in to_remove:
            del d_edges[node1][node2]
            del d_edges[node2][node1]
            nb_removed += 1
    # save weights distribution in log file
    for i,v in enumerate(vector_weights):
        log_file.write(str(i/20) + '	' + str(v) + '\n')
       
    log_file.write('\n-> ' + str(nb_removed) + ' edges removed\n\n')
    
    # sort all edges of each node by values (used to break ties in the label propagation)
    for node1 in d_edges:
        d_sorted = dict()
        for key, value in sorted(d_edges[node1].items(), key=lambda x: x[1], reverse=True):
            d_sorted[key] = value
        # save it in the form of tuple
        d_edges[node1] = d_sorted    
    
    return d_nodes, d_edges


def multithread_lcc(d_nodes, n_threads):

    # define number of threads (limit to 3 threads to avoid high memory consumption)
    if n_threads > 3:
        nb_thr = 3
    else:
        nb_thr = n_threads

    # split node dict in a list of lists
    list_nodes = [x for x in d_nodes]
    new_list_of_lists = [list_nodes[i::nb_thr] for i in range(nb_thr)]  
        
    # start multithreading
    print(' compute lcc for each node')
    files_start = zip(new_list_of_lists, itertools.repeat(all_edges), itertools.repeat(limit_degree))    
    pool = Pool(nb_thr) 
    tmp_res = pool.starmap_async(calculate_lcc, files_start, chunksize=1)
    results_2 = tmp_res.get()
    pool.close() 
    pool.join()   
    
    # get all results together
    d_lcc = dict()
    for l in results_2:
        for t in l:
            d_lcc[t[0]] = t[1]

    return  d_lcc 


def get_limit_lcc(d_species):
    ld = 2 * len(d_species)
    nb_max = ld * (ld - 1)
    if ld < 10:
        ld = 10
        nb_max = 90
    return ld, nb_max


def calculate_lcc(l, all_edges, limit_degree):
    out = list()
    for k in l:
        degree = len(all_edges[k])
        if degree < 4:
            lcc = 0    
        else:
            if degree > limit_degree:
                nb_max = limit_nb_max
                tmp_l = sorted([x for x in all_edges[k]])[:limit_degree]
            else:
                nb_max = degree * (degree - 1)
                tmp_l = [x for x in all_edges[k]]
            
            nb_found = 0
            for i, neigb_1 in enumerate(tmp_l):
                for n in range(i+1, len(tmp_l)):
                    neigb_2 = tmp_l[n]
                    if neigb_1 in all_edges[neigb_2]:
                        nb_found += 1
            lcc = round( 2 * nb_found / nb_max, 5)
        out.append((k,lcc))
    return out


def load_search_outputs(dir_, str_):
    
    # dict of number of hits (used for detection of spurious hits)
    d_other_hits = collections.defaultdict(set)
    
    # list pickle files
    p = dir_.glob('*')
    tmp_l = [x for x in p if x.is_file() and str_ in str(x.parts[-1])]
    # load pickle files
    for file_path in tmp_l:
        d = dict()
        with open(file_path, 'rb') as content:
            gc.disable()  # disable garbage collector
            d.update(pickle.load(content))
            gc.enable()   # enable garbage collector again
        
        # add output search to edges
        for query, t in d.items():
            for t2 in t:
                if t2[0] != query:
                    # add info to network if hit is present in the network
                    if t2[0] in all_edges[query]:
                        # first time we add this target
                        if type(all_edges[query][t2[0]]) == float: 
                            all_edges[query][t2[0]] = (all_edges[query][t2[0]], t2[1], t2[2])
                        else:
                            start = min([all_edges[query][t2[0]][1], t2[1]])
                            end   = max([all_edges[query][t2[0]][2], t2[2]])
                            # save it
                            all_edges[query][t2[0]] = (all_edges[query][t2[0]][0], start, end)
                    # otherwise add it to other hits
                    else:
                        d_other_hits[query].add(t2[0])
        
    # check all edges
    for node1, d in all_edges.items():
        for node2, x in d.items():
            if type(x) == float:
                all_edges[node1][node2] = (x, 0, 0) 
    
    return d_other_hits           


# --------------------- #

def analyse_cc(l_cc):
    
    print(' apply LPA and corrections:')    
    final_list  = list()
    d_chimeric  = dict()
    
    for ll in l_cc:
        # check number of and species nodes in cc -> apply LPA method if nb_species > 1 and nb_nodes > 4
        nb_species = fast_count_species(ll)
        nb_nodes   = len(ll)

        if nb_nodes < 4 or nb_species == 1:
            # save list
            final_list.append(ll)
                
        else:
            
            nodes_in_cc = {k:all_lcc[k] for k in ll}
            
            # sort the node dict by lcc (decreasing order) + replace lcc value by label
            tmp_nodes_in_cc = dict()
            for key, value in sorted(nodes_in_cc.items(), key=lambda x: x[1], reverse=True):
                tmp_nodes_in_cc[key] = key
        
            # apply LPA to identify communities
            node2label = label_propagation(tmp_nodes_in_cc)
                             
            # reconstruct communities (i.e. labels present in values of node2label)
            label2node = create_label_dict(node2label)
                   
            # build list of communities (list of list of nodes)
            tmp_com = [l for l in label2node.values()]
            
            if len(tmp_com) == 1:
                # save community  
                for l in tmp_com:
                    final_list.append(l)
                
            else:
                # detect gene-fusions and modify communities accordingly
                tmp2_com, chimeric = detect_chimeric_proteins(tmp_com)
                
                # save community  
                for l in tmp2_com:
                    final_list.append(l)
                # save chimeric proteins
                for k in chimeric:
                    d_chimeric[k] = list()
    
    # print results
    print('      _ ' + str(len(l_cc)) + ' connected components')
    print('      _ ' + str(len(final_list)) + ' communities')
    print('      _ ' + str(len(d_chimeric)) + ' chimeric proteins')    
    
    # save to log file
    log_file.write('#network analysis:\n' + str(len(l_cc)) + ' connected components\n' + str(len(final_list)) + ' communities\n' + str(len(d_chimeric)) + ' chimeric proteins')
    
    return final_list, d_chimeric

              
def label_propagation(tmp_nodes_in_cc):

    # loop until no difference in node labels between 2 generations
    modif = True
    while modif:
        
        modif = False 
        # prepare previous label dict()
        previous_labels = dict(tmp_nodes_in_cc)
        
        # start label propagation ... nodes 1 by 1 (highest lcc first)
        for node in tmp_nodes_in_cc:
            
            # get all (ortho-para) values for each label
            all_neighbors_labels = collections.defaultdict(float)
                
            for v in all_edges[node]:
                label = tmp_nodes_in_cc[v]
                all_neighbors_labels[label] += all_edges[node][v][0]
                        
            new_label = max(all_neighbors_labels.items(), key=operator.itemgetter(1))[0]            
            
            # update its label with new_label
            if tmp_nodes_in_cc[node] != new_label:
                tmp_nodes_in_cc[node] = new_label
                modif = True  
        
    return tmp_nodes_in_cc
    
                   
def create_label_dict(in_d):
    out_d = dict()
    for k, v in in_d.items():
        out_d.setdefault(v, []).append(k)  
    return out_d

     
def detect_chimeric_proteins(ll_com):

    s_chimeric = set()

    # prepare dict limit nodes per OG (nb nodes * limit ratio) AND node 2 OG id
    limit_nodes_per_OG = dict()
    node_2_OG          = dict()
    for i,l in enumerate(ll_com):
        limit_nodes_per_OG[i] = len(l) * chimeric_edges
        for k in l:
             node_2_OG[k] = i
    
    # each node one by 1
    for node1 in node_2_OG:
        # count nb of connected OGs (other than its own OG)
        d_found  = collections.defaultdict(list)
        for node2 in all_edges[node1]:
            # node2 in node_2_OG and 
            if node2 in node_2_OG and all_edges[node1][node2][2] != 0:
                og_2 = node_2_OG[node2]
                d_found[og_2].append(node2)
                            
        if len(d_found) > 1:
            
            # test nb of species in connected OGs
            l_found = [i for i in d_found]
            for i in l_found:
                nb_sp = len({prot_2_sp[x] for x in d_found[i]})
                nb_edges = len(d_found[i])
                if nb_sp < chimeric_species or nb_edges < limit_nodes_per_OG[i]:
                    del d_found[i]
                        
            if len(d_found) > 1 and node_2_OG[node1] in d_found:    # test if OG of that node is still there
                
                # get start end of its own OG
                all_start = [all_edges[node1][x][1] for x in d_found[node_2_OG[node1]] if all_edges[node1][x][2] != 0]
                ref_start = statistics.median(all_start)
                all_end   = [all_edges[node1][x][2] for x in d_found[node_2_OG[node1]] if all_edges[node1][x][2] != 0]
                ref_end   = statistics.median(all_end)
                            
                # remove its own OG
                del d_found[node_2_OG[node1]]
                
                # test overlap with each OG
                for i, l in d_found.items():                    
                    all_start  = [all_edges[node1][x][1] for x in l if all_edges[node1][x][2] != 0]
                    test_start = statistics.median(all_start)
                    all_end    = [all_edges[node1][x][2] for x in l if all_edges[node1][x][2] != 0]
                    test_end   = statistics.median(all_end)
                                         
                    if test_start > ref_start:
                        overlap = ref_end - test_start
                    else:
                        overlap = test_end - ref_start
                    
                    if overlap <= 0:
                        # save it
                        s_chimeric.add(node1)
                        # add it to the other community
                        ll_com[i].append(node1)
        
    return ll_com, s_chimeric


def fast_count_species(l):
    s = set()
    for name in l:
        s.add(prot_2_sp[name])
        if len(s) > 1:
            break
    return len(s)
    

def count_species(l):
    s = {prot_2_sp[name] for name in l}
    return len(s)


def remove_spurious_hits(l_com):
    nb_removed = 0
    l2_com = list()
   
    for com in l_com:
        if len(com) > (min_nb_hits + 1):
            ref_s = set(com)
            # check all nodes
            new_com = list()
            for node in com:
                nb_hits = 0
                # check in network
                if node in all_edges:
                    for k in all_edges[node]:
                        if all_edges[node][k][2] != 0 and k in ref_s:
                            nb_hits += 1
                # check in other hits
                if node in other_hits:
                    for k in other_hits[node]:
                        if k in ref_s:
                            nb_hits += 1
                # verdict
                if nb_hits < min_nb_hits:
                    nb_removed += 1
                else:
                    new_com.append(node)
            # save new_com
            l2_com.append(new_com) 
        else:
            l2_com.append(com)
            
    print('      _ ' + str(nb_removed) + ' spurious hits removed') 
    return l2_com


def save_outputs(l_com, d_chimeric, d_species):    
    
    ## load original and combined names
    original_name = utils.get_pickle(Path('dir_step1') / 'original_names.pic')
    combined_prot = utils.get_pickle(Path('dir_step1') / 'combined_names.pic')
           
    # create vector nb_species as index and nb_OG as value
    vector_sp = [0] * (len(d_species) + 1)
    
    # create dict for species counts
    nb_per_sp = collections.defaultdict(int)

    # create dict for table OGs
    table_og = dict()
    
    ## save the lists of OGs and get fusion info and get OG info
    classified = set()
    d_OG = collections.defaultdict(list)
    c = 0
    OGs_in_network = dict()
    file_list_OGs = open(out_dir / 'orthologous_groups.txt', 'w+')
    file_list_OGs.write('#OG_name	protein_names\n')
    for l in l_com:
        nb_species = count_species(l)
        # keep OG if more than 1 species
        if nb_species > 1:
            vector_sp[nb_species] += 1
            # create name OG
            c += 1
            name_OG = 'OG_' + str(c)
            # create OG vector
            #table_og[name_OG] = {x:0 for x in d_species}
            table_og[name_OG] = {x:list() for x in d_species}
            # save old names
            OGs_in_network[name_OG] = l
            # prepare full OG with combined proteins and save it
            l2 = list()
            for k in l: 
                if k in combined_prot:
                    for k2 in combined_prot[k]:
                        l2.append(original_name[k2])
                else:
                    l2.append(original_name[k])
            file_list_OGs.write(name_OG + '	' + ' '.join(l2) + '\n')
            # update dict per species and vector OG
            for k in l: 
                sp = prot_2_sp[k]
                if k in combined_prot:
                    nb_per_sp[sp] += len(combined_prot[k])
                    #table_og[name_OG][sp] += len(combined_prot[k])
                    for k2 in combined_prot[k]:                    
                        table_og[name_OG][sp].append(original_name[k2])
                        classified.add(original_name[k2])
                else:
                    nb_per_sp[sp] += 1
                    #table_og[name_OG][sp] += 1
                    table_og[name_OG][sp].append(original_name[k])
                    classified.add(original_name[k])
            # check if gene fusion in OG -> save name OG for each gene-fusion
            for k in l:
                if k in d_chimeric:
                    d_chimeric[k].append(name_OG)
            # count number of edges corresponding to this OG and calculate the clustering coefficient
            nb_edges = 0
            s = set(l)
            for node in l:
                for node2 in all_edges[node]:
                    if node2 in s:
                        nb_edges += 1
            clustering_coefficient = nb_edges / (len(s) * (len(s) - 1))                   
            # save OG info
            d_OG[name_OG] = [str(nb_species), str(len(l)), str(len(l2)), str(round(clustering_coefficient,4))]
       
    ## save dict old names
    utils.save_pickle(out_dir / 'OGs_in_network.pic', OGs_in_network)
          
    ## save gene-fusions
    file_fusions = open(out_dir / 'chimeric_proteins.txt', 'w+')
    file_fusions.write('#species_file	protein_name	nb_OG_fused	list_fused_OGs\n')
    for k,l in d_chimeric.items():
        file_fusions.write(d_species[str(prot_2_sp[k])] + '	' + original_name[k] + '	' + str(len(l))  + '	' + ' '.join(l) + '\n')
    
    ## save unclassified protein names
    s_all = {x for x in original_name.values()}
    d_inverse = {y:x for x,y in original_name.items()}
    unclassified = s_all.difference(classified);
    file_unclassified = open(out_dir / 'unclassified_proteins.txt', 'w+')
    for name in unclassified:
        file_unclassified.write(name + '\n')
    
    ## save statistics for each OG
    file_stats_each_OG = open(out_dir / 'statistics_per_OG.txt', 'w+')
    file_stats_each_OG.write('#OG_name	nb_species	nb_reduced_prot	nb_all_prot	clustering_coefficient\n')
    for k,l in d_OG.items():
        file_stats_each_OG.write(k + '	' + '	'.join(l) + '\n')

    ## save table OG counts
    file_stats_each_OG = open(out_dir / 'table_OGs_protein_counts.txt', 'w+')
    file_stats_each_OG.write('#OG_name	' + '	'.join(x for x in d_species.values()) + '\n')
    for og_name, d in table_og.items():
        file_stats_each_OG.write(og_name + '	' + '	'.join(str(len(x)) for x in d.values()) + '\n')

    ## save table OG names
    file_stats_each_OG = open(out_dir / 'table_OGs_protein_names.txt', 'w+')
    file_stats_each_OG.write('#OG_name	' + '	'.join(x for x in d_species.values()) + '\n')
    for og_name, d in table_og.items():
        file_stats_each_OG.write(og_name + '	' + '	'.join(' '.join(x) for x in d.values()) + '\n')
    
    ## calculate nb total prot per species, and then % assigned
    total_per_sp = collections.defaultdict(int)
    for sp in prot_2_sp.values():
        total_per_sp[sp] += 1
    perc_per_sp = {sp: (100 * nb_per_sp[sp] / total) for sp, total in total_per_sp.items()}
    
    ## save statistics for each species
    file_stats_each_species = open(out_dir / 'statistics_per_species.txt', 'w+')
    file_stats_each_species.write('#species	perc_prot_assigned nb_prot_assigned\n')
    for sp, perc in perc_per_sp.items():
        file_stats_each_species.write(d_species[sp] + '	' + str(round(perc,1)) + '	' + str(nb_per_sp[sp]) + '\n')
    
    ## save OG stats: nb sp VS nb OGs
    file_stats_OGs_sp = open(out_dir / 'statistics_nb_OGs_VS_nb_species.txt', 'w+')
    file_stats_OGs_sp.write('#nb_species	nb_OGs\n')
    for i,v in enumerate(vector_sp):
        file_stats_OGs_sp.write(str(i) + '	' + str(v) + '\n')

