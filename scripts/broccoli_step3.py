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
from statistics import mean
import itertools
import gzip
import pickle
from multiprocessing import Pool as ThreadPool 
from scripts import utils
try:
    from ete3 import PhyloTree
except:
    sys.exit("\n            ERROR: the ete3 library is not installed\n\n")




def step3_orthology_network(minh, shar, overl, nbsp, nt):

    # convert the parameters to global variables (horrible hack)
    global min_nb_hits, limit_shared, limit_overlap, limit_nb_sp, nb_threads
    min_nb_hits, limit_shared, limit_overlap, limit_nb_sp, nb_threads = minh, shar, overl, nbsp, nt
    
    print('\n --- STEP 3: phylomes\n')
    print(' ## parameters')
    print(' min nb hits    : ' + str(min_nb_hits))
    print(' fusion shared  : ' + str(limit_shared))
    print(' fusion nb sp   : ' + str(limit_nb_sp))
    print(' fusion overlap : ' + str(limit_overlap))
    print(' threads        : ' + str(nb_threads))
    print('\n ## analysis')

    ## create output directory (or empty it if it already exists)
    utils.create_out_dir('./dir_step3')

    ## get all species
    all_species = utils.get_pickle('./dir_step1/species_index.pic')
        
    ## check directory
    files_trees = pre_checking('./dir_step2/')
            
    ## create log file
    global log_file
    log_file = open('./dir_step3/log_step3.txt', 'w+')
    
    ## load all search outputs
    print(' load similarity search outputs')
    global search_outputs
    search_outputs = utils.get_multi_pickle('./dir_step2/', '_output.pic')
    
    global prot_2_sp
    ## load prot_name 2 species dict (string version)
    prot_2_sp = utils.get_pickle('./dir_step2/prot_str_2_species.pic')
       
    ## extract ortho-para from dir_step2 (takes some time)
    ortho_pairs, para_pairs = extract_ortho_para(files_trees)

    ## load prot_name 2 species dict (integer version)
    prot_2_sp = utils.get_pickle('./dir_step2/prot_int_2_species.pic')
      
    ## build network
    global all_nodes, all_edges 
    all_nodes, all_edges = build_network(ortho_pairs, para_pairs)
    
    # free memory
    ortho_pairs, para_pairs = False, False  
        
    ## define maximum number of edges to consider for lcc (= 2 * nb_species or 10 if less species) -> improve speed
    global limit_degree, limit_nb_max
    limit_degree, limit_nb_max = get_limit_lcc(all_species)

    ## calculate the local clustering coefficient of each node
    global all_lcc
    all_lcc = multithread_lcc(all_nodes, nb_threads)
            
    ## get all connected_components
    list_cc = utils.get_connected_components(all_edges, all_nodes)
    
    ## analyse each connected_components
    network_nodes, nb_node_removed, gene_fusions, communities = multithread_analyse_cc(list_cc, nb_threads)
    
    ## save OG lists, fusions and stats
    save_outputs(communities, gene_fusions, all_edges, all_species) 
        
    print('')   
    
   
def pre_checking(directory):
       
    # check if directory exists
    if not os.path.isdir(directory):
        sys.exit("\n            ERROR STEP 3: the directory ./dir_step2 does not exist.\n\n")
   
    # list the input _blast_ortho.pic and _trees.pic pickle files (only names of _trees.pic files are returned)
    list_1, list_2 = list(), list()
    for file in os.listdir(directory):
        if file.endswith('_blast_ortho.pic'):
            list_1.append(file)
        elif file.endswith('_trees.pic'):
            list_2.append(file)
    list_2.sort()
    
    # print error
    if len(list_1) != len(list_2):
        sys.exit("\n            ERROR STEP 3: the number of *_trees.pic and *_blast_ortho.pic are different\n\n")
    return list_2

# --------------------- #

def extract_ortho_para(l_trees):
        
    print(' extract ortho from NO tree files')
    # load pickle files
    tmp_d = utils.get_multi_pickle('./dir_step2/', '_blast_ortho.pic')

    # extract ortho
    d_ortho = collections.defaultdict(int)
    for l in tmp_d.values():
        for sub in itertools.combinations(l, 2):
            d_ortho[ sub[0]+'-'+sub[1] ] += 1
    
    # free memory
    tmp_d = None
               
    # extract ortho from tree files
    print(' extract ortho from tree files')    
    pool = ThreadPool(nb_threads) 
    tmp_res = pool.map_async(extract_ortho_from_trees, l_trees, chunksize=1)
    results_2 = tmp_res.get()
    pool.close() 
    pool.join()    
    
    # unpack ortho and save them
    print(' load ortho extracted from tree files')
    ortho_para = list()
    for t in results_2:
        xx = pickle.loads(t[0])
        for pair, nb in xx.items():
            d_ortho[pair] += nb 
        # save ortho@para
        ortho_para += t[1]
        
    # free memory
    results_2 = None
    
    # remove ortho found only once and create para dict()
    print(' remove ortho found only once')
    d_ortho = {k:v for k, v in d_ortho.items() if v != 1}
    
    # extract para from tree files
    print(' extract para from tree files')    
    d_para = extract_para_from_trees(ortho_para, d_ortho)
    
    return d_ortho, d_para


def extract_ortho_from_trees(filename):
    
    # prepare output variables
    out_ortho = collections.defaultdict(int)
    l_ortho_para = list()
    
    # load dict of trees
    tmp_d = utils.get_pickle('./dir_step2/' + filename)
    
    # analyse trees 1 by 1
    for ref_leaf, newick in tmp_d.items():
        
        # load tree and get all leaves
        tree = PhyloTree(newick)
        all_leaves = {leaf.name for leaf in tree}
                                   
        # get all leaves from last interesting nodes      
        ref_node = tree.search_nodes(name = ref_leaf)[0]
        ortho = custom_species_overlap(ref_node)
        
        # add ref_leaf to ortho in case no good node selected
        if len(ortho) == 0:
            ortho.add(ref_leaf)
        
        # save ortho
        xx = list(ortho)    
        xx.sort()    
        for sub in itertools.combinations(xx, 2):
            out_ortho[sub[0] + '-' + sub[1]] += 1
        
        # get para
        para = all_leaves - ortho
        
        # save ortho@para if there is a paralogous group
        if para:
           l_ortho_para.append(' '.join(ortho) + '@' + ' '.join(para))
        
    return [pickle.dumps(out_ortho), l_ortho_para]


def custom_species_overlap(node):
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
        if len(sister_leaves)==0:
            current = current.up
            continue
        # Gets species at both sides of event
        #sister_spcs        = set([n.species for n in sister_leaves])
        sister_spcs        = set(prot_2_sp[leaf.name] for leaf in sister_leaves)
        overlaped_spces    = len(browsed_spcs & sister_spcs)
        all_spcs           = len(browsed_spcs | sister_spcs)
        sp_only_in_sister  = len(sister_spcs - browsed_spcs)
        sp_only_in_browsed = len(browsed_spcs - sister_spcs)

        # Updates browsed species
        browsed_spcs   |= sister_spcs
        browsed_leaves |= sister_leaves
        sister_leaves  = set([])
        
        if all_spcs == 1 or overlaped_spces == 0 or (overlaped_spces == 1 and sp_only_in_sister >= 2 and sp_only_in_browsed >= 2):
            last_good = set(n.name for n in browsed_leaves)
        
        # And keep ascending
        current = current.up
    return last_good


def extract_para_from_trees(l_ortho_para, ortho_d):
    # create para dict()
    out_para = {x:0 for x in ortho_d}
    
    # analyse trees 1 by 1
    for st in l_ortho_para:
        ortho, para = st.split('@')
        l_ortho = ortho.split(' ')
        l_para = para.split(' ')
        
        for name1 in l_para:
            for name2 in l_ortho:
                if name1 < name2:
                    combined_name = name1 + '-' + name2
                else:
                    combined_name = name2 + '-' + name1
                # save it if in ortho
                try:
                    out_para[combined_name] += 1
                except:
                    pass           
    return out_para

# --------------------- #

def build_network(d_ortho_pairs, d_para_pairs):
     
    d_edges = collections.defaultdict(dict)
    d_nodes = dict()
    nb_edges = 0
    
    print(' build network')
    for k, nb_ortho in d_ortho_pairs.items():
        # build edge and nodes if nb ortho superior to nb para
        if nb_ortho > d_para_pairs[k]:
            # extract names and convert them to integers
            insert = k.split('-') 
            name1, name2 = int(insert[0]), int(insert[1])
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
    log_file.write('network size:\n' + str(len(d_nodes)) + ' nodes\n' + str(nb_edges) + ' edges\n\n')
        
    # get maximum connected value for each node
    max_ortho = dict()
    for node in d_nodes:
        max_ortho[node] = max(d_edges[node].items(), key=operator.itemgetter(1))[1]
        
    # convert values to ratio nb_ortho / max_ortho
    for node1 in d_edges:
        to_remove = set()
        for node2, nb_ortho in d_edges[node1].items():
            d_edges[node1][node2] = nb_ortho / max_ortho[node1]
    
    # sort all edges of each node by values (used to break ties in the label propagation)
    for node1 in d_edges:
        d_sorted = dict()
        for key, value in sorted(d_edges[node1].items(), key=lambda x: x[1], reverse=True):
            d_sorted[key] = value
        # save it
        d_edges[node1] = d_sorted
    
    return d_nodes, d_edges

# --------------------- #

def multithread_lcc(d_nodes, n_threads):
    # split node dict in a list of lists
    list_nodes = [x for x in d_nodes]
    new_list_of_lists = [list_nodes[i::n_threads] for i in range(n_threads)]  
    
    # start multithreading
    print(' compute lcc for each node')
    pool = ThreadPool(n_threads) 
    tmp_res = pool.map_async(calculate_lcc, new_list_of_lists, chunksize=1)
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


def calculate_lcc(l):
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

# --------------------- #

def multithread_analyse_cc(l_cc, n_threads):
    
    print(' analyse each connected components')
    # start multithreading
    pool = ThreadPool(n_threads) 
    tmp_res = pool.map_async(analyse_cc, l_cc, chunksize=1)
    results_2 = tmp_res.get()
    pool.close() 
    pool.join()                
    
    # get all results together
    final_list    = list()
    GF            = dict()
    network_nodes = set()
    nb_removed    = 0
    for l in results_2:
        # get the communities
        for l2 in l[0]:
            final_list.append(l2)
        # get number of nodes removed
        nb_removed += len(l[1])
        # get the gene-fusions
        for k in l[2]:
            GF[k] = list()
    
    # print results
    print('      _ ' + str(len(l_cc)) + ' connected components')
    print('      _ ' + str(len(final_list)) + ' communities')
    print('      _ ' + str(nb_removed) + ' nodes removed')
    print('      _ ' + str(len(GF)) + ' gene fusions')    
    
    # save to log file
    log_file.write('network analysis:\n' + str(len(l_cc)) + ' connected components\n' + str(len(final_list)) + ' communities\n' + str(nb_removed) + ' nodes removed\n' + str(len(GF)) + ' gene fusions')
     
    return network_nodes, nb_removed, GF, final_list
    
    
def analyse_cc(ll):
    
    final_list      = list()
    nb_corrected    = 0
    nodes_removed   = set()
    gene_fusions    = dict()
    
    # check number of and species nodes in cc -> apply LPA method if nb_species > 1 and nb_nodes > 4
    nb_species = fast_count_species(ll)
    nb_nodes   = len(ll)
    
    if nb_nodes < 4 or nb_species == 1:
        # save list
        final_list.append(ll)
                
    else:
                        
        # sort nodes in list to get stable results
        ll.sort()
        
        # create dictionary of nodes in cc (name as key and lcc as value)
        nodes_in_cc = {k:all_lcc[k] for k in ll}
            
        # apply LPA to identify communities
        node2label = label_propagation(nodes_in_cc, all_edges)
                             
        # reconstruct communities (i.e. labels present in values of node2label)
        label2node = create_label_dict(node2label)
                   
        # build list of communities (list of list of nodes)
        tmp_com = [l for l in label2node.values()]
            
        if len(tmp_com) == 1:
            # save community  
            for l in tmp_com:
                final_list.append(l)
                
        else:
            # remove spurious hits
            tmp2_com, nodes_removed = remove_false_positives(tmp_com, all_edges)
        
            # detect gene-fusions and modify communities accordingly
            tmp3_com, gene_fusions = detect_chimeric_proteins(tmp2_com, all_edges)
            
            # save community  
            for l in tmp3_com:
                final_list.append(l)
                    
    return [final_list, nodes_removed, gene_fusions]

                
def label_propagation(d_nodes, d_edges):

    # sort the node dict by lcc (decreasing order) + replace lcc value by label
    tmp_nodes_in_cc = dict()
    for key, value in sorted(d_nodes.items(), key=lambda x: x[1], reverse=True):
        tmp_nodes_in_cc[key] = key
    
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
                
            for v in d_edges[node]:
                label = tmp_nodes_in_cc[v]
                all_neighbors_labels[label] += d_edges[node][v]
                        
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


def remove_false_positives(ll_com, d_edges):    
    nodes_removed = set()
    # each communities 1 by 1
    for i,l in enumerate(ll_com):
        set_node = set(l)
        modif = False
        # each node 1 by 1
        for node in l:
            # count how many hits for this node belong to this community
            nb_found = 0
            for t in search_outputs[node]:
                if t[0] in set_node:
                    nb_found += 1
                    if nb_found == min_nb_hits:
                        break
            # remove node if it didn't reach the limit
            if nb_found < min_nb_hits:
                set_node.remove(node)
                nodes_removed.add(node)
                modif = True
        # update the community if modification
        if modif:
            ll_com[i] = list(set_node)
    return ll_com, nodes_removed
     

def detect_chimeric_proteins(ll_com, d_edges):
    # prepare dict limit nodes per OG (nb nodes * limit ratio) AND node 2 OG id
    limit_nodes_per_OG = dict()
    node_2_OG          = dict()
    for i,l in enumerate(ll_com):
        limit_nodes_per_OG[i] = len(l) * limit_shared
        for k in l:
             node_2_OG[k] = i
    
    # isolated nodes shared by several OGs
    shared_nodes = dict()
    nodes_in_cc = {x for l in ll_com for x in l}
    for node1 in nodes_in_cc:
        counts = collections.defaultdict(int)
        for node2 in d_edges[node1]:
            if node2 in nodes_in_cc:          # the node might not exist anymore if spurious hit
                counts[node_2_OG[node2]] += 1
        connected_OGs = {x for x,v in counts.items() if v > limit_nodes_per_OG[x]}
        if len(connected_OGs) > 1:
            shared_nodes[node1] = connected_OGs
    
    # check each shared node 1 by 1
    fusions = set()
    modif = False
    for node, conn_OGs in shared_nodes.items():
        # initialise variable
        d_found = collections.defaultdict(set)
        d_start = dict()
        d_end   = dict()
        # analyse search output of this protein
        for l in search_outputs[node]:
            node2, start, end = l[0], l[1], l[2]
            # do not take node2 if is in shared_nodes (possible multiple gene-fusions that would screw the analyses)
            if node2 in nodes_in_cc and node2 not in shared_nodes:
                OG = node_2_OG[node2]
                # only consider the node if it belongs to one of the connected OG
                if OG in conn_OGs:
                    d_found[OG].add(node2)
                    if OG not in d_start:
                        d_start[OG] = start
                        d_end[OG]   = end
                    else:
                        if start < d_start[OG]:
                            d_start[OG] = start
                        if end > d_end[OG]:
                            d_end[OG] = end
        # select connected OGs if more or equal to limit_sp
        l_good_og = list()
        for og, s in d_found.items():
            nb_sp = count_species(s)
            if nb_sp >= limit_nb_sp:
                l_good_og.append(og)
        # compare overlap between OGs
        for i, og1 in enumerate(l_good_og):
            for n in range(i+1,len(l_good_og)):
                og2 = l_good_og[n]
                # calculate the overlap of hits
                if d_start[og1] >= d_start[og2]:
                    overlap = d_end[og2] - d_start[og1]
                else:
                    overlap = d_end[og1] - d_start[og2]
                # if overlap, save fusion and add it to the 2 OGs
                if overlap < limit_overlap:
                    fusions.add(node)
                    ll_com[og1].append(node)
                    ll_com[og2].append(node)
                    modif = True
        
    # remove redundant nodes in communities if gene-fusions have been added
    if modif:
        for i,l in enumerate(ll_com):
            ll_com[i] = list(set(l))
    
    return ll_com, fusions


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


def save_outputs(l_com, d_fusions, d_edges, d_species):    
    
    ## load original and combined names
    original_name = utils.get_pickle('./dir_step1/original_names.pic')
    combined_prot = utils.get_pickle('./dir_step1/combined_names.pic')
            
    # create vector nb_species as index and nb_OG as value
    vector_sp = [0] * (len(d_species) + 1)

    ## save the lists of OGs and get fusion info and get OG info
    d_OG = collections.defaultdict(list)
    c = 0
    OGs_in_network = dict()
    file_list_OGs = open('./dir_step3/orthologous_groups.txt', 'w+')
    file_list_OGs.write('#OG_name	protein_names\n')
    for l in l_com:
        nb_species = count_species(l)
        # keep OG if more than 1 species
        if nb_species > 1:
            vector_sp[nb_species] += 1
            # create name OG
            c += 1
            name_OG = 'OG_' + str(c)
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
            # check if gene fusion in OG -> save name OG for each gene-fusion
            for k in l:
                if k in d_fusions:
                    d_fusions[k].append(name_OG)
            # count number of edges corresponding to this OG and calculate the clustering coefficient
            nb_edges = 0
            s = set(l)
            for node in l:
                for node2 in d_edges[node]:
                    if node2 in s:
                        nb_edges += 1
            clustering_coefficient = nb_edges / (len(s) * (len(s) - 1))                   
            # save OG info
            d_OG[name_OG] = [str(nb_species), str(len(l)), str(len(l2)), str(round(clustering_coefficient,4))]
    
    
    ## save dict old names
    utils.save_pickle('./dir_step3/OGs_in_network.txt', OGs_in_network)
          
    ## save gene-fusions
    file_fusions = open('./dir_step3/chimeric_proteins.txt', 'w+')
    file_fusions.write('#species_file	original_name	nb_OG_fused	list_fused_OGs\n')
    for k,l in d_fusions.items():
        file_fusions.write(d_species[str(prot_2_sp[k])] + '	' + original_name[k] + '	' + str(len(l))  + '	' + ' '.join(l) + '\n')
        
    ## save statistics for each OG
    file_stats_each_OG = open('./dir_step3/statistics_each_OG.txt', 'w+')
    file_stats_each_OG.write('#OG_name	nb_species	nb_reduced_prot	nb_all_prot	clustering_coefficient\n')
    for k,l in d_OG.items():
        file_stats_each_OG.write(k + '	' + '	'.join(l) + '\n')
    
    ## save OG stats: nb sp VS nb OGs
    file_stats_OGs_sp = open('./dir_step3/statistics_nb_OGs_VS_nb_species.txt', 'w+')
    file_stats_OGs_sp.write('#nb_species	nb_OGs\n')
    for i,v in enumerate(vector_sp):
        file_stats_OGs_sp.write(str(i) + '	' + str(v) + '\n')

