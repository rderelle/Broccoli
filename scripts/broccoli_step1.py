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
import itertools
from multiprocessing import Pool as ThreadPool 
from scripts import utils



def step1_kmer_clustering(dir, ext, ms, lk, ma, nt):
   
    # convert the parameters to global variables (horrible hack)
    global directory, extension, min_seq, length_kmer, min_aa, nb_threads
    directory, extension, min_seq, length_kmer, min_aa, nb_threads = dir, ext, ms, lk, ma, nt
    
    print('\n --- STEP 1: kmer clustering\n')
    print(' # parameters')
    print(' proteomes dir : ' + directory)
    print(' kmer size     : ' + str(length_kmer))
    print(' min size seq  : ' + str(min_seq))
    print(' min nb aa     : ' + str(min_aa))

    ## create output directory (delete it first if already exists)
    utils.create_out_dir('./dir_step1')
    
    ## check directory and files
    print('\n # check input files')
    global dict_files, list_files, list_start
    dict_files, list_files, list_start = pre_checking(directory, extension)
    
	## analyse each fasta file (multithreading)
    print ('\n # kmer clustering\n ' + str(len(list_files)) + ' proteomes on ' + str(nb_threads) + ' threads')
    pool = ThreadPool(nb_threads) 
    tmp_res = pool.map_async(process_file, list_files, chunksize=1)
    results_2 = tmp_res.get()
    pool.close() 
    pool.join()
    
    ## create log files
    log_file = open('./dir_step1/log_step1.txt', 'w+')
    log_file.write('#index	file_name	nb_initial	nb_short	nb_final\n')
    
    ## save log file and combine other info
    combined = dict()
    names    = dict()
    nb_final = 0
    for l in results_2:
        log_file.write('	'.join(l[:5]) + '\n')
        names.update(l[5])
        combined.update(l[6])
        nb_final += int(l[4])
    
    ## save pickle files
    utils.save_pickle('./dir_step1/combined_names.pic', combined)
    utils.save_pickle('./dir_step1/original_names.pic', names)
    utils.save_pickle('./dir_step1/species_index.pic', dict_files)
    
    print(' -> ' + str(nb_final) + ' proteins saved for the next step')
    print ('')
 
     
def pre_checking(directory, ext):
    
    # check if input directory exists
    if not os.path.isdir(directory):
        sys.exit('\n            ERROR STEP 1: the directory \' ' + directory + ' \' does not exist.\n\n')
    
    # list the input files inside that directory and (i) count and (ii) check the name
    else: 
        l_files  = list()
        d_nb_seq = dict()
        s_names  = set()
        for file in os.listdir(directory):
            if file.endswith(ext):
                # save file to list
                l_files.append( (file, os.path.getsize(directory + file)) )
                nb_seq = 0
                # check names
                content = open(directory + file, "r")
                for line in content:
                    if line.startswith('>'):
                        nb_seq += 1
                        name = line.split(' ')[0]
                        if name in s_names:
                            print('problem : the protein name \' ' + name.replace('\n','') + ' \' from file ' + file + ' already exists')
                        else:
                             s_names.add(name)
                # save nb seq
                d_nb_seq[file] = nb_seq
    
    # sort the list by (reverse) size
    l_files.sort(key=lambda x: x[1], reverse=True)
    # re-populate list with just filenames
    l_start = list()
    counter = 0
    for i,t in enumerate(l_files):
        l_files[i] = t[0] 
        l_start.append(counter)  
        counter += d_nb_seq[t[0]]
                
    # test if files
    if len(l_files) == 0:
        sys.exit('\n            ERROR STEP 1: there is no input fasta file (*' + ext + ') in this directory.\n\n')
    
    # convert to dictionary
    d_files = {str(i):k for i,k in enumerate(l_files)}

    # print log
    print(' ' + str(len(l_files)) + ' input files')
    print(' ' + str(counter) + ' sequences')
    
    return d_files, l_files, l_start

# --------------------- #

def process_file(filename):
    
    ## get back info
    index = list_files.index(filename)
    counter = list_start[index]
    
    ## extract fasta sequences with new IDs
    all_names, all_seq = create_dict_seq(filename, directory, counter)
    initial = len(all_names)
    
    ## remove short seq
    to_remove = list()
    for name, seq in all_seq.items(): 
        if len(seq) < min_seq:
            to_remove.append(name)
    for name in to_remove:
        del all_names[name], all_seq[name]
    
    ## extract kmers
    kmer_2_id = dict()
    all_kmers = collections.defaultdict(list)
    count = 0
    for name, seq in all_seq.items(): 
        # extract k-mers 
        for i in range(len(seq) - length_kmer + 1):
        #for i in range(round((len(seq) - length_kmer + 1)/2)):
            #j = i * 2
            #kmer = seq[j:j + length_kmer]
            kmer = seq[i:i + length_kmer]
            # check number of distinct aa
            distinct_aa = set(kmer)
            # only consider kmer with at least min_aa distinct aa, and no 'X' amd no '*'
            if len(distinct_aa) >= min_aa and 'X' not in distinct_aa and '*' not in distinct_aa:
                if kmer in kmer_2_id:
                    id = kmer_2_id[kmer]
                else:
                    id = count
                    kmer_2_id[kmer] = id
                    count += 1
                # add the protein name to the kmer
                all_kmers[id].append(name)
    
    # simplify 1: remove kmers present in only 1 seq
    kd1 = list()
    for k in all_kmers:
        if len(all_kmers[k]) > 1:
           kd1.append(all_kmers[k])
    
    # simplify 2: remove kmers found in the same list of seq
    data = sorted(kd1)
    kd2 = [k for k, v in itertools.groupby(data)]
    
    # build edges
    all_edges = collections.defaultdict(dict)
    all_nodes = set()
    for l in kd2:
        for i,name_1 in enumerate(l):
            all_nodes.add(name_1)
            for n in range(i+1, len(l)):
                all_edges[name_1][l[n]] = ''
                all_edges[l[n]][name_1] = ''
    
    # get connected components
    list_cc = utils.get_connected_components(all_edges, all_nodes)
    
    # add isolated seq (i.e. seq without kmer shared with another seq)
    for name in all_names:
        if name not in all_nodes:
            list_cc.append([name])
    
    # save output
    all_combined = dict()
    out_file = open('./dir_step1/' + str(index) + '.fas', 'w+')
    for l in list_cc:
        if len(l) == 1:
            name = l[0]
            out_file.write('>' + str(name) + '\n' + all_seq[name] + '\n')
        else:
            # select longest seq
            max = 0
            ref = ''
            for k in l:
                seq_length = len(all_seq[k])
                if seq_length > max:
                    max = seq_length
                    ref = k
            # save it
            out_file.write('>' + str(ref) + '\n' + all_seq[ref] + '\n')
            # save combined names
            all_combined[ref] = tuple(l)
    out_file.close()
    
    return [str(index), filename, str(initial), str(len(to_remove)), str(len(list_cc)), all_names, all_combined]

# --------------------- #

def create_dict_seq(file, dir, c):
    d_name  = dict()
    d_seq   = dict()
    with open(dir + file) as fasta_content:
        for name_seq, seq in utils.read_fasta(fasta_content):
            d_name[c] = name_seq.split(' ')[0]
            d_seq[c]  = seq
            c += 1
    return d_name, d_seq

