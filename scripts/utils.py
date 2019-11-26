import os
import pickle
import gc


''' to create dir (and remove previous version if exists) ''' 
def create_out_dir(dir):
    if os.path.isdir(dir):
        os.system('rm -R ' + dir)
    os.mkdir(dir)


''' to parse FASTA files '''
def read_fasta(fasta_content):
    name, seq = None, []
    for line in fasta_content:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name.replace('>',''), ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name.replace('>',''), ''.join(seq))


''' save dict on disk using pickle '''
def save_pickle(file_name, d):
    with open(file_name, 'wb') as handle:
        pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)    


''' load pickle dict from disk '''
def get_pickle(file_name):
    d = dict()
    with open(file_name, 'rb') as content:
        gc.disable()  # disable garbage collector
        d.update(pickle.load(content))
        gc.enable()   # enable garbage collector again
    return d


''' load multiple pickle dict from disk (input is directory and a string matching the end of the file names) '''
def get_multi_pickle(dir_, str_):
    d = dict()
    for file in os.listdir(dir_):
        if file.endswith(str_):
            with open(dir_ + file, 'rb') as content:
                gc.disable()  # disable garbage collector
                d.update(pickle.load(content))
                gc.enable()   # enable garbage collector again
    return d


''' to get connected components from network '''
def get_connected_components(edges, nodes):
    l = list()
    seen = set()
    for v in nodes:
        if v not in seen:
            c = set(BFS_algo(edges, v))
            l.append(list(c))
            seen.update(c)
    return l

def BFS_algo(edges, source):
    seen = set()
    nextlevel = {source}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in seen:
                yield v
                seen.add(v)
                nextlevel.update(edges[v])
