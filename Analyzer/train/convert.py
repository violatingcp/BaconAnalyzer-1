import ROOT as root
import numpy as np
import root_numpy as rnp 
from keras.preprocessing import sequence

def read_branches(filenames, tree, branches, cut, treename = "events", xargs = ()):
    if not(filenames or treename) or (filenames and tree):
        PError("root_interface.read_branches", "Exactly one of filenames and tree should be specified!")
        return None
    branches_ = list(set(branches)) # remove duplicates
    if filenames:
        return rnp.root2array(filenames = filenames, 
                              treename = treename, 
                              branches = branches_, 
                              selection = cut, 
                              *xargs)
    else:
        return rnp.tree2array(tree = tree, 
                              branches = branches_, 
                              selection = cut, 
                              *xargs)
def read_tree(tree, branches, cut = None, xargs = ()):
    return read_branches(filenames = None, 
                         tree = tree,  
                         branches = branches,  
                         cut = cut, 
                         xargs = xargs)
    
#def convert(tree, branches = ['jet_pt','part_pt','part_eta','part_phi','part_pdgId','part_d0'], astype=np.float16):
def convert(tree,npad=1, branches = ['jet_pt','part_pt'], astype=np.float16):
    print "reading"
    struct_arr = read_tree(tree, branches)
    arr = []
    print "formating"
    for pEvent in range(min(len(struct_arr),100000)):
        #tmparr = [ [ struct_arr[f][pEvent][j].astype(astype) for j in range(1) ] for f in branches]
        tmparr = [ [ struct_arr[f][pEvent][j].astype(astype) for j in range(len(struct_arr[f][pEvent])) ] for f in branches]
        if npad > 1:
            tmparr = sequence.pad_sequences(tmparr, maxlen=npad)
        arr.append(tmparr)
    print "done"
    return np.array(arr)

def process_file(infile,npad, *args, **kwargs):
    f = root.TFile(infile)
    t = f.Get('Events')
    arr = convert(t,npad, *args,**kwargs)
    return arr
    #np.save(infile.replace('.root','.npy'), arr)


if __name__=='__main__':
    from sys import argv
    process_file(argv[1])
