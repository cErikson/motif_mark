# -*- coding: utf-8 -*-
#!/bin/python
"""
Project:Motif_Mark
@author: Christian Erikson
Motif_mark searches for a user defined set of motifs in a fasta file of genes/mRNAs. 
Using a aho-Corasick search tree the program is able to create a JSON file or SVG plot that defines the location of diffrent classes of motifs in relation to exons.
"""
##### Debug #####
testing=True

###### Imports #####
import argparse as arg
import pdb
import sys



##### Args #####
if __name__ == "__main__" and testing != True:
    parser = arg.ArgumentParser()
    # Positional mandatory arguments
    parser.add_argument("fasta", help="fasta file of genes/mRNAs, with exons in upper case", type=str)
    parser.add_argument("motifs", help="tsv of motif_seq, motif_type, optional_plot_color", type=str)
    # Optional arguments
    #parser.add_argument("-d", "--indel", help="Gap Penalty", type=int, default=None)
    #parser.add_argument("-a", "--alignment", help="Output Alignment", action='store_true')
    # Parse arguments
    ARGS = parser.parse_args()
    
##### Testing #####
else:   # Else test
    sys.stderr.writelines("!!!!!___RUNNING_IN_TESTING_MODE_WITH_TEST_ARGS___!!!!!\n")
    class test_args(object):
        motifs='/home/christian/lab/bgmp/motif_mark/test_motifs.txt'
    ARGS = test_args()


##### Classes #####
class motif_rec:
    
    def __init__(self, name, seq, tree, keep_seq=True):
        self.name = name
        self.leng = len(seq)
        self.motifs = self.add_motifs(seq, tree)
        self.exons = self.add_exons(seq)
        self.seq = seq if keep_seq is True else None

    def add_motifs(self, seq, tree):
        '''
        take the seq and the aho tree and return {motif_type_A:{match_seq_X:[pos1, pos2], match_seq_Y:[pos1]}, motif_type_B:...}
        Also look for exons, since were iterating over the seq.
        '''
        results={} # inti the results 
        for i, atom in enumerate(seq): # for each letter in the seq 
            if True in tree: # if there is a match in the tree
                for hit in tree[True]: # for each hit
                    motif = results.setdefault(hit[1],{}) # add/move to motif type 
                    pos = motif.setdefault(hit[0],[]) # add/move to the motif seq 
                    pos.append(i) # Add the postion it was found
            tree=tree.get(atom.lower(), tree[False]) # Move along an edge, if not avalible move to failure node.
        return results
                    
    def add_exons(self, seq):
        'take seq with exons in uppercase return [(exon_start,exon_stop), ...]. zero-based inclusive '
        i = 0 # int index
        exons=[] # 
        leng=len(seq)-1
        while i <= leng: # while not at end 
            if seq[i].isupper(): # if a rise 
                rise=i # save index
                while seq[i].isupper() and i < leng: # while not a fall
                    i+=1 
                exons.append((rise,i if i==leng else i-1)) # a fall happened or ended while in exon, save the exon 
            i+=1
        return(exons)
        
    def export_as_dict(self):
        pass
        
        
##### DEFS #####

def grow_aho_tree(motif_file):
    '''
    Grow the Aho tree database from the motif patterns 
    '''
    with open(motif_file, 'r') as fh: # open the motif file
        root = {False:None} # plant the tree root
        for line in fh: # for each line in  the motif file
            path, motif = line.lower().split()[0:2] # grab the pattern and motif type 
            branch = root # move to the root
            for edge in path: # for each charater in patten 
                branch = branch.setdefault(edge, {False:None}) # move along the branch in tree and add the edges in the pattern
            branch[True] = [(path, motif)] # finally at the end of the branch,save the pattern and motif type. 
    
    queue=[]
    for k,v in root.items():
        if k != True and k != False: # grab all the nodes off root
            v[False]=root # define failures as root
            queue.append(v) # add to work queue
    # Failure links 
    while len(queue) > 0: # while there is work to be done
        wnode = queue.pop(0) # get next up
        for key, nnode in [[k,v] for k,v in wnode.items() if k != True and k != False]: # for the nodes off the working node 
            queue.append(nnode) # add them to queue      
            fnode = wnode[False]
            while fnode != None and not key in fnode: # while there is a failure node and the edge isn't in failure 
                fnode = fnode[False] # move to the fnode and check conditions
            nnode[False] = fnode[key] if fnode else root # add that node as nnodes failure. 
            if True in nnode[False]:
                nnode.setdefault(True, [])
                nnode[True] += nnode[False][True] if True in nnode[False] else [] # and add the things that the failure node found
    return root # return the tree        

def yield_fasta(fasta):
    '''Make a generator that yields (seq_header, seq) for each entry''' 
    fhs=open(fasta, 'r')
    l=fhs.readline().strip()
    while l != '':
        seq = ''
        if l.startswith('>'):   #if header
            header=l  #save header
            l=fhs.readline().strip()
            while l.startswith('>') == False and l != '':
                seq+=l
                l=fhs.readline().strip() # read seq
            yield header, seq # yield data

##### MAIN #####
mot
