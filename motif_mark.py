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
    parser.add_argument("-j", "--json", help="json file output", type=str, default=None)
    parser.add_argument("-p", "--plot", help="SVG plot output", type=str, default=None)
    parser.add_argument("-x", "--plotx", help="SVG x size", type=int, default=1000)
    parser.add_argument("-y", "--ploty", help="SVG y size", type=int, default=100)
    # Parse arguments
    ARGS = parser.parse_args()
    
##### Testing #####
else:   # Else test
    sys.stderr.writelines("!!!!!___RUNNING_IN_TESTING_MODE_WITH_TEST_ARGS___!!!!!\nCHANGE_TESTING_TO_FALSE\n")
    class test_args(object):
        motifs='/home/christian/lab/bgmp/motif_mark/test_motifs.txt'
        fasta='/home/christian/lab/bgmp/motif_mark/test_genes.fa'
        json='/home/christian/lab/bgmp/motif_mark/example.json'
        plot='/home/christian/lab/bgmp/motif_mark/example.svg'
        plotx=4000
        ploty=1000
    ARGS = test_args()


##### Classes #####
class motif_rec:
    'Motif record class'
    def __init__(self, name, seq, tree, keep_seq=True):
        self.name = name.strip('>')
        self.leng = len(seq)
        self.motifs = self.add_motifs(seq, tree)
        self.exons = self.add_exons(seq)
        self.seq = seq if keep_seq is True else None
    
    def add_motifs(self, seq, root):
        '''
        take the seq and the aho tree and return {motif_type_A:{match_seq_X:[pos1, pos2], match_seq_Y:[pos1]}, motif_type_B:...}
        '''
        seq=seq.lower() # sanitize
        results={} # inti the results 
        node = root # don't forget home
        for i in range(len(seq)): 
            while node != None and not seq[i] in [x for x in node.keys() if x != False and x != True]: # while we can't find the edge we are looking for, and there is a failure node
                node = node[False] # move to failure node
            if node == None: # if there is no fallback
                node = root # move to root
                continue
            node = node[seq[i]] # we found an edge to move on, lets go there.
            if True in node: # if we found a motif here
                for hit in node[True]: # for each hit
                    motif = results.setdefault(hit[1],{}) # add/move to motif type 
                    pos = motif.setdefault(hit[0],[]) # add/move to the motif seq 
                    pos.append(i-len(hit[0])+1) # Add the postion it was found
        return results
    
    def add_exons(self, seq):
        'take seq with exons in uppercase return [(exon_start,exon_stop), ...]. zero-based inclusive '
        i = 0 # int index
        exons=[] 
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

def expand_degen(seq):
    from itertools import product
    d={'A': 'A', 'B': 'CGT', 'C': 'C', 'D': 'AGT', 'G': 'G', 'H': 'ACT', 'K': 'GT', 'M': 'AC', 'N': 'GATC', 'R': 'AG', 'S': 'CG', 'T': 'T', 'V': 'ACG', 'W': 'AT', 'X': 'GATC', 'Y': 'CT'}
    return list(map("".join, product(*map(d.get, seq)))) #get the degen translation 

def grow_aho_tree(motif_file):
    '''
    Grow the Aho search tree from the motif patterns 
    '''
    #!!!!! ADD_DEGENERATES
    with open(motif_file, 'r') as fh: # open the motif file
        root = {False:None} # plant the tree root
        for line in fh: # for each line in the motif file
            path, motif = line.lower().split()[0:2] # grab the pattern and motif type 
            branch = root # move to the root
            for edge in path: # for each charater in patten 
                branch = branch.setdefault(edge, {False:None}) # move along the branch in tree and add the edges in the pattern
            branch[True] = [(path, motif)] # finally at the end of the branch,save the pattern and motif type. 
    
    queue=[]
    for k,v in root.items():
        if k != True and k != False: # grab all the nodes off root
            v[False]=root # define their failures as root
            queue.append(v) # add to work queue
    # Generate Failure links 
    while len(queue) > 0: # while there is work to be done
        wnode = queue.pop(0) # get next up
        for key, nnode in [[k,v] for k,v in wnode.items() if k != True and k != False]: # for the nodes off the working node 
            queue.append(nnode) # add them to queue      
            fnode = wnode[False] # give the failure node an easy name
            while fnode != None and not key in fnode: # while there is a failure node and the edge isn't in the failure 
                fnode = fnode[False] # move to the fnode 
            nnode[False] = fnode[key] if fnode is not None else root # add that node as nnodes failure. 
            if True in nnode[False]: # if there are matches in the failure node 
                nnode.setdefault(True, []) # make match set
                nnode[True] += nnode[False][True] if True in nnode[False] else [] # and add the things that was found in the failure node
    return root # return the tree        


def motifs2colors(motifs, sat=1,value=1,alpha=1):
    'makes a color dict useing hue for motif classes. Has ajustable saturation, value and alpha'
    from colorsys import hsv_to_rgb
    classes=set()
    for a in motifs: # for each motif
        for b in a.motifs.keys(): # get classes
            classes.add(b) # save classes
    HSV = [(x*1.0/len(classes), sat, value) for x in range(len(classes))] # generate colors by dividing hue
    return {k:hsv_to_rgb(*c)+(alpha,) for k,c in zip(classes, HSV)} # and convert to RGB and add to dict

def plot_genes(genes, output, devx=8000, devy=1000):
    import cairo as cr
    max_leng=max(x.leng for x in genes)
    num_genes=len(genes)
    # style ajustments for plotting
    line_col, line_wid = (0,0,0,1), 0.1 #base line
    exon_col, exon_wid = (0,0,0,1), 0.20 #exon boxes
    col_dict, mot_wid=  motifs2colors(motifs, value=.75, alpha=.5), 0.30 # motif box style
    seq_col,seq_font=(1,1,1,1),"Mono" # sequence text style
    lab_col, lab_font, lab_size=(0,0,0,1), "Mono", 0.25 # gene labels style
    leg_level, leg_start, leg_lab = 1-1/(num_genes+1)*.5, 0.01, 'Motifs: ' # y, x, title

    surface = cr.SVGSurface(output, devx, devy)
    ctx = cr.Context (surface)
    ctx.scale (devx, devy) # Normalizing the canvas
    for i,obj in enumerate(genes): # for each gene obj
         ctx.set_source_rgba(*line_col); ctx.set_line_width (line_wid/num_genes); ctx.move_to(0, (i+1)*1/(num_genes+1)); ctx.line_to(obj.leng/max_leng, (i+1)*1/(num_genes+1)); ctx.stroke() # draw base line across screen
         for e in obj.exons: #for each exon
             ctx.set_source_rgba(*exon_col); ctx.rectangle(e[0]/max_leng, (i+1)*1/(num_genes+1)-exon_wid/num_genes/2, (e[1]+1)/max_leng-e[0]/max_leng, exon_wid/num_genes); ctx.fill() # draw exon box
         for mk,mv in obj.motifs.items(): # for each motif type
             ctx.set_source_rgba(*col_dict[mk]) # set color
             for sk, sv in mv.items(): # for each motif seq
                 m_leng=len(sk)
                 for x in sv: # for each location
                     ctx.rectangle(x/max_leng, (i+1)*1/(num_genes+1)-mot_wid/num_genes/2, m_leng/max_leng, mot_wid/num_genes); ctx.fill() # add motif mark
         ctx.set_source_rgba(*seq_col); ctx.select_font_face(seq_font, cr.FONT_SLANT_NORMAL, cr.FONT_WEIGHT_BOLD); ctx.set_font_matrix(cr.Matrix(1/max_leng, 0, 0, line_wid/num_genes, 0, 0)) # set seq font and normilize 
         for x, letter in enumerate(obj.seq): # for each base
             xbearing, ybearing, width, height, xadvance, yadvance = (ctx.text_extents(letter)) # grab letters dims
             ctx.move_to(x/max_leng + width/2, (i+1)*1/(num_genes+1) - ybearing - height / 2); ctx.show_text(letter) # move to bases position, draw letter
         ctx.set_source_rgba(*lab_col); ctx.select_font_face(lab_font, cr.FONT_SLANT_NORMAL, cr.FONT_WEIGHT_BOLD); ctx.set_font_matrix(cr.Matrix(1/(num_genes+1)*lab_size*(devy/devx), 0, 0, 1/(num_genes+1)*lab_size, 0, 0)) # set label font and normilize
         ctx.move_to(0, (i+1)*1/(num_genes+1)-(exon_wid/num_genes)); ctx.show_text(obj.name) # draw text
    ctx.move_to(leg_start,leg_level); ctx.show_text(leg_lab);  leg_start+=ctx.text_extents(leg_lab)[4]
    for k,v in col_dict.items():
        ctx.move_to(leg_start,leg_level); ctx.set_source_rgba(*v); ctx.show_text(k)
        leg_start+=ctx.text_extents('  '+k)[4]
    surface.finish() # Output to PNG

##### MAIN #####
sys.stderr.write('Building Aho-Corasick search tree\n')
stree=grow_aho_tree(ARGS.motifs)
sys.stderr.write('Searching with tree\n')
motifs=[motif_rec(name, seq, stree) for name, seq in yield_fasta(ARGS.fasta)]
if ARGS.json is not None:
    sys.stderr.write('Writing to JSON at {}\n'.format(ARGS.json))
    from json import dump
    with open(ARGS.json, 'w') as fh:
        dump({o.name:vars(o) for o in motifs}, fh)
if ARGS.plot is not None:
    sys.stderr.write('Writing plot to {}\n'.format(ARGS.plot))
    plot_genes(motifs, ARGS.plot, ARGS.plotx, ARGS.ploty)
sys.stderr.write('DONE!\n')

