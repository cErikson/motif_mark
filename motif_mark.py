#!/bin/python

"""
Project:Motif_Mark
@author: Christian Erikson
Motif_mark searches for a user defined set of motifs in a fasta file of genes/mRNAs. 
Using a aho-Corasick search tree the program is able to create a JSON file or SVG plot that defines the location of diffrent classes of motifs in relation to exons.
"""

##### Debug #####
testing=False

###### Imports #####
import argparse as arg
import sys
from json import loads
from argparse import RawTextHelpFormatter


##### Args #####
if __name__ == "__main__" and testing != True:
    parser = arg.ArgumentParser(description='Motif_mark searches for a user defined set of motifs in a fasta file of genes/mRNAs. \n Using a aho-Corasick search tree the program is able to create a JSON file or SVG plot that defines the location of diffrent classes of motifs in relation to exons.',
                                epilog=""" 
    This is the default params for plotting, you can control the plotter by passing in these key:value pairs through --styles.
    Note: Passing functions will not work. JSON format.\n
     {"line_col":(0,0,0,1),             #base line: color
      "line_wid":0.15,                  # base line: width
      "exon_col":(0,0,0,1),             # exon box: color
      "exon_wid":0.25,                  # exon box: width
      "col_dict": motifs2colors(motifs, value=.75, alpha=.75), # motif box: color mappings. {"Motif_ID1": [R,G,B,A], "Motif_ID2":[1,.075,.5,.8]}
      "mot_wid":0.35,                   # motif box: width
      "seq_col":(1,1,1,1),              # sequence text: color
      "seq_font":"Mono",                # sequence text: font
      "lab_col":(0,0,0,1),              # gene labels: color
      "lab_font":"Mono",                # gene labels: font
      "lab_size": 0.25,                 # gene labels: size
      "lab_start":0.01,                 # gene labels:  x_offset
      "axis_col":(.0,0,0,1),            # axis: color,
      "axis_wid":0.1,                   # axis: width,
      "axis_tick": 10**((math.log10(max_leng)//1)-1), # axis: tick every n bp
      "axis_text_size": 0.3,            # axis: size
      "axis_text_col":(0,0,0,1),        # axis: tick col
      "axis_font":"Mono",               # axis: tick font
      "leg_box_y":0.1,                  # legend: area y
      "leg_start":0.01,                 # legend: x_offset
      "leg_lab":"Motifs: ",             # legend: leg_title
      "leg_col":(0,0,0,1),              # legend: color
      "leg_text_size":0.3,              # legend: size
      "leg_font": "Mono""               # legend: font
      }""", formatter_class=RawTextHelpFormatter)
    # Positional mandatory arguments
    parser.add_argument("fasta", help="fasta file of genes/mRNAs, with exons in upper case", type=str)
    parser.add_argument("motifs", help="tsv of motif_seq \t motif_id", type=str)
    parser.add_argument("mol", help="Molecule type: RNA or DNA", type=str)
    # Optional arguments
    parser.add_argument("-j", "--json", help="json file: output.json", type=str, default=None)
    parser.add_argument("-p", "--plot", help="SVG plot: output.svg", type=str, default=None)
    parser.add_argument("-x", "--plotx", help="SVG x size", type=int, default=1000)
    parser.add_argument("-y", "--ploty", help="SVG y size", type=int, default=100)
    parser.add_argument('--mod_color', type=loads, default=None, help='Specify custom colors. Usage: a python dict where keys are the ID and values a list of RGBA bound by 0-1, must use ". \n "{"Motif_ID1": [R,G,B,A], "Motif_ID2":[1,.075,.5,.8]}"')
    parser.add_argument('--styles', type=loads, default=None, help='Specify custom plotting parameters . Usage: "{"style": value, "lab_col":[1,0,0,1]}"')
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
        plotx=8000
        ploty=1000
        mol='RNA'
        mod_color={"RBFOX":[1,0,0,1]}
        styles={'axis_tick':50}
    ARGS = test_args()

if ARGS.json == None and ARGS.plot == None:
    sys.exit('No output setting provided by user, seek help. Exiting')

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
        seq=seq.upper() # sanitize
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

def expand_degen(motif_file, mol='DNA'):
    from itertools import product
    if mol.upper() == 'DNA':
        d={'A': 'A', 'B': 'CGT', 'C': 'C', 'D': 'AGT', 'G': 'G', 'H': 'ACT', 'K': 'GT', 'M': 'AC', 'N': 'GATC', 'R': 'AG', 'S': 'CG', 'T': 'T', 'V': 'ACG', 'W': 'AT', 'X': 'GATC', 'Y': 'CT'} 
    elif mol.upper() == 'RNA':
        d={'A': 'A', 'B': 'CGU', 'C': 'C', 'D': 'AGU', 'G': 'G', 'H': 'ACU', 'K': 'GU', 'M': 'AC', 'N': 'GAUC', 'R': 'AG', 'S': 'CG', 'U': 'U', 'V': 'ACG', 'W': 'AU', 'X': 'GAUC', 'Y': 'CU'}
    else:
        raise ValueError("Mol didn't match DNA or RNA" )
    expanded=[]
    with open(motif_file, 'r') as fh: # open the motif file
       for line in fh: # for each line in the motif file
          path, motif = line.split()[0:2] # grab the pattern and motif type 
          expanded.extend([[e,motif] for e in list(map("".join, product(*map(d.get, path.upper())))) ])
    return expanded #get the degen translation 

def grow_aho_tree(motifs):
    '''
    Grow the Aho search tree from the motif patterns 
    '''         
    #grow tree
    root = {False:None} # plant the tree root
    for path, motif in motifs : # for each line in the motif file
        branch = root # move to the root
        for edge in path: # for each charater in patten 
            branch = branch.setdefault(edge, {False:None}) # move along the branch in tree and add the edges in the pattern
        branch[True] = [(path, motif)] # finally at the end of the branch,save the pattern and motif type. 
    #prep failure link work
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
    from math import log10
    max_leng=max(x.leng for x in genes)
    num_genes=len(genes)
    # style ajustments for plotting
    style={'line_col':(0,0,0,1), 'line_wid':0.15,  #base line
          'exon_col':(0,0,0,1), 'exon_wid':0.25, #exon boxe
          'col_dict': motifs2colors(motifs, value=.75, alpha=.75), 'mot_wid':0.35, # motif box style
          'seq_col':(1,1,1,1), 'seq_font':"Mono", # sequence text style
          'lab_col':(0,0,0,1), 'lab_font':"Mono", 'lab_size': 0.25, 'lab_start':0.01, # gene labels style
          'axis_col':(.0,0,0,1), 'axis_wid':0.1, 'axis_tick': 10**((log10(max_leng)//1)-1), 'axis_text_size': 0.3, 'axis_text_col':(0,0,0,1),'axis_font':"Mono",
          'leg_box_y':0.1, 'leg_start':0.01, 'leg_lab':"Motifs: ", 'leg_col':(0,0,0,1), 'leg_text_size':0.3, 'leg_font': "Mono"# legend area y, y, x, title
          }
    if ARGS.styles is not None:
        style.update(ARGS.styles)
        
    # Custom colors 
    if ARGS.mod_color is not None:
        style['col_dict'].update(ARGS.mod_color)

    surface = cr.SVGSurface(output, devx, devy)
    ctx = cr.Context (surface)
    ctx.scale (devx, devy*(1-style['leg_box_y'])) # Normalizing the canvas
    for i,obj in enumerate(genes): # for each gene obj
         ctx.set_source_rgba(*style['line_col']); ctx.set_line_width (style['line_wid']/num_genes); ctx.move_to(0, (i+1)*1/(num_genes+1)); ctx.line_to(obj.leng/max_leng, (i+1)*1/(num_genes+1)); ctx.stroke() # draw base line across screen
         for e in obj.exons: #for each exon
             ctx.set_source_rgba(*style['exon_col']); ctx.rectangle(e[0]/max_leng, (i+1)*1/(num_genes+1)-style['exon_wid']/num_genes/2, (e[1]+1)/max_leng-e[0]/max_leng, style['exon_wid']/num_genes); ctx.fill() # draw exon box
         for mk,mv in obj.motifs.items(): # for each motif type
             ctx.set_source_rgba(*style['col_dict'][mk]) # set color
             for sk, sv in mv.items(): # for each motif seq
                 m_leng=len(sk)
                 for x in sv: # for each location
                     ctx.rectangle(x/max_leng, (i+1)*1/(num_genes+1)-style['mot_wid']/num_genes/2, m_leng/max_leng, style['mot_wid']/num_genes); ctx.fill() # add motif mark
         ctx.set_source_rgba(*style['seq_col']); ctx.select_font_face(style['seq_font'], cr.FONT_SLANT_NORMAL, cr.FONT_WEIGHT_BOLD); ctx.set_font_matrix(cr.Matrix(1/max_leng, 0, 0, style['line_wid']/num_genes, 0, 0)) # set seq font and normilize 
         for x, letter in enumerate(obj.seq): # for each base
             xbearing, ybearing, width, height, xadvance, yadvance = (ctx.text_extents(letter)) # grab letters dims
             ctx.move_to(x/max_leng + width/2, (i+1)*1/(num_genes+1) - ybearing - height / 2); ctx.show_text(letter) # move to bases position, draw letter
         ctx.set_source_rgba(*style['lab_col']); ctx.select_font_face(style['lab_font'], cr.FONT_SLANT_NORMAL, cr.FONT_WEIGHT_BOLD); ctx.set_font_matrix(cr.Matrix(1/(num_genes+1)*style['lab_size']*(devy/devx), 0, 0, 1/(num_genes+1)*style['lab_size'], 0, 0)) # set label font and normilize
         ctx.move_to(0+style['lab_start'], (i+1)*1/(num_genes+1)-(style['exon_wid']/num_genes)); ctx.show_text(obj.name) # draw text
    ctx.scale (1, 1+style['leg_box_y']) # Normalizing to full canvas
    ctx.set_source_rgba(*style['axis_col']); ctx.set_line_width (style['axis_wid']/num_genes); ctx.move_to(0, 1-style['leg_box_y']); ctx.line_to(1, 1-style['leg_box_y']); ctx.stroke() # draw base line across screen
    ctx.set_source_rgba(*style['axis_text_col']); ctx.select_font_face(style['axis_font'], cr.FONT_SLANT_NORMAL, cr.FONT_WEIGHT_BOLD); ctx.set_font_matrix(cr.Matrix(style['leg_box_y']*style['axis_text_size']*(devy/devx), 0, 0, style['leg_box_y']*style['axis_text_size'], 0, 0)) # set axis font and normilize 
    for x in range(int(max_leng//style['axis_tick'])+1): # for each tick
             xbearing, ybearing, width, height, xadvance, yadvance = (ctx.text_extents(letter)) # grab letters dims
             ctx.move_to((style['axis_tick']*x)/max_leng, 1-style['leg_box_y']/3*2); ctx.show_text(str(int(style['axis_tick']*x))) # move to tick position, draw letter
    ctx.select_font_face(style['leg_font'], cr.FONT_SLANT_NORMAL, cr.FONT_WEIGHT_BOLD); ctx.set_font_matrix(cr.Matrix(style['leg_box_y']*style['leg_text_size']*(devy/devx), 0, 0, style['leg_box_y']*style['leg_text_size'], 0, 0)) # set axis font and normilize 
    ctx.set_source_rgba(*style['leg_col']);ctx.move_to(style['leg_start'],1-style['leg_box_y']/3+ctx.text_extents(style['leg_lab'])[4]/2); ctx.show_text(style['leg_lab']);  style['leg_start']+=ctx.text_extents(style['leg_lab'])[4] #start legend
    for k,v in style['col_dict'].items():
        ctx.move_to(style['leg_start'],1-style['leg_box_y']/3+ctx.text_extents(style['leg_lab'])[4]/2); ctx.set_source_rgba(*v); ctx.show_text(k) # wtite legend
        style['leg_start']+=ctx.text_extents('  '+k)[4] # find next position
    surface.finish() # Output to PNG

##### MAIN #####
sys.stderr.write('Building Aho-Corasick search tree\n')
stree=grow_aho_tree(expand_degen(ARGS.motifs, ARGS.mol))
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

