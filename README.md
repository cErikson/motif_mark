An implementation of an Aho-Corasick search tree to find RNA/DNA motifs in relation to exons. And generate plots of the motifs in relation to the exon, or export as JSON.

**Positional Arguments**

- fasta: FASTA file of genes/mRNAs, with exons in upper case.

- motifs: A TSV of columns motif\_seq \\t motif\_identifier. Motifs may use IUPAC degenerate sequences. Motif identifier be non unique.

- mol: Molecule type: RNA or DNA

****Optional arguments****

- \-j, --json: JSON file output path/name.json

- \-p, --plot: SVG plot output path/name.svg

- \-x, --plotx: SVG x size in px, default 1000.

- \-y, --ploty, SVG y size in px, default 1000.
