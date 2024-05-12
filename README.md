# Flexmap
Flexmap is a datastructure for k-mer based bioinformatics. Flexmap is a multimap, where multiple values can be stored for a single key. Keys are k-mers of DNA of size 15 or shorter which are stored with all. Additionally with the positions, flexmap stores additional fixed size flanking regions around the k-mer. Flexmap is the data structure used in protal, a taxonomic profiler for metagenomics with a proprietary read alignment (not published yet). The principle of the flexmap is described roughly described in the following image.

![flexmer](https://github.com/4less/flexmap/assets/15218435/4e7e63b2-5c24-400c-9be4-5837b3d14c61)

