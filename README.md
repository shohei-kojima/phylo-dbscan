# About
This is a simple script doing clustering by DBSCAN of phylogenetic tree (.nwk).

# Requirement
Python 3.7 or later.

# Usage
```
# Options
-h, --help     show this help message and exit
-i str         Input phylogenetic tree (must be newick format).
-o str         Output file name.
-t float       Threshold for distance to make clusters.
-m int         Minimum count of members to make clusters. Default = 1.
-v, --version  show program's version number and exit
```

```
# typical usage
python main.py \
-i file.nwk \
-o clustering_result.txt \
-t 0.001 \
-m 5
```
