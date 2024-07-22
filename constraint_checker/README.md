## Constraint Checking

```
Failed setting feasible start rates/dates after 10 attempts. Aborting.
```
The error above is produced when data submitted to a treePL analysis possess some issue(s) that prevents the program from generating some initial rates/dates that can be optimized through penalized likehood.

From experience dealing with troubleshooting many such issues, one (or more) of the following is invariably the culprit:

1. Inconsistent constraints. For example, calibrating a descendant node with a larger minimum age than that specified for one of its ancestral nodes.

2. Having multiple (perhaps conflicting) constraints applied to the same node. This almost certainly involves a problem involving the MRCA statement. In trees of many hundreds or thousands of tips, two different pairs of tips can easily have the same MRCA. In a slightly different manner, the relationships in the tree may not coincide exactly with expectations; for example, if a (say) family is not recovered to be monophyletic in the present tree, then selection of a different representative taxon within that expected group may define a different MRCA node. [Aside: the code [here](https://gist.github.com/josephwb/f3d35f8833a07f71002af7726b12652b) may help diagnose the last situation].

3. The tree is unrooted. The meaning of 'ancestral' and 'descendant' nodes is ambiguous when a tree is unrooted. The code below can tell you if your tree is unrooted, but (because taxonomy is not involved) does not know the tree _shoiuld_ be rooted.


To help deal with these issues, I've put together some R code that you can use to vet your data. Example data (with built-in issues) is also provided.

### How To Use

Put your data in the present directory; for now we will use the example data provided. In R, set your working directory to the present directory (I'll use my own system path here):

```
setwd("/home/josephwb/Work/Phylogenetics/treePL/constraint_checker");
```
(Yes, I use semicolons in R, and I will never _not_ do so ;)). Next, source the functions file:
```
source("check_fossils.R");
```
This will make two functions available:
```
check_constraints_consistent(phy, cnstrnts);
```
which does the actual constraint checking, and
```
writeTreePLConfig(phyname, nsites, ndconstrnts, nthreads=8, writeIt=TRUE);
```
which can write out a treePL config file (if desired).

Let's run it with the example data, shall we? Let us read in our data. First, the tree:
```
phy <- read.tree("phy.tre");
```
The tree includes node labels (matches 'Constraint.name' in constraint table, below), which can be visualized in FigTree so the exmaple can be easily followed.

Now load the csv file containing information on the names, positions, and ages of constraints:
```
fossils <- read.csv("constraint_table.csv", stringsAsFactors=FALSE, na.strings=c("", " ", "NA"));
```
We can look at the structure of the constraint table by just printing it out:
```
fossils;
   Constraint.name   mrca_L   mrca_R  Min  Max                                      info
1             Root  taxon_1 taxon_39 70.1 86.5                                      <NA>
2     CladeA_crown  taxon_1 taxon_14 26.1 31.3      younger than descendant node FossilW
3          FossilW  taxon_5 taxon_14 34.7 52.3    older than ancestral node CladeA_crown
4     CladeB_crown taxon_15 taxon_25 43.7 57.6                                      <NA>
5      CladeC_stem taxon_15 taxon_39 49.4 56.2                                      <NA>
6     CladeC_crown taxon_26 taxon_39 38.9 43.7                                      <NA>
7          FossilX taxon_26 taxon_38 27.9 34.6      younger than descendant node FossilY
8          FossilY taxon_26 taxon_35 35.1 37.2         older than ancestral node FossilX
9          FossilZ taxon_26 taxon_34 27.3 29.2 dates same node as FossilQ but is younger
10         FossilQ taxon_27 taxon_33 31.9 34.6   dates same node as FossilZ but is older
```
The table _must_ contain the columns (in any order):
'Constraint.name', 'mrca_L', 'mrca_R', 'Min', and 'Max'. 'mrca_L' and 'mrca_R' are two taxa which share a common ancestor at some node in the tree; the 'L' (left) and 'R' (right) are arbitrary, and can be thought of instead as 'tip1' and 'tip2' in any order. The `info` column present in this example is not necessary, and here just explains how introduced issues should be interperted. 

Ok, let's try it out:
```
res <- check_constraints_consistent(phy=phy, cnstrnts=fossils);
[1] "Error. The following constraints map to the same node:"
[1] "   'FossilZ' (minage=27.3)."
[1] "   'FossilQ' (minage=31.9)."
[1] "Problem with constraint 'FossilW' (minage = 34.7)."
[1] "  Ancestor 'CladeA_crown' has a younger minage (26.1)."
[1] "Problem with constraint 'FossilY' (minage = 35.1)."
[1] "  Ancestor 'FossilX' has a younger minage (27.9)."
[1] "Problem with constraint 'FossilQ' (minage = 31.9)."
[1] "  Ancestor 'FossilX' has a younger minage (27.9)."
```
As can be seen, several issues are found. Let's look at the results:
```
res;
     constraint node  min  max     rtax     ltax  keep             notes
1          Root   40 70.1 86.5 taxon_39  taxon_1  TRUE              good
2  CladeA_crown   41 26.1 31.3 taxon_14  taxon_1 FALSE       invalidated
3       FossilW   45 34.7 52.3 taxon_14  taxon_5  TRUE              good
4  CladeB_crown   55 43.7 57.6 taxon_25 taxon_15  TRUE              good
5   CladeC_stem   54 49.4 56.2 taxon_39 taxon_15  TRUE              good
6  CladeC_crown   65 38.9 43.7 taxon_39 taxon_26  TRUE              good
7       FossilX   66 27.9 34.6 taxon_38 taxon_26 FALSE       invalidated
8       FossilY   67 35.1 37.2 taxon_35 taxon_26  TRUE              good
9       FossilZ   68 27.3 29.2 taxon_34 taxon_26 FALSE redundant,younger
10      FossilQ   68 31.9 34.6 taxon_33 taxon_27  TRUE              good
```
If you like, you can write the results to file for your records:
```
write.csv(res, row.names=FALSE, quote=FALSE, file="fossil_analysis_results.csv");
```
If you are happy with the results, we can clean up the table a bit:
```
# prune the invalidated constraints
res <- res[res$keep,];
# and pitch columns that serve no future purpose
res <- within(res, rm(keep, notes));
```
If desired, you can even write a treePL config file using the updated constraints. *NOTE:* make sure to use variable values that match your data! Below I use made up values for the tree name and number of sites. Note also that some generic default settings are, er, set, so make sure you update these to what is approproriate (e.g., smooth, thorough, log_pen, cvstart, etc.; you can see what those (and other) settings are [here](https://github.com/blackrim/treePL/wiki/Run-Options)).
```
writeTreePLConfig(phyname="My_awesome_phylogeny.tre", nsites=131313, ndconstrnts=res, nthreads=8);
```
For the example data, the following config file is produced:
```
## Tree Information ##

treefile = My_awesome_phylogeny.tre
outfile = My_awesome_phylogeny.d8s.tre
cvoutfile = My_awesome_phylogeny.tre.cv.out


## Temporal Calibrations ##

mrca = Root taxon_39 taxon_1
max = Root 86.5
min = Root 70.1

mrca = FossilW taxon_14 taxon_5
max = FossilW 52.3
min = FossilW 34.7

mrca = CladeB_crown taxon_25 taxon_15
max = CladeB_crown 57.6
min = CladeB_crown 43.7

mrca = CladeC_stem taxon_39 taxon_15
max = CladeC_stem 56.2
min = CladeC_stem 49.4

mrca = CladeC_crown taxon_39 taxon_26
max = CladeC_crown 43.7
min = CladeC_crown 38.9

mrca = FossilY taxon_35 taxon_26
max = FossilY 37.2
min = FossilY 35.1

mrca = FossilQ taxon_33 taxon_27
max = FossilQ 34.6
min = FossilQ 31.9


## Run Settings ##

# Make sure the settings below match your data!

numsites = 131313
smooth = 10
thorough
log_pen
nthreads = 8
plsimaniter = 200000
#cv
#randomcv
cvstart = 1000
cvstop = 0.01
prime
```