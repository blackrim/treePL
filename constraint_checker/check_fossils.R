require(ape);
require(phangorn);

# purpose: given a table of fossil constraints and a tree:
# 1. check that fossil ages are consistent with topology (i.e., descendant not older than ancestors)
# 2. check if fossils map to same node (get rid of younger one)
# 3. produce output files for dating (e.g., treePL)

## functions are at the top of the file, commands (and example data) are below

################################################################################

# map each constraint in the tree. then check if any ancestors are younger than descendants
# this does not currently check max ages
check_constraints_consistent <- function (phy, cnstrnts) {

    rtax <- cnstrnts$mrca_R;
    ltax <- cnstrnts$mrca_L;
    
    # map mrcas to nodes in the tree
    nds <- NULL; # mrca nodes
    for (i in 1:length(cnstrnts[,1])) {
        nds <- c(nds, getMRCA(phy, c(ltax[i], rtax[i])));
    }
    
    # set up df to return
    df <- as.data.frame(cbind(constraint=cnstrnts$Constraint.name, node=nds, min=cnstrnts$Min,
                              max=cnstrnts$Max, rtax=rtax, ltax=ltax), stringsAsFactors=FALSE);
    class(df[,2]) <- "numeric";
    class(df[,3]) <- "numeric";
    class(df[,4]) <- "numeric";
    
    # decisions and notes. default to accept
    keep <- rep(TRUE, length(df[,1]));
    notes <- rep("good", length(df[,1]));
    
    ## error-checking
    # first, check for constraints mapping to the same node
    if (any(duplicated(nds))) {
        d.nds <- as.numeric(names(which(table(nds) > 1)));
        for (i in 1:length(d.nds)) {
            idx <- which(df$node == d.nds[i]);
            i.ages <- df$min[idx];
            print("Error. The following constraints map to the same node:");
            for (j in 1:length(idx)) {
                print(paste0("   '", cnstrnts$Constraint.name[idx[j]], "' (minage=", cnstrnts$Min[idx[j]],
                             ")."));
            }
            toDrop <- which(i.ages != max(i.ages)); # *** have not considered identical mins yet
            keep[idx[toDrop]] <- "FALSE";
            notes[idx[toDrop]] <- "redundant,younger";
        }
    }
    
    # finally, check for consistency (i.e. ancestors are not younger than descendants)
    for (i in 1:length(nds)) {
        pnds <- Ancestors(phy, df$node[i], "all");
        # throw out nodes without constraints
        pnds <- pnds[pnds %in% df$node];
        if (length(pnds) > 0) {
            i.age <- df$min[i];
            p.ages <- df$min[match(pnds, df$node)];
            if (any(i.age >= p.ages)) {
                print(paste0("Problem with constraint '", df$constraint[i], "' (minage = ", i.age, ")."));
                # give details
                idx <- which(p.ages <= i.age);
                for (ii in 1:length(idx)) {
                    iidx <- which(df$node == pnds[idx[ii]]);
                    print(paste0("  Ancestor '", df$constraint[iidx],
                                 "' has a younger minage (", df$min[iidx], ")."));
                    keep[iidx] <- FALSE;
                    notes[iidx] <- "invalidated";
                }
            }
        }
    }
    
    df <- cbind(df, keep, notes);
    class(df$keep) <- "logical";
    
    if (all(df$keep)) {
        print("No issues detected!");
    }
    return(df);
}

# write a treePL config file to pwd for an individual source tree
# ndconstrnts is a dataframe with cols: mrca1, mrca2, age, cname, node
# ndconstrnts is a dataframe with cols: constraint, node, min, max, rtax, ltax
writeTreePLConfig <- function (phyname, nsites, ndconstrnts, nthreads=8, writeIt=TRUE) {
    
    treeInfo <- paste(
        "## Tree Information ##\n",
        paste0("treefile = ", phyname),
        paste0("outfile = ", sub(".tre", ".d8s.tre", phyname)),
        paste0("cvoutfile = ", phyname, ".cv.out\n"), sep="\n");
    
    cals <- "\n\n## Temporal Calibrations ##\n";
    # calibrations
    for (j in 1:length(ndconstrnts[,1])) {
        # calib name
        calName <- as.character(ndconstrnts$constraint[j]);
        indiv <- paste0("\nmrca = ", calName, " ", ndconstrnts$rtax[j], " ", ndconstrnts$ltax[j], "\n");
        
        # max constraint
        indiv <- paste0(indiv, paste0("max = ", calName, " ", ndconstrnts$max[j], "\n"));
        
        # min constraint
        indiv <- paste0(indiv, paste0("min = ", calName, " ", ndconstrnts$min[j], "\n"));
        cals <- paste0(cals, indiv);
    }
    # run conditions
    runTreePL <- paste(
        "\n\n## Run Settings ##\n",
        "# Make sure the settings below match your data!\n",
        paste0("numsites = ", nsites),
        "smooth = 10", # what is a reasonable general value?
        "thorough",
        "log_pen",
        paste0("nthreads = ", nthreads),
        "plsimaniter = 200000",
        "#cv",
        "#randomcv",
        "cvstart = 1000",
        "cvstop = 0.01",
        "prime", sep="\n");
    
    toWrite <- paste0(treeInfo, cals, runTreePL);
    print(toWrite);
    if (writeIt) { # set to false for development
        write(toWrite, file=paste0(sub(".tre", ".treepl.config", phyname)));
    }
}

################################################################################

## dummy example data
# tree. includes node labels (matches constraint table, below), which can be visualized in FigTree
phy <- read.tree("phy.tre");
# constraint data. must minimally have column names (in any order):
# 'Constraint.name', 'mrca_L', 'mrca_R', 'Min', 'Max'
# can have extra columns. this example includes 'info' which 
fossils <- read.csv("constraint_table.csv", stringsAsFactors=FALSE, na.strings=c("", " ", "NA"));

# let's run it!
res <- check_constraints_consistent(phy=phy, cnstrnts=fossils);
# look at the results
res;
# write a copy of the results for records
write.csv(res, row.names=FALSE, quote=FALSE, file="fossil_analysis_results.csv");

# if you are happy with the results, we can prune the invalidated constraints
res <- res[res$keep,];
# and pitch columns that serve no future purpose
res <- within(res, rm(keep, notes));

# write a treePL config file. make sure to use variable values that match your data!
writeTreePLConfig(phyname="My_awesome_phylogeny.tre", nsites=131313, ndconstrnts=res, nthreads=8);

