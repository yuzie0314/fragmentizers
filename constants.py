#!/usr/bin/env python3
# -*- coding: utf-8 -*-
maxMolSize = 50
preprocessMols = True
useHs = False # If True, include hydrogen atoms in fragment definition. Default is False.
keepFragIdsPickle = True
numThreads = 1 # Number of CPU threads to use for catalog generation. Default is 1.


# verbose: If True, print progress updates during catalog generation. Default is False.
# pklFile: If specified, save the generated catalog as a binary pickle file with the specified filename. Default is None.
# minHAC: Minimum heavy atom count for fragments to be included in the catalog. Default is 3.
# maxHAC: Maximum heavy atom count for fragments to be included in the catalog. Default is 8.
# maxNumFragments: Maximum number of fragments to include in the catalog. Default is 1000.
# pathLength: Maximum number of bonds to include in fragment definition. Default is 3.
# tolerance: Similarity threshold for fragment clustering. Default is 0.3.
