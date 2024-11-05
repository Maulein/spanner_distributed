# Sampling Based Spanner Algorithm
**Generating compact k-spanners in distributed setup**

Copyright (c) 2024 Maulein
This code is part of our work titled “Scalable algorithms for compact spanners on real-world graphs"

This work is built using the Parallel Boost Grpah Library (PBGL)

//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

**Organization**

The repository contains code files for our work titled “Scalable algorithms for compact low-stretch spanners on real-world graphs". The implementations are built using the PBGL. 

Our work includes implementations of distriubted graph algorithm for producing graph spanners under the variants: Baseline, Hybrid Spanner and Sampling Based Algorithms using the adjacency list format and CSR format. The code files (spanner.cpp, bfs1hop.hpp, visitors.hpp) for these variants are located under the appropriate directories. 

**Compilation and Running Code** 

Bring up the PBGL library (refer [https://www.boost.org/doc/libs/1_82_0/libs/graph_parallel/doc/html/index.html])

For a specific implemetation, copy the files from the corresponding folder to PATH/boost_1_82_0/

Compile spanner.cpp and run (refer PBGL setup for details)
