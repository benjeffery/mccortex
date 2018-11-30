/*========================================================================

  File: REA.h

  Program REA - Version 1.1 - July 1999

  ========================================================================

  This module contains an implementation of the Recursive Enumeration
  Algorithm (REA) that enumerates (by increasing weight) the N shortest
  paths in weighted graphs. The algorithm is described in: 

    "Computing the K Shortest Paths: a New Algorithm and an
    Experimental Comparison" by Victor Jimenez and Andres Marzal,
    3rd Workshop on Algorithm Engineering, London, July 1999.
    To be published by Springer-Verlag in the LNCS series.

  The sets of candidate paths are implemented using binary heaps.

  ========================================================================

    Copyright (C) 1999 Victor Jimenez and Andres Marzal

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (file COPYING); if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    You can contact the authors at:

           Victor Jimenez and Andres Marzal
           Dept. de Informatica, Universitat Jaume I, Castellon, Spain
           {vjimenez,amarzal}@inf.uji.es

    Updated information about this software will be available at the
    following www address:

           http://terra.act.uji.es/REA
 
  ========================================================================= */


#include <db_graph.h>

#ifndef _REA_H_INCLUDED

#define COST_TYPE uint32_t     /* You can change this type but you should verify
                             its use in printf's, scanf's and castings */

#define INFINITY_COST 1000000000 /* Instead of INT_MAX, to avoid overflow  */
                                  /* when adding some cost to INFINITY_COST */


/* The following data structures represent a graph in memory providing
   access to incoming and outgoing arcs for each REANode, and also allow
   to represent multiple shortest paths from the initial REANode to each REANode. 
*/

typedef struct Path {
  COST_TYPE   cost;      /* Path cost                                         */
  struct REANode *lastREANode; /* REANode in which this path ends                      */
  struct Path *backPath; /* Prefix path, ending in a predecessor REANode         */
  struct Path *nextPath; /* Next path in the list of computed paths from the  */
                         /* initial REANode to the REANode in which this paths ends */
} Path;   

typedef struct Heap {
#ifdef DEBUG
  int dimension; /* Allocated space                                         */
#endif
  int64_t size;      /* Current number of elements                              */
  Path **elt;    /* Elements are pointers to paths, allocated by CreateHeap */
} Heap;

typedef struct REANode {
  uint64_t name;          /* An integer in the range 1..numNodes */
  Path   *bestPath;     /* First path in the list of computed paths     */
                        /* from the initial REANode to this REANode           */
  Heap   heap;          /* Set of candidate paths (binary heap)         */
} REANode;

Path *CreatePath (REANode *REANode, Path *backPath, COST_TYPE cost);
uint64_t rea_name(dBNode node, uint64_t ht_size);
dBNode get_db_node(uint64_t rea_node, uint64_t ht_size);
void REA(dBGraph* db_graph, dBNode start, dBNode target, uint8_t *kmer_mask, uint64_t num_paths);

#define _REA_H_INCLUDED
#endif
