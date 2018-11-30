/*========================================================================

  File: REA.c

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

#include <db_node.h>
#include "global.h"

#ifdef DEBUG
#include <assert.h>
#endif

#include "REA.h"
#include "dijkstra.h"
#include "db_graph.h"
#include "ctx_alloc.h"



/*============================================================================
  GLOBAL VARIABLES
  ============================================================================*/

/* A buffer for heap elements is dinamically allocated by the main function  */
/* with size the number of arcs in the graph (which is the worst case bound) */
/* New heaps are initialized using this buffer by the function CreateHeap    */
Path **heapsBuffer;
Path **firstFreeHeap;

/* For experimental purposes we avoid to use the malloc system for new paths */
/* New paths are taken from this buffer by the function CreatePath           */
Path *pathsBuffer;
uint64_t MAX_PATHS = 0;
uint64_t numberUsedPaths = 0;

uint64_t rea_name(dBNode node, uint64_t ht_size) {
    return node.key + node.orient * ht_size;
}

dBNode get_db_node(uint64_t rea_node, uint64_t ht_size) {
    dBNode n = {
            .key = rea_node % ht_size,
            .orient = rea_node / ht_size
    };
    return n;
}

/*============================================================================*/
Path *CreatePath(REANode *node, Path *backPath, COST_TYPE cost) {

    Path *newPath;

#ifdef DEBUG
    assert (node != NULL);
#endif

    newPath = &pathsBuffer[numberUsedPaths];

    if (++numberUsedPaths > MAX_PATHS) {
        die("Exceeded MAX_PATHS");
    }

    newPath->backPath = backPath;
    newPath->cost = cost;
    newPath->lastREANode = node;

    return newPath;

}


/*============================================================================*/
//void PrintPath(Path *path) {
//    /* Prints the path in reverse order (with the final node first). */
//    Path *backPath;
//
//#ifdef DEBUG
//    assert (path != NULL);
//#endif
//    if (path->cost < INFINITY_COST) {
//        for (backPath = path; backPath != NULL; backPath = backPath->backPath)
//            printf("%i-", backPath->lastREANode->name);
//        printf(" \t(Cost: %i)", path->cost);
//    }
//
//}


#define COST(_i_) (H->elt[_i_]->cost)

/*============================================================================*/
void CreateHeap(Heap *H, int dimension) {
    /* Allocates memory for the heap elements from the heaps buffer
       and initializes the number of element inside the heap to 0 */
#ifdef DEBUG
    assert (dimension >= 0);
    H->dimension = dimension;
#endif
    H->size = 0;
    H->elt = firstFreeHeap;
    firstFreeHeap += dimension;
}


/*============================================================================*/
void PreInsertInHeap(Heap *H, Path *elt) {
    /* Inserts a new element in the heap in time O(1)
       without preserving the heap property */
#ifdef DEBUG
    assert (H->size < H->dimension);
#endif

    H->elt[H->size++] = elt;
}


/*============================================================================*/
void BuildHeap(Heap *H) {
    /* Obtains the heap property in time linear with the heap size */

    register int64_t i;
    for (i = H->size / 2; i >= 0; i--) {
        register Path *item = H->elt[i];
        register int64_t j = 2 * i + 1;
        while (j < H->size) {
            if (j < H->size - 1 && COST(j + 1) < COST(j)) j++;
            if (item->cost <= COST(j)) break;
            H->elt[(j - 1) / 2] = H->elt[j];
            j = 2 * j + 1;
        }
        H->elt[(j - 1) / 2] = item;
    }
}


/*============================================================================*/
Path *DeleteBestInHeap(Heap *H) {
    /* Returns the element with minimum cost in the heap and deletes
       it from the heap, restoring the heap property in worst case time
       logarithmic with the heap size */
    Path *best;

    if (H->size == 0) return (NULL);
    best = H->elt[0];
    H->elt[0] = H->elt[--H->size];
    {
        register Path *item = H->elt[0];
        register int j = 1;
        while (j < H->size) {
            if (j < H->size - 1 && COST(j + 1) < COST(j)) j++;
            if (item->cost <= COST(j)) break;
            H->elt[(j - 1) / 2] = H->elt[j];
            j = 2 * j + 1;
        }
        H->elt[(j - 1) / 2] = item;
    }
    return (best);
}

/*============================================================================*/
Path *BestInHeap(Heap *H) {
    /* Returns the element with minimum cost in the heap without deleting it */

    if (H->size == 0) return (NULL);
    return H->elt[0];
}

/*============================================================================*/
void UpdateBestInHeap(Heap *H, Path *elt) {
    /* Deletes the element with minimum cost in the heap and inserts a new one
   preserving the heap property, in worst case time logarithmic with the
   heap size */

    register Path *item = elt;
    register int j = 1;
    while (j < H->size) {
        if (j < H->size - 1 && COST(j + 1) < COST(j)) j++;
        if (item->cost <= COST(j)) break;
        H->elt[(j - 1) / 2] = H->elt[j];
        j = 2 * j + 1;
    }
    H->elt[(j - 1) / 2] = item;
}


/*============================================================================*/
void InitializeCandidateSet(REANode *rea_node, REANode *rea_nodes, dBGraph *db_graph, uint64_t ht_size) {
    /* The set of candidates is  initialized with the best path from
       each  predecessor node, except  the  one from which  the best
       path at the current node  comes. If the current node is the
       initial node in  the graph, all its predecessor nodes (if some
       exists) provide a candidate */

    Path *path,
            *newCand;
    dBNode nodes[4];
    Nucleotide next_bases[4];
    Edges edges;
    uint8_t number_edges;

    dBNode db_node = get_db_node(rea_node, ht_size);

    edges = db_node_get_edges_union(db_graph, db_node.key);
    number_edges = db_graph_next_nodes(
            db_graph,
            db_node_get_bkey(db_graph, db_node.key),
            db_node.orient, //REVERSE relative to the starting kmer
            edges, (dBNode *) &nodes, next_bases);


    CreateHeap(&rea_node->heap, number_edges);

    for (uint8_t i = 0; i < number_edges; ++i) {
        REANode source = rea_nodes[rea_name(nodes[i], ht_size)];
        if ((path = source.bestPath) != rea_node->bestPath->backPath) {
            /* It is important to compare pointers and not costs or node names,
               because  several candidates could came from the same node (if
               the graph has parallel arcs) and/or with the same cost. */
            newCand = CreatePath(rea_node, path, path->cost + 1);
            PreInsertInHeap(&rea_node->heap, newCand);
        }
    }
    BuildHeap(&rea_node->heap);
}


/*============================================================================*/
Path *NextPath(Path *path, REANode *rea_nodes, dBGraph *db_graph, uint64_t ht_size) {
    /* Central routine of the Recursive Enumeration Algorithm: computes the
       next path from the initial node to the same node in which the argument
       path ends, assuming that it has not been computed before.
    */

    REANode *node = NULL;
    Path *backPath = NULL,
            *nextBackPath = NULL,
            *bestCand = NULL,
            *newCand = NULL;
    COST_TYPE arcCost;

#ifdef DEBUG
    assert (path != NULL);
    assert (path->nextPath == NULL);
    assert (path->lastREANode != NULL);
#endif

    node = path->lastREANode;

#ifdef DEBUG
    printf ("\nComputing next path at node %i\n", node->name);
#endif

    if (node->heap.elt == NULL)
        /* This is done here instead of in Dijkstra function, so that it is */
        /* only done for nodes in which alternative paths are required.     */
        InitializeCandidateSet(node, rea_nodes, db_graph, ht_size);

    if ((backPath = path->backPath) != NULL) {
        nextBackPath = backPath->nextPath;
        if (nextBackPath == NULL)
            nextBackPath = NextPath(backPath, rea_nodes, db_graph, ht_size);
        if (nextBackPath != NULL) {
            arcCost = path->cost - backPath->cost;
            newCand = CreatePath(node, nextBackPath, nextBackPath->cost + arcCost);
        }
    }

    if (newCand == NULL)
        bestCand = DeleteBestInHeap(&node->heap);
    else {
        bestCand = BestInHeap(&node->heap);
        if (bestCand != NULL && bestCand->cost < newCand->cost)
            UpdateBestInHeap(&node->heap, newCand);
        else
            bestCand = newCand;
    }

    if (bestCand == NULL)
        return NULL;

    /* Adds the best candidate to the list of paths ending in the same node */
    bestCand->nextPath = NULL;
    path->nextPath = bestCand;

#ifdef DEBUG
    printf ("\nNext path at node %i:\t", node->name);
    PrintPath (bestCand);
#endif

    return bestCand;

}

void REA(dBGraph *db_graph, dBNode start, dBNode target, uint8_t *kmer_mask, uint64_t max_length) {
    Path *path;
    REANode *rea_nodes;
    uint64_t ht_size = db_graph->ht.capacity;

    /********** Allocates memory for heaps ************************************/
    /* Worst case is number of edges, but that is very unlikely to happen */
//    heapsBuffer = (Path **) ctx_malloc(ht_size * sizeof(*heapsBuffer));
//    firstFreeHeap = heapsBuffer;
    uint32_t *forward_distances = ctx_calloc(ht_size * 2, sizeof(uint32_t));
    uint32_t *reverse_distances = ctx_calloc(ht_size * 2, sizeof(uint32_t));
//    pathsBuffer = ctx_malloc(sizeof(Path) * ht_size * 4);
//    MAX_PATHS = ht_size * 4;
    /* Initializes for each node the counters of input and output arcs, */
    /* the heap of candidate paths, the best path and the node name     */
    for (uint64_t i = 0; i < ht_size * 2; ++i) {
        forward_distances[i] = INFINITY_COST;
        reverse_distances[i] = INFINITY_COST;
    }

    /****************** Computes the shortest path tree ***********************/
    status("Running forward Dijkstra BFS....");
    Dijkstra(db_graph, forward_distances, start, FORWARD);
    status("Shortest: %d", forward_distances[rea_name(target, ht_size)]);
    status("Running reverse Dijkstra BFS....");
    Dijkstra(db_graph, reverse_distances, target, REVERSE);
    status("Shortest: %d", reverse_distances[rea_name(start, ht_size)]);

    for (uint64_t i = 0; i < ht_size; i++) {
        if (forward_distances[i] + reverse_distances[i] < max_length ||
            forward_distances[i + ht_size] + reverse_distances[i + ht_size] < max_length) {
            bitset_set(kmer_mask, i);
        }
    }

    ctx_free(forward_distances);
    ctx_free(reverse_distances);
#ifdef DEBUG
        for (node = graph.node, numberREANodes = graph.numberREANodes; numberREANodes != 0;
             node++, numberREANodes--) {
          printf ("\nBest Path for FinalREANode=%i: ", node->name);
          PrintPath (node->bestPath);
        }
#endif
        /******************** Computes the K shortest paths ***********************/
//    uint64_t i = 2;
//    uint64_t last  = 0;
//    path = finalNode.bestPath;
//    while (path->cost < max_length && path != NULL) {
//        path = NextPath(path, rea_nodes, db_graph, ht_size);
//        if (last < path->cost) {
//            last = path->cost;
//            status("%lu", path->cost);
//        }
//        i++;
//    }
//    status("%lu paths", i);

        /************ Prints the computed paths and time counters ******************/
//    i = 1;
//    path = graph.finalREANode->bestPath;
//    while (i <= numberPaths && path != NULL && path->cost < INFINITY_COST) {
//        printf("\nN=%i:\t", i);
//        if (showPaths == 1) PrintPath(path);
//        printf(" \t(CumulatedSeconds: %.2f)", (float) cumulatedSeconds[i - 1]);
//        path = path->nextPath;
//        i++;
//    }
//    printf("\nTotalExecutionTime: %.2f\n", (float) cumulatedSeconds[i - 2]);
//    return (0);

    }
