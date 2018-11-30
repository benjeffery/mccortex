/*========================================================================

  File: dijkstra.c

  Program REA - Version 1.1 - July 1999

  ========================================================================

  This module contains an implementation of Dijkstra's algorithm for
  computing the shortest path from the initial node to every other node
  in a graph.

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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef DEBUG
#include <assert.h>
#endif

#include "REA.h"
#include "db_node.h"


/* The following data structures and fuctions implement a binary
   heap including the possibility of improving the priority of an
   element that is in the heap. The are used by Dijkstra function. 
*/
typedef struct {
    COST_TYPE cost;      /* Priority                               */
    uint64_t name;      /* Identifier in the range 0..dimension-1 */
    int64_t position;  /* Position in the array of elements      */
} DijkstraHeapElement;

typedef struct {
#ifdef DEBUG
    int                 dimension; /* Allocated space            */
#endif
    int64_t size;      /* Current number of elements */
    DijkstraHeapElement *info;     /* Information fields, indexed by name */
    DijkstraHeapElement **elt;     /* Pointers to the information fields  */
} DijkstraHeap;


#define NOT_IN_HEAP -1

#define PRIORITY(_i_) (heap->elt[_i_]->cost)


/*============================================================================*/
void CreateDijkstraHeap(DijkstraHeap *heap, uint64_t dimension) {

    uint64_t i;
#ifdef DEBUG
    heap->dimension = dimension;
#endif

    heap->size = 0;
    heap->elt = ctx_malloc(dimension * sizeof(heap->elt[0]));
    heap->info = ctx_malloc(dimension * sizeof(heap->info[0]));
    for (i = 0; i < dimension; i++)
        heap->info[i].position = NOT_IN_HEAP;
}


/*============================================================================*/
void FreeDijkstraHeap(DijkstraHeap *heap) {

    ctx_free(heap->info);
    ctx_free(heap->elt);
}


/*============================================================================*/
void InsertInDijkstraHeap(DijkstraHeap *heap, COST_TYPE cost, uint64_t name) {

    int64_t position, parent;
    DijkstraHeapElement *elt;

#ifdef DEBUG
    assert (heap->size < heap->dimension);
#endif

    heap->size++;
    position = heap->size;
    parent = position >> 1; /* position/2 */
    while (parent >= 1 && PRIORITY(parent) > cost) {
        heap->elt[position] = heap->elt[parent];
        heap->elt[position]->position = position;
        position = parent;
        parent = position >> 1; /* position/2 */
    }
    elt = heap->info + name;
    heap->elt[position] = elt;
    elt->position = position;
    elt->cost = cost;
    elt->name = name;

}


/*============================================================================*/
COST_TYPE DeleteBestInDijkstraHeap(DijkstraHeap *heap, uint64_t *name) {

    COST_TYPE bestCost;
    DijkstraHeapElement *item;
    int64_t position, parent;

#ifdef DEBUG
    assert (heap->size > 0);
#endif

    *name = heap->elt[1]->name;
    heap->elt[1]->position = NOT_IN_HEAP;
    bestCost = PRIORITY(1);

    item = heap->elt[heap->size];
    heap->size--;

    position = 2;
    while (position <= heap->size) {
        if (position <= heap->size - 1 && PRIORITY(position + 1) < PRIORITY(position)) position++;
        if (item->cost <= PRIORITY(position)) break;
        parent = position >> 1;   /* position/2 */
        heap->elt[parent] = heap->elt[position];
        heap->elt[parent]->position = parent;
        position = position << 1; /* 2*position */
    }
    parent = position >> 1;  /* position/2 */
    heap->elt[parent] = item;
    heap->elt[parent]->position = parent;

    return (bestCost);
}


/*============================================================================*/
int BelongsToDijkstraHeap(DijkstraHeap *heap, uint64_t name) {

    return (heap->info[name].position != NOT_IN_HEAP);
}


/*============================================================================*/
void DecreaseCostInDijkstraHeap(DijkstraHeap *heap, COST_TYPE cost, uint64_t name) {

    int64_t parent, position;

    position = heap->info[name].position;

#ifdef DEBUG
    assert (position != NOT_IN_HEAP);
    assert (cost < PRIORITY(position));
#endif

    parent = position >> 1;  /* position/2 */
    while (parent >= 1 && PRIORITY(parent) > cost) {
        heap->elt[position] = heap->elt[parent];
        heap->elt[position]->position = position;
        position = parent;
        parent = position >> 1;  /* position/2 */
    }
    heap->elt[position] = heap->info + name;
    heap->elt[position]->position = position;
    heap->elt[position]->cost = cost;

}

/*============================================================================*/
void Dijkstra(dBGraph* db_graph, uint32_t* distances, dBNode start, Orientation direction) {
    /* Computes the tree with the best path from the initial node to
       every node in the graph. Works for graphs with positive weigths.
       rea_node->bestPath was initialized by caller for every rea_node.
    */

    COST_TYPE bestCost;
    uint64_t indexNode;
    uint8_t numberEdges;
    dBNode nodes[4];
    dBNode db_node;
    Nucleotide next_bases[4];
    Edges edges;
    DijkstraHeap heap;
    uint64_t ht_size = db_graph->ht.capacity;

    CreateDijkstraHeap(&heap, ht_size * 2);

    uint64_t initialNode = rea_name(start, ht_size);
#ifdef DEBUG
    assert (initialNode != NULL);
#endif

    distances[initialNode] = 0;
    indexNode = initialNode;
    InsertInDijkstraHeap(&heap, 0, indexNode);

    while (heap.size > 0) {
        bestCost = DeleteBestInDijkstraHeap(&heap, &indexNode);
        db_node = get_db_node(indexNode, ht_size);

        edges = db_node_get_edges_union(db_graph, db_node.key);
        numberEdges = db_graph_next_nodes(
                db_graph,
                db_node_get_bkey(db_graph, db_node.key),
                db_node.orient ^ direction, //Forward or back relative to the starting kmer
                edges, (dBNode*)&nodes, next_bases);


//        char k[7];
//        BinaryKmer bk = db_node_get_bkey(db_graph, db_node.key);
//        if (db_node.orient != start.orient) bk = binary_kmer_reverse_complement(bk, 7);
//        binary_kmer_to_str(bk, 7, k);
//        status("%d %d  %s\n", bestCost, numberEdges, k);

        for (uint8_t i = 0; i < numberEdges; ++i) {
            db_node = nodes[i];
            db_node.orient = db_node.orient ^ direction;
            uint64_t dest = rea_name(db_node, ht_size);
            if (bestCost + 1 < distances[dest]) {
                indexNode = dest;
                if (BelongsToDijkstraHeap(&heap, indexNode))
                    DecreaseCostInDijkstraHeap(&heap, bestCost + 1, indexNode);
                else
                    InsertInDijkstraHeap(&heap, bestCost + 1, indexNode);
                distances[dest] = bestCost + 1;
            }
        }
    }
    FreeDijkstraHeap(&heap);
}
