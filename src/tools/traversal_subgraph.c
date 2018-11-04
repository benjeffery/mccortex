#include "global.h"
#include "subgraph.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "prune_nodes.h"
#include "seq_loading_stats.h"
#include "util.h"

typedef struct {
    dBNode nodes[4];
    unsigned char current;
    unsigned char num;
} StackEntry;

typedef struct {
    const dBGraph *const db_graph; // graph we are operating on
    uint8_t *const kmer_mask; // bitset of visited kmers
    SeqLoadingStats stats;
    BinaryKmer start_kmer, end_kmer;
    dBNode start_node, end_node;
} TraversalSubgraphBuilder;


static void traversal_subgraph_builder_alloc(TraversalSubgraphBuilder *builder,
                                             uint8_t *kmer_mask,
                                             const dBGraph *graph) {
    TraversalSubgraphBuilder tmp = {.db_graph = graph,
            .kmer_mask = kmer_mask,
    };

    memcpy(builder, &tmp, sizeof(TraversalSubgraphBuilder));
    seq_loading_stats_init(&builder->stats);
}

static void traversal_subgraph_builder_dealloc(TraversalSubgraphBuilder *builder) {
//    db_node_buf_dealloc(&builder->nbufs[0]);
//    db_node_buf_dealloc(&builder->nbufs[1]);
//    db_node_buf_dealloc(&builder->unitig_buf);
}


static void get_start_from_read(read_t *r1, read_t *r2,
                                uint8_t qoffset1, uint8_t qoffset2, void *ptr) {
    TraversalSubgraphBuilder *builder = (TraversalSubgraphBuilder *) ptr;
    const dBGraph *db_graph = builder->db_graph;
    if (r1->seq.end < db_graph->kmer_size) {
        die("Start sequence is shorter than kmer size");
    }

    //Cut off the string at kmer size
    r1->seq.b[db_graph->kmer_size] = "\0";
    builder->start_kmer = binary_kmer_from_str(r1->seq.b, db_graph->kmer_size);
    builder->start_node = db_graph_find(db_graph, builder->start_kmer);
    if (builder->start_node.key == HASH_NOT_FOUND) {
        die("Start kmer not found in graph");
    }
}

static void get_end_from_read(read_t *r1, read_t *r2,
                              uint8_t qoffset1, uint8_t qoffset2, void *ptr) {
    TraversalSubgraphBuilder *builder = (TraversalSubgraphBuilder *) ptr;
    const dBGraph *db_graph = builder->db_graph;
    if (r1->seq.end < db_graph->kmer_size) {
        die("End sequence is shorter than kmer size");
    }
    builder->end_kmer = binary_kmer_from_str(r1->seq.b + r1->seq.end - db_graph->kmer_size,
                                             db_graph->kmer_size);
    builder->end_node = db_graph_find(db_graph, builder->end_kmer);
    if (builder->end_node.key == HASH_NOT_FOUND) {
        die("End kmer not found in graph");
    }

}

void traverse_and_mark(dBGraph *db_graph, size_t nthreads, size_t max_depth,
                       bool invert,
                       size_t stack_mem, uint8_t *kmer_mask,
                       seq_file_t **startfiles, size_t num_start_files,
                       seq_file_t **endfiles, size_t num_end_files) {
    TraversalSubgraphBuilder builder;
    traversal_subgraph_builder_alloc(&builder, kmer_mask, db_graph);

    read_t r1;
    if (seq_read_alloc(&r1) == NULL)
        die("Out of memory");

    seq_parse_se_sf(startfiles[0], 0, &r1, get_start_from_read, &builder);
    seq_parse_se_sf(endfiles[0], 0, &r1, get_end_from_read, &builder);
    seq_read_dealloc(&r1);

    StackEntry *stack = ctx_malloc(sizeof(StackEntry) * max_depth);
    uint8_t *seen_forward = ctx_calloc(roundup_bits2bytes(db_graph->ht.capacity), 1);
    uint8_t *seen_reverse = ctx_calloc(roundup_bits2bytes(db_graph->ht.capacity), 1);
    uint8_t *keep_forward = ctx_calloc(roundup_bits2bytes(db_graph->ht.capacity), 1);
    uint8_t *keep_reverse = ctx_calloc(roundup_bits2bytes(db_graph->ht.capacity), 1);

    size_t pos = 0;
    Edges edges = db_node_get_edges_union(db_graph, builder.start_node.key);
    Nucleotide next_bases[4];
    stack[pos].current = 0;
    stack[pos].num = 1;
    stack[pos].nodes[0] = builder.start_node;
    bitset_set(builder.end_node.orient == FORWARD ? keep_forward : keep_reverse, builder.end_node.key);
    do {
        bool go_back = false;
        StackEntry* entry = &stack[pos];
        //Have we exhausted the branches of the parent node?
        if (entry->current < entry->num) {
            dBNode node = entry->nodes[entry->current];

//            char k[7];
//            BinaryKmer bk = db_node_get_bkey(db_graph, node.key);
//            if (node.orient != builder.start_node.orient)
//                bk = binary_kmer_reverse_complement(bk, 7);
//            binary_kmer_to_str(bk, 7, k);
//            printf("%s", k);

            if(!(bitset_get(node.orient == FORWARD ? seen_forward : seen_reverse, node.key)) && !(
                    node.key == builder.end_node.key && node.orient == builder.end_node.orient)) {
                bitset_set(node.orient == FORWARD ? seen_forward : seen_reverse, node.key);

//                printf(" A");
                edges = db_node_get_edges_union(db_graph, node.key);
                //Create the next step in the traverse
                ++pos;
                stack[pos].current = 0;
                stack[pos].num = db_graph_next_nodes(
                        db_graph,
                        db_node_get_bkey(db_graph, node.key),
                        node.orient, //Forward relative to the starting kmer
                        edges, &stack[pos].nodes, next_bases);
                //Next time we visit this node take the next branch
                ++(entry->current);
            } else {
                //Node seen before, or end kmer go back up
//                printf(" SEEN");
                go_back = true;
            }
            //If we got the end kmer, or one that is on a path to the end kmer (i.e. marked to keep, then keep this path)
            if (bitset_get(node.orient == FORWARD ? keep_forward : keep_reverse, node.key)) {
//                printf(" KEEP");
                for (size_t i = pos-1; i > 0; i--) {
                    node = stack[i].nodes[stack[i].current-1]; //-1 here as current is pointing to the unexplored branch
                    if (bitset_get(node.orient == FORWARD ? keep_forward : keep_reverse, node.key)) {
                        break;
                    }
                    bitset_set(node.orient == FORWARD ? keep_forward : keep_reverse, node.key);
                }
                go_back = true;
            }
//            printf("\n");
        } else {
            //We have exhausted this node, go back up
            go_back = true;
        }
        if (go_back) {
            //Pop the stack and mark the node as unvisited
            --pos;
            dBNode node = stack[pos].nodes[stack[pos].current-1];
            bitset_del(node.orient == FORWARD ? seen_forward : seen_reverse, node.key);
        }
        //Stop when we get back to the root node
    } while (pos > 0 && pos < (max_depth - 1));

    bitset_set(builder.start_node.orient == FORWARD ? keep_forward : keep_reverse, builder.start_node.key);

    size_t num_bytes = roundup_bits2bytes(db_graph->ht.capacity);
    for(size_t i = 0; i < num_bytes; i++)
        kmer_mask[i] = keep_forward[i] | keep_reverse[i];

    ctx_free(keep_forward);
    ctx_free(keep_reverse);
    ctx_free(seen_forward);
    ctx_free(seen_reverse);
    ctx_free(stack);
    traversal_subgraph_builder_dealloc(&builder);
}
