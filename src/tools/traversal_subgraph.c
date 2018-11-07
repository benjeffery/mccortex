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

    signed long pos = 0;
    StackEntry* entry;
    Edges edges;
    Nucleotide next_bases[4];
    dBNode node;
    int keep, seen, ended_keep = 0, ended_depth = 0, ended_seen = 0;
    stack[pos].current = 0;
    stack[pos].num = 1;
    stack[pos].nodes[0] = builder.start_node;
    bitset_set(builder.end_node.orient == FORWARD ? keep_forward : keep_reverse, builder.end_node.key);

    unsigned long visited = 0;
    while (pos >= 0){
        ++visited;
        if (visited % 100000000 == 0) status("Visited: %lu depth:%ld", visited, pos);

        entry = &stack[pos];
        //Have we exhausted the branches of the parent node?
        if (entry->current >= entry->num) {
            --pos;
            //Remove node from seen so it can be visited by other paths
//            dBNode del_node = stack[pos].nodes[stack[pos].current-1];
//            bitset_del(del_node.orient == FORWARD ? seen_forward : seen_reverse, del_node.key);
            continue;
        }

        node = entry->nodes[entry->current];
        keep = bitset_get(node.orient == FORWARD ? keep_forward : keep_reverse, node.key);
        seen = bitset_get(node.orient == FORWARD ? seen_forward : seen_reverse, node.key);

//        char k[61];
//        BinaryKmer bk = db_node_get_bkey(db_graph, node.key);
//        if (node.orient != builder.start_node.orient) bk = binary_kmer_reverse_complement(bk, 61);
//        binary_kmer_to_str(bk, 61, k);
//        printf("%zu %d %d %d %s\n", pos, entry->current, seen, keep, k);

        if (keep) {
            ++ended_keep;
            //We hit a node that we are already keeping - keep all our path!
            status("Reached a target at depth: %zu\n", pos);
            for (size_t i = pos - 1; i > 0; i--) {
                node = stack[i].nodes[stack[i].current - 1]; //-1 here as current is pointing to the unexplored branch
                if (bitset_get(node.orient == FORWARD ? keep_forward : keep_reverse, node.key)) {
                    break;
                }
                bitset_set(node.orient == FORWARD ? keep_forward : keep_reverse, node.key);
            }
        }
        if (seen) {
            ++ended_seen;
        }
        if (!seen && !keep) {
            //Fresh node - mark as seen and explore children
            bitset_set(node.orient == FORWARD ? seen_forward : seen_reverse, node.key);
            ++pos;
            stack[pos].current = 0;
            if (pos < max_depth) {
                edges = db_node_get_edges_union(db_graph, node.key);
                stack[pos].num = db_graph_next_nodes(
                        db_graph,
                        db_node_get_bkey(db_graph, node.key),
                        node.orient, //Forward relative to the starting kmer
                        edges, &stack[pos].nodes, next_bases);
            } else {
                stack[pos].num = 0;
                ++ended_depth;
            }
        }
        //Explore the next child next time we visit
        ++(entry->current);
    };
    status("Traversal ended reasons: %d seen, %d keep, %d depth", ended_seen, ended_keep, ended_depth);

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
