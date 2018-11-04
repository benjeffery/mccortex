#include "global.h"
#include "subgraph.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "prune_nodes.h"
#include "seq_loading_stats.h"
#include "util.h"

typedef struct
{
    const dBGraph *const db_graph; // graph we are operating on
    uint8_t *const kmer_mask; // bitset of visited kmers
    SeqLoadingStats stats;
    BinaryKmer start_kmer, end_kmer;
} TraversalSubgraphBuilder;

static void traversal_subgraph_builder_alloc(TraversalSubgraphBuilder *builder,
                                             uint8_t *kmer_mask,
                                             const dBGraph *graph)
{
    TraversalSubgraphBuilder tmp = {.db_graph = graph,
                                    .kmer_mask = kmer_mask,
                                    };

    memcpy(builder, &tmp, sizeof(TraversalSubgraphBuilder));
    seq_loading_stats_init(&builder->stats);
}

static void traversal_subgraph_builder_dealloc(TraversalSubgraphBuilder *builder)
{
//    db_node_buf_dealloc(&builder->nbufs[0]);
//    db_node_buf_dealloc(&builder->nbufs[1]);
//    db_node_buf_dealloc(&builder->unitig_buf);
}


static void get_start_from_read(read_t *r1, read_t *r2,
                             uint8_t qoffset1, uint8_t qoffset2, void *ptr) {
    TraversalSubgraphBuilder *builder = (TraversalSubgraphBuilder*)ptr;
    const dBGraph *db_graph = builder->db_graph;
    if(r1->seq.end > db_graph->kmer_size) {
        die("Start sequence is shorter than kmer size");
    }

    //Cut off the string at kmer size
    r1->seq.b[db_graph->kmer_size] = "\0";
    builder->start_kmer = binary_kmer_from_str(r1->seq.b, db_graph->kmer_size);
}

static void get_end_from_read(read_t *r1, read_t *r2,
                                uint8_t qoffset1, uint8_t qoffset2, void *ptr) {
    TraversalSubgraphBuilder *builder = (TraversalSubgraphBuilder*)ptr;
    const dBGraph *db_graph = builder->db_graph;
    if(r1->seq.end > db_graph->kmer_size) {
        die("End sequence is shorter than kmer size");
    }
    builder->end_kmer = binary_kmer_from_str(r1->seq.b + r1->seq.end - db_graph->kmer_size,
                                             db_graph->kmer_size);
}

void traverse_and_mark(dBGraph *db_graph, size_t nthreads, size_t max_depth,
                       bool invert,
                       size_t stack_mem, uint8_t *kmer_mask,
                       seq_file_t **startfiles, size_t num_start_files,
                       seq_file_t **endfiles, size_t num_end_files)
{
    TraversalSubgraphBuilder builder;
    traversal_subgraph_builder_alloc(&builder, kmer_mask, db_graph);

    read_t r1;
    if(seq_read_alloc(&r1) == NULL)
        die("Out of memory");

    seq_parse_se_sf(startfiles[0], 0, &r1, get_start_from_read, &builder);
    seq_parse_se_sf(endfiles[0], 0, &r1, get_end_from_read, &builder);
    seq_read_dealloc(&r1);

    //As a test mark these two kmers
    dBNode node = db_graph_find(db_graph, builder.start_kmer);
    if(node.key != HASH_NOT_FOUND) {
        bitset_set(kmer_mask, node.key);
    }
    node = db_graph_find(db_graph, builder.end_kmer);
    if(node.key != HASH_NOT_FOUND) {
        bitset_set(kmer_mask, node.key);
    }

    traversal_subgraph_builder_dealloc(&builder);
}
