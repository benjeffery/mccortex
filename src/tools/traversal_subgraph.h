#ifndef TRAVERSAL_SUBGRAPH_H_
#define TRAVERSAL_SUBGRAPH_H_

#include "db_graph.h"
#include "seq_file/seq_file.h"

void traverse_and_mark(dBGraph *db_graph, size_t nthreads, size_t max_depth,
                         bool invert,
                         size_t stack_mem, uint8_t *kmer_mask,
                         seq_file_t **startfiles, size_t num_start_files,
                         seq_file_t **endfiles, size_t num_end_files
                         );


#endif /* TRAVERSAL_SUBGRAPH_H_ */
