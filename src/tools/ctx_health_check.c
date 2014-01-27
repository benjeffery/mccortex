#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "graph_file_filter.h"
#include "path_format.h"
#include "path_file_filter.h"
#include "graph_paths.h"

static const char usage[] =
"usage: "CMD" healthcheck [options] <graph.ctx>\n"
"  Load a graph into memory to check it is valid.\n"
"\n"
"  Options:\n"
"   -m <mem>       Memory to use\n"
"   -n <kmers>     Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   -p <in.ctp>    Load paths\n"
"   --noedgecheck  Don't check kmer edges\n";


int ctx_health_check(CmdArgs *args)
{
  cmd_accept_options(args, "pmn", usage);
  int argi, argc = args->argc;
  char **argv = args->argv;
  if(argc != 1) print_usage(usage, NULL);

  size_t i;
  boolean do_edge_check = true;

  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++) {
    if(strcmp(argv[argi],"--noedgecheck")) do_edge_check = false;
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  if(argi+1 < argc) print_usage(usage, "Too many arguments");
  if(argi+1 > argc) print_usage(usage, "Too few arguments");

  char *ctx_path = argv[argi];

  if(!do_edge_check && args->num_ctp_files == 0) {
    print_usage(usage, "--noedgecheck and no path files (-p in.ctp) "
                       "- nothing to check.");
  }

  //
  // Open Graph file
  //
  GraphFileReader file = INIT_GRAPH_READER;
  graph_file_open(&file, ctx_path, true); // true => errors are fatal
  size_t ncols = graph_file_outncols(&file);

  //
  // Open path files
  //
  size_t num_pfiles = args->num_ctp_files;
  PathFileReader pfiles[num_pfiles];
  size_t path_max_mem = 0, path_max_usedcols = 0;

  for(i = 0; i < num_pfiles; i++) {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], args->ctp_files[i], true);
    path_max_mem = MAX2(path_max_mem, pfiles[i].hdr.num_path_bytes);
    path_max_usedcols = MAX2(path_max_usedcols, path_file_usedcols(&pfiles[i]));
  }

  // Decide on memory
  size_t extra_bits_per_kmer, kmers_in_hash, graph_mem;
  size_t path_mem = 0, tmppathsize = 0, total_mem;

  extra_bits_per_kmer = sizeof(Edges) * ncols * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        file.hdr.num_of_kmers, false, &graph_mem);

  // Path Memory
  if(num_pfiles) {
    tmppathsize = paths_merge_needs_tmp(pfiles, num_pfiles) ? path_max_mem : 0;
    path_mem = path_max_mem + tmppathsize;
  }

  total_mem = path_mem + graph_mem;
  cmd_check_mem_limit(args, graph_mem);

  // Create db_graph
  dBGraph db_graph;
  size_t num_edge_cols = do_edge_check ? ncols : 1;
  db_graph_alloc(&db_graph, file.hdr.kmer_size, ncols, num_edge_cols, kmers_in_hash);

  // Only need one edge per colour if doing edge check
  if(do_edge_check) {
    db_graph.col_edges = calloc2(db_graph.ht.capacity * ncols, sizeof(Edges));
  } else {
    size_t nwords = roundup_bits2bytes(db_graph.ht.capacity)*ncols;
    db_graph.node_in_cols = calloc2(nwords, sizeof(uint64_t));
    db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  }

  // Paths
  uint8_t *path_store = NULL, *tmp_pdata = NULL;

  if(num_pfiles > 0) {
    db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
    memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

    path_store = malloc2(path_max_mem);
    path_store_init(&db_graph.pdata, path_store, path_max_mem, path_max_usedcols);

    // Temorary memory to load paths into
    if(tmppathsize) tmp_pdata = malloc2(tmppathsize);
  }

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  graph_load(&file, gprefs, NULL);

  // Load path files
  if(num_pfiles) {
    paths_format_merge(pfiles, num_pfiles, false, tmp_pdata, tmppathsize, &db_graph);
    if(tmp_pdata) free(tmp_pdata);
  }

  if(do_edge_check) {
    status("Running edge check...");
    db_graph_healthcheck(&db_graph);
  }

  if(num_pfiles) {
    status("Running path check...");
    graph_paths_check_all_paths(&db_graph);
  }

  status("All looks good!");

  if(num_pfiles) {
    free((void*)db_graph.kmer_paths);
    free(path_store);
  }
  if(db_graph.node_in_cols) free(db_graph.node_in_cols);
  free(db_graph.col_edges);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
