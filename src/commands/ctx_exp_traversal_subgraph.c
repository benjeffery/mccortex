#include "global.h"
#include "commands.h"
#include "seq_reader.h"
#include "hash_table.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "file_util.h"
#include "graph_writer.h"
#include "db_node.h"
#include "traversal_subgraph.h"
#include "prune_nodes.h"


const char exp_traversal_subgraph_usage[] =
        "usage: "CMD" traversalsubgraph [options] <in.ctx>[:cols] [in2.ctx ...]\n"
        "\n"
        "  Loads graphs (in.ctx) and dumps a graph (out.ctx) that contains all kmers on routes\n"
        "  between kmers in <start.fa> and <end.fa>.  Maintains number of colours / covgs etc.\n"
        "\n"
        "  -h, --help             This help message\n"
        "  -q, --quiet            Silence status output normally printed to STDERR\n"
        "  -f, --force            Overwrite output files\n"
        "  -o, --out <out.ctx>    Save output graph file [required]\n"
        "  -m, --memory <mem>     Memory to use\n"
        "  -n, --nkmers <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
        "  -t, --threads <T>      Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
        //
        "  -s, --start <start.fa> Read in a start seed file [required]\n"
        "  -e, --end <end.fa>     Read in an end seed file [required]\n"
        "  -d, --depth <N>         Maximum search depth [required]\n"
        "  -v, --invert           Dump kmers not in subgraph\n"
        "\n";

static struct option longopts[] =
        {
// General options
                {"help",    no_argument,       NULL, 'h'},
                {"force",   no_argument,       NULL, 'f'},
                {"out",     required_argument, NULL, 'o'},
                {"memory",  required_argument, NULL, 'm'},
                {"nkmers",  required_argument, NULL, 'n'},
                {"threads", required_argument, NULL, 't'},
// command specific
                {"start",   required_argument, NULL, 's'},
                {"end",     required_argument, NULL, 'e'},
                {"depth",   required_argument, NULL, 'd'},
                {"invert",  no_argument,       NULL, 'v'},
                {NULL, 0,                      NULL, 0}
        };

int ctx_exp_traversal_subgraph(int argc, char **argv) {
    size_t nthreads = 0;
    struct MemArgs memargs = MEM_ARGS_INIT;
    const char *out_path = NULL;
    size_t i, j, max_depth = 3000;
    bool invert = false;

    seq_file_t *tmp_startfile;
    SeqFilePtrBuffer startfilebuf;
    seq_file_ptr_buf_alloc(&startfilebuf, 16);

    seq_file_t *tmp_endfile;
    SeqFilePtrBuffer endfilebuf;
    seq_file_ptr_buf_alloc(&endfilebuf, 16);

    // Arg parsing
    char cmd[100], shortopts[100];
    cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
    int c;

    while ((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
        cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
        switch (c) {
            case 0: /* flag set */ break;
            case 'h':
                cmd_print_usage(NULL);
                break;
            case 'f':
                cmd_check(!futil_get_force(), cmd);
                futil_set_force(true);
                break;
            case 'o':
                cmd_check(!out_path, cmd);
                out_path = optarg;
                break;
            case 't':
                cmd_check(!nthreads, cmd);
                nthreads = cmd_uint32_nonzero(cmd, optarg);
                break;
            case 'm':
                cmd_mem_args_set_memory(&memargs, optarg);
                break;
            case 'n':
                cmd_mem_args_set_nkmers(&memargs, optarg);
                break;
            case 's':
                if ((tmp_startfile = seq_open(optarg)) == NULL)
                    die("Cannot read --start file %s", optarg);
                seq_file_ptr_buf_add(&startfilebuf, tmp_startfile);
                break;
            case 'e':
                if ((tmp_endfile = seq_open(optarg)) == NULL)
                    die("Cannot read --end file %s", optarg);
                seq_file_ptr_buf_add(&endfilebuf, tmp_endfile);
                break;
            case 'd':
                cmd_check(!max_depth, cmd);
                max_depth = cmd_uint32(cmd, optarg);
                break;
            case 'v':
                cmd_check(!invert, cmd);
                invert = true;
                break;
            case ':': /* BADARG */
            case '?': /* BADCH getopt_long has already printed error */
                // cmd_print_usage(NULL);
                die("`"
                            CMD
                            " tranversalsubgraph -h` for help. Bad option: %s", argv[optind - 1]);
            default:
                abort();
        }
    }

    // Defaults
    if (nthreads == 0) nthreads = DEFAULT_NTHREADS;

    if (startfilebuf.len == 0) cmd_print_usage("Require at least one --start file");
    if (endfilebuf.len == 0) cmd_print_usage("Require at least one --end file");
    if (optind >= argc) cmd_print_usage("Require input graph files (.ctx)");

    size_t num_gfiles = (size_t) (argc - optind);
    char **gfile_paths = argv + optind;

    size_t total_cols;

    // Open graph files
    GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
    size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

    total_cols = graph_files_open(gfile_paths, gfiles, num_gfiles,
                                  &ctx_max_kmers, &ctx_sum_kmers);

    //
    // Decide on memory
    //
    size_t bits_per_kmer, kmers_in_hash, graph_mem;
    size_t num_of_stack_entries, stack_mem, total_mem;
    char graph_mem_str[100], stack_mem_str[100], num_of_stack_entries_str[100];

    bits_per_kmer = sizeof(BinaryKmer) * 8 +
                    ((sizeof(Edges) + sizeof(Covg)) * 8 + 1);

    kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                          memargs.mem_to_use_set,
                                          memargs.num_kmers,
                                          memargs.num_kmers_set,
                                          bits_per_kmer,
                                          ctx_max_kmers, ctx_sum_kmers,
                                          false, &graph_mem);

    graph_mem = hash_table_mem(kmers_in_hash, bits_per_kmer, NULL);
    bytes_to_str(graph_mem, 1, graph_mem_str);

    if (graph_mem >= memargs.mem_to_use)
        die("Not enough memory for graph (requires %s)", graph_mem_str);

    // Stack entries
    stack_mem = memargs.mem_to_use - graph_mem;
    num_of_stack_entries = stack_mem / (sizeof(size_t/*STACK RECORD*/) * 2);
    ulong_to_str(num_of_stack_entries, num_of_stack_entries_str);
    bytes_to_str(stack_mem, 1, stack_mem_str);

    status("[memory] stack entries: %s (%s)\n", stack_mem_str, num_of_stack_entries_str);

    if (stack_mem < max_depth)
        die("Not enough memory for the graph search to the specifed depth(set -m <mem> higher)");

    // Don't need to check, but it prints out memory
    total_mem = graph_mem + stack_mem;
    cmd_check_mem_limit(memargs.mem_to_use, total_mem);

    //
    // Open output file
    //

    // Print to stdout unless --out <out> is specified
    if (out_path == NULL) out_path = "-";
    futil_create_output(out_path);

    // Create db_graph
    // multiple colours may be useful later in pulling out multiple colours
    dBGraph db_graph;
    db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, 1, 1,
                   kmers_in_hash, DBG_ALLOC_EDGES | DBG_ALLOC_COVGS);

    uint8_t *kmer_mask = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

    //
    // Load graphs
    //
    GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);

    StrBuf intersect_gname;
    strbuf_alloc(&intersect_gname, 1024);

    if (total_cols > db_graph.num_of_cols) {
        graphs_load_files_flat(gfiles, num_gfiles, gprefs, NULL);
    } else {
        for (i = 0; i < num_gfiles; i++)
            graph_load(&gfiles[i], gprefs, NULL);
    }

    // Create header
    for (i = 0; i < num_gfiles; i++) {
        for (j = 0; j < file_filter_num(&gfiles[i].fltr); j++) {
            size_t fromcol = file_filter_fromcol(&gfiles[i].fltr, j);
            graph_info_make_intersect(&gfiles[i].hdr.ginfo[fromcol], &intersect_gname);
        }
    }

    hash_table_print_stats(&db_graph.ht);

    char subgraphstr[] = "subgraph:{";
    strbuf_insert(&intersect_gname, 0, subgraphstr, strlen(subgraphstr));
    strbuf_append_char(&intersect_gname, '}');

    traverse_and_mark(&db_graph, nthreads, max_depth,
                        invert,
                        stack_mem, kmer_mask,
                        startfilebuf.b, startfilebuf.len,
                        endfilebuf.b, endfilebuf.len);

    if(invert) {
        status("Inverting selection...");
        size_t num_bytes = roundup_bits2bytes(db_graph.ht.capacity);
        for(i = 0; i < num_bytes; i++)
            kmer_mask[i] = ~kmer_mask[i];
    }

    status("Pruning untouched nodes...");
    prune_nodes_lacking_flag(nthreads, kmer_mask, &db_graph);


    for (i = 0; i < startfilebuf.len; i++) seq_close(startfilebuf.b[i]);
    seq_file_ptr_buf_dealloc(&startfilebuf);
    for (i = 0; i < endfilebuf.len; i++) seq_close(endfilebuf.b[i]);
    seq_file_ptr_buf_dealloc(&endfilebuf);

    ctx_free(kmer_mask);
    hash_table_print_stats(&db_graph.ht);

    // Dump nodes that were flagged
    Edges *intersect_edges = NULL;
    bool kmers_loaded = true;
    bool colours_loaded = (total_cols <= db_graph.num_of_cols);

    if (!colours_loaded) {
        // Need to reload graph colours - therefore construct edge intersection set
        intersect_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));
        for (i = 0; i < db_graph.ht.capacity; i++)
            intersect_edges[i] = db_node_get_edges_union(&db_graph, i);
    }

    graph_writer_merge_mkhdr(out_path, gfiles, num_gfiles,
                             kmers_loaded, colours_loaded,
                             intersect_edges, intersect_gname.b,
                             false, &db_graph);

    ctx_free(intersect_edges);
    strbuf_dealloc(&intersect_gname);
    for (i = 0; i < num_gfiles; i++) graph_file_close(&gfiles[i]);
    ctx_free(gfiles);

    db_graph_dealloc(&db_graph);

    return EXIT_SUCCESS;
}

