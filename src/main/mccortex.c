#include "global.h"
#include <time.h>

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "hash.h"

// To add a new command to mccortex31 <cmd>:
// 0. create a file src/commands/ctx_X.c
// 1. add function and usage declaration to src/commands/commands.h
// 2. add entry to cmdobjs below
// 3. that's it!  Write the function I guess...

typedef struct
{
  const char *cmd, *blurb, *usage, *optargs, *reqargs;
  int minargs, maxargs; // counts AFTER standard args taken
  int hide; // set hide to >0 to remove from listings
  int (*func)(int argc, char **argv);
} CtxCmd;

CtxCmd cmdobjs[] = {
{
  .cmd = "build", .func = ctx_build, .hide = false,
  .blurb = "construct cortex graph from FASTA/FASTQ/BAM",
  .usage = build_usage
},
{
  .cmd = "sort", .func = ctx_sort, .hide = false,
  .blurb = "sort the kmers in a graph file",
  .usage = sort_usage
},
{
  .cmd = "index", .func = ctx_index, .hide = false,
  .blurb = "index a sorted cortex graph file",
  .usage = index_usage
},
{
  .cmd = "view", .func = ctx_view, .hide = false,
  .blurb = "text view of a cortex graph file (.ctx)",
  .usage = view_usage
},
{
  .cmd = "pview", .func = ctx_pview, .hide = false,
  .blurb = "text view of a cortex link file (.ctp)",
  .usage = pview_usage
},
{
  .cmd = "check", .func = ctx_health_check, .hide = false,
  .blurb = "load and check graph (.ctx) and path (.ctp) files",
  .usage = health_usage
},
{
  .cmd = "clean", .func = ctx_clean, .hide = false,
  .blurb = "clean errors from a graph",
  .usage = clean_usage
},
{
  .cmd = "join", .func = ctx_join, .hide = false,
  .blurb = "combine graphs, filter graph intersections",
  .usage = join_usage
},
{
  .cmd = "unitigs", .func = ctx_unitigs, .hide = false,
  .blurb = "pull out unitigs in FASTA, DOT or GFA format",
  .usage = unitigs_usage
},
{
  .cmd = "subgraph", .func = ctx_subgraph, .hide = false,
  .blurb = "filter a subgraph using seed kmers",
  .usage = subgraph_usage
},
{
  .cmd = "reads", .func = ctx_reads, .hide = false,
  .blurb = "filter reads against a graph",
  .usage = reads_usage
},
{
  .cmd = "contigs", .func = ctx_contigs, .hide = false,
  .blurb = "assemble contigs for a sample",
  .usage = contigs_usage
},
{
  .cmd = "inferedges", .func = ctx_infer_edges, .hide = false,
  .blurb = "infer graph edges between kmers before calling `thread`",
  .usage = inferedges_usage
},
{
  .cmd = "thread", .func = ctx_thread, .hide = false,
  .blurb = "thread reads through cleaned graph to make links",
  .usage = thread_usage,
},
{
  .cmd = "correct", .func = ctx_correct, .hide = false,
  .blurb = "error correct reads",
  .usage = correct_usage
},
{
  .cmd = "pjoin", .func = ctx_pjoin, .hide = false,
  .blurb = "merge link files (.ctp)",
  .usage = pjoin_usage
},
{
  .cmd = "bubbles", .func = ctx_bubbles, .hide = false,
  .blurb = "find bubbles in graph which are potential variants",
  .usage = bubbles_usage
},
{
  .cmd = "breakpoints", .func = ctx_breakpoints, .hide = false,
  .blurb = "use a trusted assembled genome to call large events",
  .usage = breakpoints_usage
},
{
  .cmd = "coverage", .func = ctx_coverage, .hide = false,
  .blurb = "print contig coverage",
  .usage = coverage_usage
},
{
  .cmd = "rmsubstr", .func = ctx_rmsubstr, .hide = false,
  .blurb = "reduce set of strings to remove substrings",
  .usage = rmsubstr_usage
},
{
  .cmd = "uniqkmers", .func = ctx_uniqkmers, .hide = false,
  .blurb = "generate random unique kmers",
  .usage = uniqkmers_usage
},
{
  .cmd = "links", .func = ctx_links, .hide = false,
  .blurb = "clean and plot link files (.ctp)",
  .usage = links_usage
},
{
  .cmd = "popbubbles", .func = ctx_pop_bubbles, .hide = false,
  .blurb = "pop bubbles in the population graph",
  .usage = pop_bubbles_usage
},
{
  .cmd = "calls2vcf", .func = ctx_calls2vcf, .hide = false,
  .blurb = "convert bubble/breakpoint calls to VCF",
  .usage = calls2vcf_usage
},
{
  .cmd = "server", .func = ctx_server, .hide = false,
  .blurb = "interactively query the graph",
  .usage = server_usage
},
{
  .cmd = "dist", .func = ctx_dist_matrix, .hide = false,
  .blurb = "make colour kmer distance matrix",
  .usage = dist_matrix_usage
},
{
  .cmd = "vcfcov", .func = ctx_vcfcov, .hide = false,
  .blurb = "coverage of a VCF against cortex graphs",
  .usage = vcfcov_usage
},
{
  .cmd = "vcfgeno", .func = ctx_vcfgeno, .hide = false,
  .blurb = "genotype a VCF after running vcfcov",
  .usage = vcfgeno_usage
},
/* Experiments */
{
  .cmd = "exp_abc", .func = ctx_exp_abc, .hide = true,
  .blurb = "run experiment on traversal properties",
  .usage = exp_abc_usage
},
{
  .cmd = "hashtest", .func = ctx_exp_hashtest, .hide = true,
  .blurb = "test hash table speed",
  .usage = exp_hashtest_usage
},
{
  .cmd = "traversalsubgraph", .func = ctx_exp_traversal_subgraph, .hide = false,
  .blurb = "filter a graph to those kmers on paths between a start and end kmer",
  .usage = exp_traversal_subgraph_usage
}
};


//
// Command listing
//

static const char options[] =
"  Type a command with no arguments to see help.\n"
"\n"
"Common Options:\n"
"  -h, --help            Help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -f, --force           Overwrite output files if they already exist\n"
"  -m, --memory <M>      Memory e.g. 1GB [default: 1GB]\n"
"  -n, --nkmers <H>      Hash entries [default: 4M, ~4 million]\n"
"  -t, --threads <T>     Limit on proccessing threads [default: 2]\n"
"  -o, --out <file>      Output file\n"
"  -p, --paths <in.ctp>  Links file to load (can specify multiple times)\n"
"\n";

static int ctxcmd_cmp(const void *aa, const void *bb)
{
  const CtxCmd *a = (const CtxCmd*)aa, *b = (const CtxCmd*)bb;
  return strcmp(a->cmd, b->cmd);
}

static void print_help(FILE *out, const char *errfmt,  ...)
  __attribute__((noreturn))
  __attribute__((format(printf, 2, 3)));

static void print_help(FILE *out, const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    fprintf(out, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(out, errfmt, argptr);
    va_end(argptr);
  }

  size_t i, j, s, maxlen = 0, n = sizeof(cmdobjs) / sizeof(CtxCmd);

  qsort(cmdobjs, n, sizeof(CtxCmd), ctxcmd_cmp);

  for(i = 0; i < n; i++)
    if(!cmdobjs[i].hide) maxlen = MAX2(maxlen, strlen(cmdobjs[i].cmd));

  fprintf(out, "\n"
"usage: "CMD" <command> [options] <args>\n"
"version: "VERSION_STATUS_STR" k=%i..%i\n"
"\n", get_min_kmer_size(), get_max_kmer_size());

  fprintf(out, "Commands:   ");
  for(i = 0, j = 0; i < n; i++) {
    if(!cmdobjs[i].hide) {
      fprintf(out, "%s%s", j > 0 ? "            " : "", cmdobjs[i].cmd);
      s = maxlen - strlen(cmdobjs[i].cmd) + 2;
      while(--s != SIZE_MAX) fputc(' ', out);
      fprintf(out, "%s\n", cmdobjs[i].blurb);
      j++;
    }
  }

  fprintf(out, "\n%s", options);
  exit(EXIT_FAILURE);
}

static const CtxCmd* ctx_get_command(const char* cmd)
{
  size_t i, n = sizeof(cmdobjs) / sizeof(CtxCmd);
  for(i = 0; i < n; i++)
    if(strcasecmp(cmdobjs[i].cmd,cmd) == 0)
      return &cmdobjs[i];
  return NULL;
}

// remove --quiet and -q flags
// ['--quiet','-qf','-q','-aq'] -> ['-f','-a']
// returns true iff any quiet flags were found
static bool remove_quiet_flags(int *argcp, char **argv)
{
  bool qfound = false;
  int i, j, argc = *argcp;
  char *p, *q;
  for(i = j = 1; i < argc; i++) {
    if(strcmp(argv[i],"--quiet") == 0 || strcmp(argv[i],"-q") == 0) qfound = true;
    else {
      // Remove q from short option
      if(argv[i][0] == '-' && argv[i][1] != '-') {
        for(p = q = argv[i]+1; *p; p++) {
          if(*p == 'q') qfound = true;
          else *q++ = *p;
        }
        *q = '\0';
      }
      argv[j++] = argv[i];
    }
  }

  *argcp = j;
  return qfound;
}

int main(int argc, char **argv)
{
  time_t start, end;
  time(&start);

  ctx_msg_out = stderr;
  cortex_init();
  cmd_init(argc, argv);

  if(argc == 1) print_help(stderr, NULL);
  const CtxCmd *cmd = ctx_get_command(argv[1]);
  if(cmd == NULL) print_help(stderr, "Unrecognised command: %s", argv[1]);

  // Once we have set cmd_usage, we can call cmd_print_usage() from anywhere
  cmd_set_usage(cmd->usage);

  // If no arguments after command, print help
  if(argc == 2) cmd_print_usage(NULL);

  // Look for -q, --quiet argument, if given silence output
  if(remove_quiet_flags(&argc, argv)) { ctx_msg_out = NULL; }

  // Print status header
  cmd_print_status_header();

  SWAP(argv[1],argv[0]);
  int ret = cmd->func(argc-1, argv+1);

  time(&end);
  cmd_destroy();

  // Warn if more allocations than deallocations
  size_t still_alloced = alloc_get_num_allocs() - alloc_get_num_frees();
  if(still_alloced) warn("%zu allocates not free'd.", still_alloced);

  char nallocs_str[50];
  ulong_to_str(alloc_get_num_allocs(), nallocs_str);
  status("[memory] We made %s allocs", nallocs_str);

  status(ret == 0 ? "Done." : "Fail.");

  // Print time taken
  double diff = difftime(end,start);
  if(diff < 60) status("[time] %.2lf seconds\n", diff);
  else {
    char timestr[100];
    seconds_to_str((size_t)diff, timestr);
    status("[time] %.2lf seconds (%s)\n", diff, timestr);
  }

  cortex_destroy();

  return ret;
}
