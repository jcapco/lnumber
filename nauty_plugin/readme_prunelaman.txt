The following outputs the number of Laman graphs with 10 vertices.
geng 10 -K2 -u

-K2 implies that tightkn = 2 in the code

static long tightkn = 2; /* Specifies k for (k,l)-tight graphs. */
static long tightkd = 1;
static long tightln = 3; /* Specifies l for (k,l)-tight graphs. */
static long tightld = 1;

(tightkd == 1 && tightld == 1 && tightln >= 0 && tightln < 2 * tightkn) is satisfied
therefore prune = prunetightpebble

-u : do not output any graphs, just generate and count them