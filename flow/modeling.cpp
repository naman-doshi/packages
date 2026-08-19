// ============================================================================
//  FLOW MODELLING CHEAT-SHEET  --  how to turn a problem INTO a flow problem
// ----------------------------------------------------------------------------
//  Not runnable code; this is the reference you actually reach for. The hard
//  part of flow problems is never the algorithm (use dinic.cpp / mcmf.cpp) --
//  it is recognising the reduction. Below are the standard ones.
//
//  ---------------------------------------------------------------------------
//  1. MAX-FLOW == MIN-CUT  (Ford-Fulkerson / Menger)
//     The value of a maximum s->t flow equals the minimum total capacity of a
//     set of edges whose removal disconnects s from t. Use dinic.cpp; recover
//     the cut with minCutSide(s) / cutEdges(s).
//
//  2. VERTEX CAPACITIES
//     A node v may pass at most c units: SPLIT it. Replace v by v_in and v_out,
//     add v_in -> v_out with cap c. All edges into v go to v_in; all edges out
//     of v leave v_out. (Vertex-disjoint paths = unit vertex capacities.)
//
//  3. MULTIPLE SOURCES / SINKS
//     Super-source SS -> each real source (cap = its supply).
//     Each real sink -> super-sink TT (cap = its demand). Then max flow SS->TT.
//
//  4. EDGE-DISJOINT vs VERTEX-DISJOINT PATHS  (Menger)
//     Max number of edge-disjoint s->t paths  = max flow with every cap = 1.
//     Max number of vertex-disjoint s->t paths = same, but split vertices (2).
//
//  ---------------------------------------------------------------------------
//  BIPARTITE / MATCHING FAMILY (see hopcroftKarp.cpp, kuhn.cpp)
//
//  5. MAXIMUM BIPARTITE MATCHING = max flow, unit caps (S->L->R->T).
//
//  6. KONIG'S THEOREM (bipartite graphs only):
//        min vertex cover      = max matching
//        max independent set   = V - max matching
//        min edge cover        = V - max matching   (no isolated vertices)
//     hopcroftKarp.cpp::minVertexCover() recovers the actual sets.
//
//  7. MINIMUM PATH COVER of a DAG (cover all vertices with fewest vertex-
//     disjoint directed paths):
//        answer = n - (max bipartite matching of the "split" graph)
//     Build bipartite graph: left copy u_L, right copy v_R; add u_L -- v_R for
//     every DAG edge u->v. For paths that may share vertices, first take the
//     transitive closure, then do the same.
//
//  8. DILWORTH'S THEOREM: in a DAG/poset, the minimum number of chains covering
//     all elements = size of the maximum antichain. Min chain cover is the
//     min path cover on the transitive closure (see 7).
//
//  ---------------------------------------------------------------------------
//  MIN-CUT MODELLING FAMILY
//
//  9. PROJECT SELECTION / MAX-WEIGHT CLOSURE
//     Pick a subset of projects to MAXIMISE (sum of profits - sum of costs),
//     where choosing a project forces choosing its prerequisites.
//        - profit node p (value w>0): S -> p with cap w.
//        - cost node   c (value w>0): c -> T with cap w.
//        - prerequisite "p needs q": p -> q with cap INF (can't cut it).
//        answer = (sum of all positive profits) - (min cut).
//     Nodes on the S-side of the min cut are the chosen set.
//
//  10. MINIMUM CUT WITH "EITHER/OR" PENALTIES (project = 2 categories):
//      Assign each item to side A or B; pay penalties for individual choices
//      and for pairs split across sides. Model each item as a node between S
//      (=A) and T (=B); individual penalties are S/T edges, pairwise penalties
//      are edges between item nodes. Minimise total penalty = min cut.
//
//  11. IMAGE SEGMENTATION / labelling: identical structure to (10) -- each
//      pixel is a node, foreground/background likelihoods are S/T edges,
//      smoothness penalties are neighbour edges.
//      (10) and (11) with an ARBITRARY 2x2 cost table per pair, plus the
//      submodularity test that says whether a model is even legal, are done
//      for you in labeling.cpp.
//
//  ---------------------------------------------------------------------------
//  MIN-COST FLOW FAMILY (see mcmf.cpp, hungarian.cpp)
//
//  12. ASSIGNMENT PROBLEM: min-cost perfect matching of an n x n cost matrix.
//      Dense -> hungarian.cpp (O(n^3)); sparse/awkward -> mcmf.cpp.
//
//  13. TRANSPORTATION / b-MATCHING: supplies s_i, demands d_j, per-unit ship
//      cost c_ij. S -> supply_i (cap s_i, cost 0); supply_i -> demand_j
//      (cap INF, cost c_ij); demand_j -> T (cap d_j, cost 0). Then min-cost
//      max-flow.
//
//  14. MIN-COST TO SEND EXACTLY K UNITS: run mcmf but stop augmenting once the
//      accumulated flow reaches K (clamp the last augmentation).
//
//  15. DISJOINT PATHS OF MIN TOTAL COST (k paths): give each edge cap 1 (or
//      the allowed multiplicity), cost = length, cap the source at k, MCMF.
//
//  16. LOWER BOUNDS ON EDGES (each edge must carry >= L): see boundedFlow.cpp.
//
//  ---------------------------------------------------------------------------
//  MORE MIN-CUT REDUCTIONS
//
//  17. MAX-WEIGHT INDEPENDENT SET / MIN-WEIGHT VERTEX COVER, BIPARTITE
//      *** Konig (6) DOES NOT GENERALISE TO WEIGHTS. *** With weights this is
//      a min cut, NOT a maximum matching:
//        S -> left u (cap w_u),  u -> v (cap INF) per conflict,  right v -> T
//        (cap w_v).  min cover = min cut;  max independent set = total - cut.
//      Code + recovery of the actual sets: weightedCover.cpp.
//
//  18. GLOBAL MIN CUT (split the graph in two, NO given s and t)
//      Not a max flow -- there is no source to fix. Stoer-Wagner, O(V^3), in
//      globalMinCut.cpp. Unweighted edge connectivity = all weights 1.
//      (If s and t ARE given, it is an ordinary max flow. Do that instead.)
//
//  19. MAXIMUM DENSITY SUBGRAPH  (maximise |E_S| / |S|)
//      Binary search a guess lambda and ask "is there S with |E_S| - lambda|S|
//      > 0" -- which is exactly max closure (9) with each edge a project worth
//      1 and each vertex a cost of lambda. density.cpp does it in exact
//      integers. NOT max clique: density is an average, which is why it is
//      polynomial at all.
//
//  20. FEASIBILITY / "CAN X STILL WIN?" (baseball elimination)
//      Ask whether a max flow SATURATES the source. Give every remaining game
//      a node (S -> game, cap = games left) feeding its two teams (cap INF),
//      and team -> T with cap = how many more wins that team may take before
//      it passes X. X survives iff every source edge saturates; the teams on
//      the source side of the cut are the group that eliminates it.
//
//  21. SCHEDULING WITH RELEASE TIMES AND DEADLINES
//      Jobs need p_j units of work, each only inside [r_j, d_j]. Cut time into
//      elementary intervals at every distinct r/d. S -> job (cap p_j), job ->
//      interval (cap = interval length) for every interval the job may use,
//      interval -> T (cap = length * machines). Feasible iff the flow
//      saturates S. Preemption is assumed -- without it this is NP-hard.
//
//  22. VERTEX CONNECTIVITY (Menger, global)
//      Split every vertex (2), give the internal edges cap 1 and all others
//      INF, then min over non-adjacent pairs of a max flow. Fixing s and
//      trying all t is enough for edge connectivity; for vertex connectivity
//      you must also vary s over k+1 vertices.
//
//  23. "EXACTLY K" / CARDINALITY-CONSTRAINED CUTS
//      A plain min cut cannot count. Either Lagrangian-relax (add a penalty
//      lambda per chosen item and binary search lambda until the count comes
//      out right -- works when the count is monotone in lambda), or, if k is
//      small, route through a gadget node of capacity k.
//
//  24. MATRIX ROUNDING
//      Round every entry of a real matrix up or down so all row and column
//      sums stay their own rounded values: a flow with LOWER BOUNDS
//      (boundedFlow.cpp) -- row -> column edges bounded by [floor, ceil].
//
//  ---------------------------------------------------------------------------
//  MODELLING CHECKLIST when you suspect a flow problem:
//    * "match / assign / pair up X with Y, each used once" -> bipartite matching.
//    * "choose a subset, dependencies/penalties, maximise profit" -> min-cut.
//    * "route/ship as much/as cheaply as possible" -> max flow / min-cost flow.
//    * "fewest paths/chains to cover everything" -> path cover (matching).
//    * "each thing used between L and C times" -> lower bounds (boundedFlow).
//    * capacity on a NODE, not an edge -> split the node.
//    * "two kinds, penalty when neighbours disagree" -> labeling.cpp.
//    * "max-weight set with no conflicting pair", weights + bipartite ->
//      weightedCover.cpp. Matching is the UNWEIGHTED answer only.
//    * "split the graph in two", no s or t named -> globalMinCut.cpp.
//    * "densest / most tightly-knit subgroup" -> density.cpp.
//    * "is this still possible" -> feasibility: does the flow saturate S (20)?
//    * mcmf.cpp too slow on a BIG graph -> mcmfDijkstra.cpp (same number of
//      augmentations, O(E log V) per path instead of O(V*E)). If the FLOW
//      VALUE itself is the problem, neither helps -- remodel or cap it.
//    Keep N small: flow is polynomial but constants bite past ~1e5 edges.
// ============================================================================
int main() { return 0; }
