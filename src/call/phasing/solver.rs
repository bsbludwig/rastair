//! Splitting a segment's sites across two haplotypes.
//!
//! The objective is minimum error correction: every fragment is charged against
//! whichever haplotype it fits worse. The optimum is NP-hard, so a maximum
//! spanning forest gives the initial orientation and single-site flips refine
//! it. The solver decides one bit per site: which of the site's two alleles
//! haplotype 1 carries.

use super::problem::{AlleleSupport, PhaseProblem, SiteAllele, SiteIndex, SiteMap, Sites};
use rustc_hash::FxHashMap;
use seqair_types::SmallVec;
use std::cmp::Reverse;
use tracing::warn;

/// Improving flips [`refine`] may make per site of a block before it gives up.
///
/// A work cap, not an oscillation guard: the descent only accepts strictly
/// lower costs, so it terminates on its own. Correcting the spanning forest
/// takes about one flip per wrongly oriented site, so the cap scales with the
/// block; the multiplier leaves room for a site to be flipped back a few times
/// as its neighbours move, while bounding each block's work to quadratic in its
/// size.
const MAX_FLIPS_PER_SITE: usize = 4;

/// Sites a chain of fragments tied together, and how they split: each site
/// and the allele haplotype 1 carries there, ascending by site.
///
/// Which haplotype is "1" is arbitrary, so the orientation is only defined up
/// to a flip of the whole block; the caller fixes it on the site it anchors
/// the block to.
pub type PhaseBlock = Vec<(SiteIndex, SiteAllele)>;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PhaseStats {
    /// Heterozygous columns the solver was offered.
    pub sites: usize,
    /// How many of them ended up in a block.
    pub phased_sites: usize,
    pub blocks: usize,
    pub largest_block: usize,
    /// Total minimum-error-correction cost of the assignment.
    pub mec_cost: u32,
}

#[derive(Debug, Clone, Default)]
pub struct Solution {
    pub blocks: Vec<PhaseBlock>,
    pub stats: PhaseStats,
}

/// Assign every site to a haplotype and group the connected ones into blocks.
#[must_use]
pub fn solve(problem: &PhaseProblem) -> Solution {
    let sites = &problem.sites;
    let Forest { blocks, mut haplotype, component } = spanning_forest(sites, &build_edges(problem));
    let mut restrictions = restrictions(problem, &component, &haplotype);
    let by_site = restrictions_by_site(sites, &restrictions);

    for block in &blocks {
        refine(block, &mut restrictions, &by_site, &mut haplotype);
    }

    let stats = PhaseStats {
        sites: sites.len(),
        phased_sites: blocks.iter().map(Vec::len).sum(),
        blocks: blocks.len(),
        largest_block: blocks.iter().map(Vec::len).max().unwrap_or(0),
        mec_cost: restrictions.iter().map(|r| mec(r.agree, r.total).unsigned_abs()).sum(),
    };
    let blocks = blocks
        .iter()
        .map(|block| block.iter().map(|&site| (site, haplotype[site])).collect())
        .collect();

    Solution { blocks, stats }
}

/// Evidence that two sites' [`SiteAllele::First`] alleles share a chromosome.
#[derive(Debug, Clone, Copy)]
struct Edge {
    a: SiteIndex,
    b: SiteIndex,
    /// Positive for "same orientation", negative for "opposite"; the magnitude
    /// is what an assignment pays for going against it.
    weight: i32,
}

/// Every site pair a fragment spans, with the votes summed, strongest first.
///
/// A fragment carrying the same allele slot at both sites votes for them
/// sharing an orientation, and the vote is worth the weaker of its two
/// observations — a link is only as good as its worse end.
///
/// The order is [`spanning_forest`]'s: strongest edge first, ties broken by the
/// site pair so that it does not depend on the hasher.
fn build_edges(problem: &PhaseProblem) -> Vec<Edge> {
    let mut votes: FxHashMap<(SiteIndex, SiteIndex), i32> = FxHashMap::default();

    for fragment in &problem.fragments {
        for (offset, first) in fragment.supports.iter().enumerate() {
            for second in fragment.supports.iter().skip(offset + 1) {
                let vote = i32::from(first.cost.min(second.cost));
                *votes.entry((first.site, second.site)).or_default() +=
                    if first.allele == second.allele { vote } else { -vote };
            }
        }
    }

    // A pair whose cis and trans evidence cancelled exactly says nothing, and
    // joining a block on it would emit a coin toss as a confident phase.
    let mut edges: Vec<Edge> = votes
        .into_iter()
        .filter(|&(_, weight)| weight != 0)
        .map(|((a, b), weight)| Edge { a, b, weight })
        .collect();
    edges.sort_unstable_by_key(|edge| (Reverse(edge.weight.unsigned_abs()), edge.a, edge.b));
    edges
}

/// Connected components and an initial orientation, from a maximum spanning
/// forest.
///
/// Edges arrive strongest first and one whose endpoints are already linked is
/// skipped, which is Kruskal's algorithm; the union-find carries the
/// orientation along, so walking the tree afterwards is unnecessary.
fn spanning_forest(sites: &Sites, edges: &[Edge]) -> Forest {
    let mut orientations = Orientations::new(sites);
    for edge in edges {
        orientations.link(edge.a, edge.b, edge.weight > 0);
    }

    let mut haplotype = SiteMap::filled(sites, SiteAllele::First);
    let mut component = SiteMap::from_fn(sites, |site| site);
    // Keyed by root, in order of first site, so the blocks come out that way.
    let mut block_of_root: FxHashMap<SiteIndex, usize> = FxHashMap::default();
    let mut blocks: Vec<Vec<SiteIndex>> = Vec::new();
    for site in sites.indices() {
        let (root, flipped) = orientations.find(site);
        if flipped {
            haplotype[site] = SiteAllele::Second;
        }
        component[site] = root;
        let next = blocks.len();
        match blocks.get_mut(*block_of_root.entry(root).or_insert(next)) {
            Some(block) => block.push(site),
            None => blocks.push(vec![site]),
        }
    }
    // A site no edge reached is its own component; it cannot be phased.
    blocks.retain(|block| block.len() > 1);
    Forest { blocks, haplotype, component }
}

/// The spanning forest's answer, before refinement.
struct Forest {
    /// Components of two or more sites, each ascending, ordered by first site.
    blocks: Vec<Vec<SiteIndex>>,
    /// Per site, the allele haplotype 1 carries.
    haplotype: SiteMap<SiteAllele>,
    /// Per site, its component's root: a key telling components apart, and
    /// defined for the unphasable singletons too.
    component: SiteMap<SiteIndex>,
}

/// One site's place in the union-find: its parent, whether it is oriented
/// against that parent, and the rank that keeps the trees shallow.
#[derive(Debug, Clone, Copy)]
struct Node {
    parent: SiteIndex,
    /// Orientation relative to `parent`: `true` means flipped against it.
    flipped: bool,
    rank: u32,
}

/// Union-find that also tracks whether a site is oriented like its component's
/// root or against it.
#[derive(Debug)]
struct Orientations(SiteMap<Node>);

impl Orientations {
    fn new(sites: &Sites) -> Self {
        // Every site starts as its own root, oriented with itself.
        Self(SiteMap::from_fn(sites, |site| Node { parent: site, flipped: false, rank: 0 }))
    }

    /// The site's component root and its orientation relative to it.
    ///
    /// No path compression: union by rank already keeps the walk logarithmic,
    /// and a segment holds at most a few thousand sites.
    fn find(&self, site: SiteIndex) -> (SiteIndex, bool) {
        let mut root = site;
        let mut flipped = false;
        while self.0[root].parent != root {
            flipped ^= self.0[root].flipped;
            root = self.0[root].parent;
        }
        (root, flipped)
    }

    /// Join two sites, `same` saying whether they take the same orientation.
    ///
    /// Sites already linked are left alone, which is what makes the accepted
    /// edges a maximum spanning forest, given they arrive strongest first.
    fn link(&mut self, a: SiteIndex, b: SiteIndex, same: bool) {
        let (root_a, flipped_a) = self.find(a);
        let (root_b, flipped_b) = self.find(b);
        if root_a == root_b {
            return;
        }

        // Orientation of the absorbed root relative to the absorbing one, such
        // that `a` and `b` end up as `same` says. Symmetric in the two roots,
        // so it holds whichever way round the rank comparison goes.
        let relative = flipped_a ^ flipped_b ^ !same;
        let (rank_a, rank_b) = (self.0[root_a].rank, self.0[root_b].rank);
        let (child, parent) = if rank_a < rank_b { (root_a, root_b) } else { (root_b, root_a) };

        self.0[child].parent = parent;
        self.0[child].flipped = relative;
        if rank_a == rank_b {
            self.0[parent].rank = self.0[parent].rank.saturating_add(1);
        }
    }
}

/// One fragment's supports inside one block: what the objective sums and what a
/// flip has to re-charge.
///
/// A fragment reaching two blocks — possible when the edge between them
/// cancelled out — is restricted to each separately, since its haplotype in one
/// says nothing about its haplotype in the other.
#[derive(Debug)]
struct Restriction {
    supports: SmallVec<AlleleSupport, 4>,
    /// Cost this fragment currently matches haplotype 1 on. Charging it to
    /// haplotype 2 costs exactly this, and to haplotype 1 the rest.
    agree: i32,
    total: i32,
}

const fn mec(agree: i32, total: i32) -> i32 {
    let disagree = total - agree;
    if agree < disagree { agree } else { disagree }
}

fn restrictions(
    problem: &PhaseProblem,
    component: &SiteMap<SiteIndex>,
    haplotype: &SiteMap<SiteAllele>,
) -> Vec<Restriction> {
    let mut restrictions = Vec::new();

    for fragment in &problem.fragments {
        let mut groups: SmallVec<(SiteIndex, SmallVec<AlleleSupport, 4>), 2> = SmallVec::new();
        for &support in &fragment.supports {
            let root = component[support.site];
            match groups.iter_mut().find(|(candidate, _)| *candidate == root) {
                Some((_, supports)) => supports.push(support),
                None => {
                    let mut supports = SmallVec::new();
                    supports.push(support);
                    groups.push((root, supports));
                }
            }
        }

        // One support inside a component charges nothing whichever haplotype
        // the fragment takes, so it cannot inform a flip either. This is also
        // what leaves out the singletons, which are in no block.
        for (_, supports) in groups.into_iter().filter(|(_, supports)| supports.len() > 1) {
            let total = supports.iter().map(|support| i32::from(support.cost)).sum();
            let agree = supports
                .iter()
                .filter(|support| haplotype[support.site] == support.allele)
                .map(|support| i32::from(support.cost))
                .sum();
            restrictions.push(Restriction { supports, agree, total });
        }
    }

    restrictions
}

fn restrictions_by_site(sites: &Sites, restrictions: &[Restriction]) -> SiteMap<Vec<usize>> {
    let mut by_site = SiteMap::filled(sites, Vec::new());
    for (index, restriction) in restrictions.iter().enumerate() {
        for support in &restriction.supports {
            by_site[support.site].push(index);
        }
    }
    by_site
}

/// What a restriction would agree on if this one site flipped, or `None` when
/// it does not reach the site.
fn agree_after_flip(
    restriction: &Restriction,
    site: SiteIndex,
    haplotype: &SiteMap<SiteAllele>,
) -> Option<i32> {
    let support = restriction.supports.iter().find(|support| support.site == site)?;
    let cost = i32::from(support.cost);
    let agrees = haplotype[site] == support.allele;
    Some(if agrees { restriction.agree - cost } else { restriction.agree + cost })
}

/// Flip single sites while that lowers the block's cost.
///
/// One flip per pass, always the most profitable one: a greedy descent on an
/// objective the spanning forest only approximated, since the forest fixes each
/// site against one neighbour and ignores every other edge. Runs to convergence
/// unless it reaches [`MAX_FLIPS_PER_SITE`] flips per site.
fn refine(
    block: &[SiteIndex],
    restrictions: &mut [Restriction],
    by_site: &SiteMap<Vec<usize>>,
    haplotype: &mut SiteMap<SiteAllele>,
) {
    let max_flips = MAX_FLIPS_PER_SITE * block.len();
    for _ in 0..max_flips {
        let mut best: Option<(SiteIndex, i32)> = None;
        for &site in block {
            let delta = flip_delta(site, restrictions, by_site, haplotype);
            if delta < 0 && best.is_none_or(|(_, previous)| delta < previous) {
                best = Some((site, delta));
            }
        }
        let Some((site, _)) = best else {
            return;
        };

        for &index in &by_site[site] {
            let Some(restriction) = restrictions.get_mut(index) else {
                continue;
            };
            if let Some(agree) = agree_after_flip(restriction, site, haplotype) {
                restriction.agree = agree;
            }
        }
        haplotype[site] = haplotype[site].flipped();
    }

    // Only reachable by falling out of the loop: converging returns above.
    warn!(
        sites = block.len(),
        flips = max_flips,
        "Phase block refinement hit its flip cap with an improving flip still available; \
         the block's phasing may be worse than it could be"
    );
}

fn flip_delta(
    site: SiteIndex,
    restrictions: &[Restriction],
    by_site: &SiteMap<Vec<usize>>,
    haplotype: &SiteMap<SiteAllele>,
) -> i32 {
    by_site[site]
        .iter()
        .filter_map(|&index| {
            let restriction = restrictions.get(index)?;
            let agree = agree_after_flip(restriction, site, haplotype)?;
            Some(mec(agree, restriction.total) - mec(restriction.agree, restriction.total))
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::phasing::{
        Cost, FragmentId, GtAllele,
        fixtures::ALT_1,
        problem::{Fragment, HetAllele, PhaseSite},
    };
    use SiteAllele::{First, Second};
    use seqair_types::{Base::*, Pos0};

    /// Sites are interchangeable to the solver — only their count and their
    /// positions matter — so the fixtures name them by index.
    fn sites(count: u32) -> Vec<PhaseSite> {
        (0..count)
            .map(|index| PhaseSite {
                pileup: index as usize,
                pos: Pos0::new(1000 + index * 10).expect("a small position"),
                alleles: [
                    HetAllele { allele: GtAllele::Ref, base: A },
                    HetAllele { allele: GtAllele::Alt(ALT_1), base: C },
                ],
                methylable: None,
            })
            .collect()
    }

    fn fragment(id: u64, supports: &[(u32, SiteAllele, u8)]) -> Fragment {
        Fragment {
            id: FragmentId::new(id).expect("non-zero test fragment id"),
            supports: supports
                .iter()
                .map(|&(site, allele, cost)| AlleleSupport {
                    site: SiteIndex::new(site),
                    allele,
                    cost: Cost(cost),
                })
                .collect(),
        }
    }

    /// A fragment on haplotype `flipped`, clean, over the given sites.
    fn clean(id: u64, truth: &[SiteAllele], span: std::ops::Range<u32>, flipped: bool) -> Fragment {
        let supports: Vec<_> = span
            .map(|site| {
                let allele = truth[site as usize];
                (site, if flipped { allele.flipped() } else { allele }, 40)
            })
            .collect();
        fragment(id, &supports)
    }

    fn solve_sites(count: u32, fragments: Vec<Fragment>) -> Solution {
        solve(&PhaseProblem { sites: sites(count).into(), fragments })
    }

    /// Each block oriented so its first site is `First`: the solver's answer
    /// is only defined up to a flip of the whole block.
    fn assignment(solution: &Solution) -> Vec<Vec<(usize, SiteAllele)>> {
        solution
            .blocks
            .iter()
            .map(|block| {
                let flip = block.first().is_some_and(|&(_, allele)| allele == Second);
                block.iter().map(|&(s, a)| (s.get(), if flip { a.flipped() } else { a })).collect()
            })
            .collect()
    }

    #[test]
    fn a_chain_of_four_sites_becomes_one_block() {
        // Three overlapping pairs, all cis, so every site takes the same
        // orientation as its neighbour.
        let solution = solve_sites(
            4,
            vec![
                fragment(1, &[(0, First, 40), (1, First, 40)]),
                fragment(2, &[(1, First, 40), (2, First, 40)]),
                fragment(3, &[(2, First, 40), (3, First, 40)]),
            ],
        );

        assert_eq!(assignment(&solution), [[(0, First), (1, First), (2, First), (3, First)]]);
        assert_eq!(solution.stats.blocks, 1);
        assert_eq!(solution.stats.phased_sites, 4);
        assert_eq!(solution.stats.largest_block, 4);
        assert_eq!(solution.stats.mec_cost, 0, "no fragment contradicts the answer");
    }

    #[test]
    fn a_trans_link_puts_the_two_sites_on_opposite_haplotypes() {
        let solution = solve_sites(2, vec![fragment(1, &[(0, First, 40), (1, Second, 40)])]);

        assert_eq!(assignment(&solution), [[(0, First), (1, Second)]]);
    }

    #[test]
    fn a_lone_noisy_fragment_is_outvoted() {
        let solution = solve_sites(
            2,
            vec![
                fragment(1, &[(0, First, 40), (1, First, 40)]),
                fragment(2, &[(0, Second, 40), (1, Second, 40)]),
                fragment(3, &[(0, First, 40), (1, Second, 40)]),
            ],
        );

        assert_eq!(
            assignment(&solution),
            [[(0, First), (1, First)]],
            "two fragments say cis, one says trans"
        );
        assert_eq!(solution.stats.mec_cost, 40, "the noisy fragment is what MEC pays for");
    }

    #[test]
    fn groups_with_no_fragment_between_them_are_separate_blocks() {
        let solution = solve_sites(
            4,
            vec![
                fragment(1, &[(0, First, 40), (1, Second, 40)]),
                fragment(2, &[(2, First, 40), (3, First, 40)]),
            ],
        );

        assert_eq!(assignment(&solution), [[(0, First), (1, Second)], [(2, First), (3, First)]]);
        assert_eq!(solution.stats.blocks, 2);
        assert_eq!(solution.stats.largest_block, 2);
    }

    #[test]
    fn a_site_no_fragment_reaches_stays_unphased() {
        let solution = solve_sites(3, vec![fragment(1, &[(0, First, 40), (2, First, 40)])]);

        assert_eq!(assignment(&solution), [[(0, First), (2, First)]]);
        assert_eq!(solution.stats.sites, 3);
        assert_eq!(solution.stats.phased_sites, 2, "site 1 is in no block");
    }

    /// Cis and trans evidence of exactly equal weight is no evidence, and
    /// phasing on it would write out a coin toss.
    #[test]
    fn a_pair_whose_evidence_cancels_is_not_linked() {
        let solution = solve_sites(
            2,
            vec![
                fragment(1, &[(0, First, 40), (1, First, 40)]),
                fragment(2, &[(0, First, 40), (1, Second, 40)]),
            ],
        );

        assert!(solution.blocks.is_empty());
        assert_eq!(solution.stats.phased_sites, 0);
    }

    /// The spanning forest fixes each site against one neighbour and ignores
    /// every other edge, so it can leave a site on the wrong haplotype; the
    /// flip pass is what notices. Here the strongest edge (0–1) is the wrong
    /// one, outweighed by the two weaker edges that agree with each other.
    #[test]
    fn refinement_corrects_a_spanning_tree_mistake() {
        let solution = solve_sites(
            3,
            vec![
                fragment(1, &[(0, First, 40), (1, Second, 40)]),
                fragment(2, &[(1, First, 30), (2, First, 30)]),
                fragment(3, &[(0, First, 30), (2, First, 30)]),
                fragment(4, &[(0, First, 30), (2, First, 30)]),
                fragment(5, &[(1, First, 30), (2, First, 30)]),
            ],
        );

        assert_eq!(
            assignment(&solution),
            [[(0, First), (1, First), (2, First)]],
            "four agreeing fragments beat the single strong trans link"
        );
    }

    /// Each of 25 leaves hangs off a backbone by one strong trans fragment the
    /// spanning forest takes, against two weaker cis fragments that together
    /// outweigh it. Every leaf needs its own corrective flip, more than a fixed
    /// cap of 20 flips per block allowed.
    #[test]
    fn refinement_runs_to_convergence_on_a_large_block() {
        const LEAVES: u32 = 25;
        let mut fragments = Vec::new();
        let mut id = 1u64;
        let mut push = |supports: &[(u32, SiteAllele, u8)]| {
            fragments.push(fragment(id, supports));
            id += 1;
        };
        // The backbone, 0-1-2, is cis and so heavily supported that no flip of
        // a backbone site ever pays.
        for _ in 0..40 {
            push(&[(0, First, 60), (1, First, 60)]);
            push(&[(1, First, 60), (2, First, 60)]);
            push(&[(0, First, 60), (2, First, 60)]);
        }
        for leaf in 3..3 + LEAVES {
            push(&[(0, First, 59), (leaf, Second, 59)]);
            push(&[(1, First, 40), (leaf, First, 40)]);
            push(&[(2, First, 40), (leaf, First, 40)]);
        }

        let solution = solve_sites(3 + LEAVES, fragments);

        let expected: Vec<_> = (0..3 + LEAVES as usize).map(|site| (site, First)).collect();
        assert_eq!(assignment(&solution), [expected]);
        assert_eq!(solution.stats.mec_cost, 59 * LEAVES, "only the trans fragments are wrong");
    }

    mod proptests {
        use super::*;
        use proptest::prelude::*;

        /// One misread fragment: where it starts, which haplotype it came from,
        /// and which of its two observations is wrong.
        type Misread = (u32, bool, usize);
        type Scenario = (Vec<bool>, Vec<Misread>);

        /// Fragments spanning a window of adjacent sites, as read pairs do.
        fn scenario() -> impl Strategy<Value = Scenario> {
            (2u32..=6).prop_flat_map(|count| {
                let truth = proptest::collection::vec(any::<bool>(), count as usize);
                let noise =
                    proptest::collection::vec((0..count - 1, any::<bool>(), 0usize..2), 0..=3);
                (truth, noise)
            })
        }

        /// Up to eight sites and up to a dozen fragments over any of them, with
        /// any alleles and costs: no structure for the forest to get right by
        /// luck, so refinement has real work to do.
        fn arbitrary_problem() -> impl Strategy<Value = (u32, Vec<Fragment>)> {
            (2u32..=8)
                .prop_flat_map(|count| {
                    let sites: Vec<u32> = (0..count).collect();
                    let fragment = (
                        proptest::sample::subsequence(sites, 2..=4.min(count as usize)),
                        proptest::collection::vec((any::<bool>(), 1u8..=60), 4),
                    );
                    (Just(count), proptest::collection::vec(fragment, 1..=12))
                })
                .prop_map(|(count, raw)| {
                    let fragments = raw
                        .into_iter()
                        .zip(1u64..)
                        .map(|((sites, observations), id)| {
                            let supports: Vec<_> = sites
                                .into_iter()
                                .zip(observations)
                                .map(|(site, (second, cost))| {
                                    (site, if second { Second } else { First }, cost)
                                })
                                .collect();
                            fragment(id, &supports)
                        })
                        .collect();
                    (count, fragments)
                })
        }

        /// MEC cost of the blocks as returned, computed from scratch: each
        /// fragment pays, per block, the cheaper of the two haplotypes.
        fn mec_from_scratch(fragments: &[Fragment], blocks: &[PhaseBlock]) -> u32 {
            let placed: FxHashMap<usize, (usize, SiteAllele)> = blocks
                .iter()
                .enumerate()
                .flat_map(|(block, phased)| {
                    phased.iter().map(move |&(site, allele)| (site.get(), (block, allele)))
                })
                .collect();

            let mut cost = 0;
            for fragment in fragments {
                let mut per_block: FxHashMap<usize, (u32, u32)> = FxHashMap::default();
                for support in &fragment.supports {
                    if let Some(&(block, allele)) = placed.get(&support.site.get()) {
                        let (agree, total) = per_block.entry(block).or_default();
                        *total += u32::from(support.cost.0);
                        if allele == support.allele {
                            *agree += u32::from(support.cost.0);
                        }
                    }
                }
                cost +=
                    per_block.values().map(|&(agree, total)| agree.min(total - agree)).sum::<u32>();
            }
            cost
        }

        proptest! {
            /// The reported cost is the cost of the blocks reported, and no single
            /// site can be flipped to lower it: refinement ran to a local
            /// optimum and its incremental bookkeeping did not drift.
            #[test]
            fn refinement_leaves_no_improving_flip((count, fragments) in arbitrary_problem()) {
                let solution = solve_sites(count, fragments.clone());
                let cost = mec_from_scratch(&fragments, &solution.blocks);
                prop_assert_eq!(solution.stats.mec_cost, cost);

                for (block, phased) in solution.blocks.iter().enumerate() {
                    for flip in 0..phased.len() {
                        let mut flipped = solution.blocks.clone();
                        if let Some((_, allele)) =
                            flipped.get_mut(block).and_then(|b| b.get_mut(flip))
                        {
                            *allele = allele.flipped();
                        }
                        let after = mec_from_scratch(&fragments, &flipped);
                        prop_assert!(
                            after >= cost,
                            "flipping site {} of block {} lowers the cost from {} to {}",
                            flip, block, cost, after
                        );
                    }
                }
            }

            /// Three clean fragments per adjacent pair are enough to recover the
            /// truth up to a global flip, whatever a handful of misreads say.
            #[test]
            fn clean_majority_recovers_the_haplotypes((truth, noise) in scenario()) {
                let truth: Vec<SiteAllele> = truth
                    .into_iter()
                    .map(|flipped| if flipped { Second } else { First })
                    .collect();
                let count = u32::try_from(truth.len()).expect("at most six sites");

                let mut fragments = Vec::new();
                let mut id = 1u64;
                for start in 0..count - 1 {
                    for copy in 0..3 {
                        fragments.push(clean(id, &truth, start..start + 2, copy % 2 == 1));
                        id += 1;
                    }
                }
                for (start, flipped, corrupt) in noise {
                    let mut fragment = clean(id, &truth, start..start + 2, flipped);
                    id += 1;
                    if let Some(support) = fragment.supports.get_mut(corrupt) {
                        support.allele = support.allele.flipped();
                        // Weaker than the clean evidence, so a pile of noise on
                        // one pair still loses.
                        support.cost = Cost(20);
                    }
                    fragments.push(fragment);
                }

                let solution = solve_sites(count, fragments);
                let block = solution.blocks.first().expect("every adjacent pair is linked");
                prop_assert_eq!(block.len(), truth.len());

                let flipped = block
                    .first()
                    .is_some_and(|&(site, allele)| truth[site.get()] != allele);
                for &(site, allele) in block {
                    let expected =
                        if flipped { truth[site.get()].flipped() } else { truth[site.get()] };
                    prop_assert_eq!(allele, expected, "site {} is on the wrong haplotype", site.get());
                }
            }
        }
    }
}
