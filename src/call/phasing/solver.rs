//! Splitting a segment's sites across two haplotypes.
//!
//! The objective is minimum error correction: every fragment is charged against
//! whichever haplotype it fits worse. The optimum is NP-hard, so a maximum
//! spanning forest gives the initial orientation and single-site flips refine
//! it. The solver decides one bit per site: which of the site's two alleles
//! haplotype 1 carries.

use super::{
    Cost,
    problem::{AlleleSupport, PhaseProblem, SiteAllele, SiteIndex, SiteMap, Sites},
};
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

/// `--phase-min-switch-quality`: below it a solved block is cut in two. The
/// quality is a staircase of fragment counts (one fragment 26 at most, two
/// about 50, three 72 and up), so every threshold from 27 to 49 asks for two
/// fragments: on NA12878 chr12 a fifth of the switch errors of 26, for a
/// sixth fewer phased pairs.
pub const DEFAULT_MIN_SWITCH_QUALITY: u8 = 30;

/// Sites a chain of fragments tied together, and how they split, ascending by
/// site.
///
/// Which haplotype is "1" is arbitrary, so the orientation is only defined up
/// to a flip of the whole block; the caller fixes it on the site it anchors
/// the block to.
pub type PhaseBlock = Vec<PhasedSite>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhasedSite {
    pub site: SiteIndex,
    /// The allele haplotype 1 carries.
    pub allele: SiteAllele,
    /// How sure the block is of this site's orientation against its
    /// neighbours, VCF `PQ`.
    pub quality: Cost,
}

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

/// Assign every site to a haplotype, group the connected ones into blocks, and
/// cut the blocks where their sides are not reliably oriented.
#[must_use]
pub fn solve(problem: &PhaseProblem, min_switch_quality: i32) -> Solution {
    let sites = &problem.sites;
    let Forest { blocks, mut haplotype, component } = spanning_forest(sites, &build_edges(problem));
    let mut restrictions = restrictions(problem, &component, &haplotype);
    let by_site = restrictions_by_site(sites, &restrictions);

    for block in &blocks {
        refine(block, &mut restrictions, &by_site, &mut haplotype);
    }

    let mut pieces: Vec<PhaseBlock> = Vec::new();
    for block in &blocks {
        let switches = switch_qualities(block, &restrictions, &by_site, &haplotype);
        let weak_cuts = switches.iter().enumerate().filter(|&(_, &q)| q < min_switch_quality);
        let ends = weak_cuts.map(|(cut, _)| cut + 1).chain([block.len()]);
        let mut start = 0;
        for end in ends {
            if let (Some(piece), Some(inner)) =
                (block.get(start..end), switches.get(start..end.saturating_sub(1)))
                && piece.len() > 1
            {
                pieces.push(phased_piece(piece, inner, &haplotype));
            }
            start = end;
        }
    }

    attach_guarded(problem, &mut pieces, &mut haplotype);

    let stats = PhaseStats {
        sites: sites.len(),
        phased_sites: pieces.iter().map(Vec::len).sum(),
        blocks: pieces.len(),
        largest_block: pieces.iter().map(Vec::len).max().unwrap_or(0),
        mec_cost: restrictions.iter().map(|r| mec(r.agree, r.total).unsigned_abs()).sum(),
    };

    Solution { blocks: pieces, stats }
}

/// Least net vote a guarded site needs from one piece to join it: two agreeing
/// fragments (a link capped at 36), the bar an unguarded link needs to survive
/// the cut. One fragment was enough for a sixth of the switch errors left.
const MIN_ATTACH: i32 = 30;

/// Add each guarded site to the one piece its fragments tie it to most
/// strongly, if that is [`MIN_ATTACH`] or more, never joining two pieces.
///
/// The fragments over the guarded site and one phased site form a [`Link`],
/// capped like a solved one, and the piece's score is the sum of its sites'
/// links, signed for the guarded site's first allele being on haplotype 1.
/// Only sites phased before this pass vote, so the result does not depend on
/// the order of attaching.
fn attach_guarded(
    problem: &PhaseProblem,
    pieces: &mut [PhaseBlock],
    haplotype: &mut SiteMap<SiteAllele>,
) {
    let mut piece_of = SiteMap::filled(&problem.sites, None);
    for (index, piece) in pieces.iter().enumerate() {
        for phased in piece {
            piece_of[phased.site] = Some(index);
        }
    }

    let mut links: FxHashMap<(SiteIndex, SiteIndex), Link> = FxHashMap::default();
    for fragment in &problem.fragments {
        for guarded in &fragment.supports {
            if problem.sites[guarded.site].guard.is_none() {
                continue;
            }
            for other in fragment.supports.iter().filter(|other| other.site != guarded.site) {
                links.entry((guarded.site, other.site)).or_default().add(guarded, other);
            }
        }
    }

    let mut scores: SiteMap<SmallVec<(usize, i32), 2>> =
        SiteMap::filled(&problem.sites, SmallVec::new());
    for (&(guarded, other), link) in &links {
        let Some(piece) = piece_of[other] else { continue };
        let weight = link.capped_weight();
        let vote = if haplotype[other] == SiteAllele::First { weight } else { -weight };
        let score = &mut scores[guarded];
        match score.iter_mut().find(|(candidate, _)| *candidate == piece) {
            Some((_, total)) => *total += vote,
            None => score.push((piece, vote)),
        }
    }

    let mut attached = SiteMap::filled(&problem.sites, false);
    for site in problem.sites.indices() {
        // The lower piece wins a tie, so the hasher's order does not decide it.
        let best = scores[site]
            .iter()
            .copied()
            .max_by_key(|&(piece, score)| (score.unsigned_abs(), Reverse(piece)));
        let Some((piece, score)) = best.filter(|(_, score)| score.abs() >= MIN_ATTACH) else {
            continue;
        };
        let Some(block) = pieces.get_mut(piece) else { continue };
        let allele = if score > 0 { SiteAllele::First } else { SiteAllele::Second };
        haplotype[site] = allele;
        let at = block.partition_point(|phased| phased.site < site);
        block.insert(at, PhasedSite { site, allele, quality: Cost::calibrated(score.abs()) });
        piece_of[site] = Some(piece);
        attached[site] = true;
    }

    veto_attached(&links, &piece_of, &attached, haplotype, pieces);
}

/// Fragments over an attached site and a site of its piece that, all going
/// against the phase written for the pair, unphase the attached site.
const MIN_VETO_FRAGMENTS: u32 = 3;

/// Take an attached site back out of its piece when the fragments reaching it
/// and one of the piece's sites go against the written phase, at least
/// [`MIN_VETO_FRAGMENTS`] of them and twice as many as agree.
///
/// Attaching sums links over the whole piece, so distant sites can outvote
/// the direct reads; at that bar the direct reads were right at a seventh of
/// the switch errors and against no correct pair sampled. Only attached
/// sites: cutting a solved block on direct reads would re-litigate the
/// solve's minimum error correction.
fn veto_attached(
    links: &FxHashMap<(SiteIndex, SiteIndex), Link>,
    piece_of: &SiteMap<Option<usize>>,
    attached: &SiteMap<bool>,
    haplotype: &SiteMap<SiteAllele>,
    pieces: &mut [PhaseBlock],
) {
    let mut vetoed = SiteMap::filled(attached, false);
    for (&(guarded, other), link) in links {
        if !attached[guarded] || piece_of[other] != piece_of[guarded] {
            continue;
        }
        let (with, against) = if haplotype[guarded] == haplotype[other] {
            (link.agreeing, link.disagreeing)
        } else {
            (link.disagreeing, link.agreeing)
        };
        if against >= MIN_VETO_FRAGMENTS && against >= 2 * with {
            vetoed[guarded] = true;
        }
    }
    for piece in pieces {
        piece.retain(|phased| !vetoed[phased.site]);
    }
}

/// A piece of a block as phased sites, each rated by the weaker of the cuts
/// beside it; `switches` are the piece's own cuts, one fewer than its sites.
fn phased_piece(
    piece: &[SiteIndex],
    switches: &[i32],
    haplotype: &SiteMap<SiteAllele>,
) -> PhaseBlock {
    piece
        .iter()
        .enumerate()
        .map(|(k, &site)| {
            let before = k.checked_sub(1).and_then(|cut| switches.get(cut));
            let after = switches.get(k);
            let weakest = before.into_iter().chain(after).copied().min().unwrap_or(0);
            PhasedSite { site, allele: haplotype[site], quality: Cost::calibrated(weakest) }
        })
        .collect()
}

/// What switching the haplotypes of every site after each cut of a block
/// would cost, one per cut: the MEC rise of the restrictions crossing it.
///
/// A restriction's left side only changes at its own supports, so it adds one
/// value per gap between them, spread over that gap's cuts by a difference
/// array; the work is linear in the supports, however far apart they sit.
fn switch_qualities(
    block: &[SiteIndex],
    restrictions: &[Restriction],
    by_site: &SiteMap<Vec<usize>>,
    haplotype: &SiteMap<SiteAllele>,
) -> Vec<i32> {
    let mut steps = vec![0i32; block.len()];
    for &site in block {
        for restriction in by_site[site].iter().filter_map(|&index| restrictions.get(index)) {
            // Supports ascend by site, so this visits each restriction once.
            if restriction.supports.first().is_none_or(|first| first.site != site) {
                continue;
            }
            let agrees = |support: &AlleleSupport| haplotype[support.site] == support.allele;
            let now = mec(restriction.agree, restriction.total);
            let (mut left_agree, mut left_total) = (0, 0);
            let mut supports = restriction.supports.iter().peekable();
            while let Some(support) = supports.next() {
                let cost = i32::from(support.cost);
                left_total += cost;
                if agrees(support) {
                    left_agree += cost;
                }
                let Some(next) = supports.peek() else { break };
                let (Ok(from), Ok(to)) =
                    (block.binary_search(&support.site), block.binary_search(&next.site))
                else {
                    continue;
                };
                let right_agree = restriction.agree - left_agree;
                let right_total = restriction.total - left_total;
                let switched = left_agree + (right_total - right_agree);
                let rise = mec(switched, restriction.total) - now;
                if let Some(step) = steps.get_mut(from) {
                    *step += rise;
                }
                if let Some(step) = steps.get_mut(to) {
                    *step -= rise;
                }
            }
        }
    }
    steps
        .iter()
        .take(block.len().saturating_sub(1))
        .scan(0, |running, step| {
            *running += step;
            Some(*running)
        })
        .collect()
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

/// A fragment's supports the solve may use: a guarded site only attaches to a
/// block afterwards ([`attach_guarded`]), so it must not link two.
fn solved_supports(
    problem: &PhaseProblem,
    supports: &[AlleleSupport],
) -> SmallVec<AlleleSupport, 4> {
    supports.iter().filter(|support| problem.sites[support.site].guard.is_none()).copied().collect()
}

/// Every site pair a fragment spans, with the votes summed, strongest first.
///
/// A fragment carrying the same allele slot at both sites votes for them
/// sharing an orientation, and the vote is worth the weaker of its two
/// observations — a link is only as good as its worse end. The sum is then
/// capped by how many fragments made it ([`Link::capped_weight`]).
///
/// The order is [`spanning_forest`]'s: strongest edge first, ties broken by the
/// site pair so that it does not depend on the hasher.
fn build_edges(problem: &PhaseProblem) -> Vec<Edge> {
    let mut votes: FxHashMap<(SiteIndex, SiteIndex), Link> = FxHashMap::default();

    for fragment in &problem.fragments {
        let supports = solved_supports(problem, &fragment.supports);
        for (offset, first) in supports.iter().enumerate() {
            for second in supports.iter().skip(offset + 1) {
                votes.entry((first.site, second.site)).or_default().add(first, second);
            }
        }
    }

    // A pair whose cis and trans evidence cancelled exactly says nothing, and
    // joining a block on it would emit a coin toss as a confident phase.
    let mut edges: Vec<Edge> = votes
        .into_iter()
        .map(|((a, b), link)| Edge { a, b, weight: link.capped_weight() })
        .filter(|edge| edge.weight != 0)
        .collect();
    edges.sort_unstable_by_key(|edge| (Reverse(edge.weight.unsigned_abs()), edge.a, edge.b));
    edges
}

/// The fragments voting on one site pair: their summed signed votes, and how
/// many said "same" and "opposite".
#[derive(Debug, Clone, Copy, Default)]
struct Link {
    weight: i32,
    agreeing: u32,
    disagreeing: u32,
}

/// Most a link may weigh when its fragments disagree: 15 % of such links are
/// wrong.
const MIXED_LINK_CAP: i32 = 20;
/// Most a link of one fragment may weigh: one fragment's link is wrong at
/// Phred 26, the observation cap.
const ONE_FRAGMENT_CAP: i32 = Cost::CAP as i32;
/// Most a link of two agreeing fragments may weigh: measured wrong at Phred 36.
const TWO_FRAGMENT_CAP: i32 = 36;
/// Most any link may weigh: no link of three or more agreeing fragments was
/// wrong in 36k, but fragment errors share causes, so the sum stops growing.
const MANY_FRAGMENT_CAP: i32 = 45;

impl Link {
    /// One fragment's vote: the weaker of its two observations, for the same
    /// allele slot at both sites or against.
    fn add(&mut self, first: &AlleleSupport, second: &AlleleSupport) {
        let vote = i32::from(first.cost.min(second.cost));
        if first.allele == second.allele {
            self.weight += vote;
            self.agreeing += 1;
        } else {
            self.weight -= vote;
            self.disagreeing += 1;
        }
    }

    /// The summed vote, capped by what fragment counts like these were
    /// measured to be worth: base qualities overstate how independent a
    /// pair's fragments are.
    fn capped_weight(self) -> i32 {
        let cap = if self.agreeing > 0 && self.disagreeing > 0 {
            MIXED_LINK_CAP
        } else {
            match self.agreeing + self.disagreeing {
                0 | 1 => ONE_FRAGMENT_CAP,
                2 => TWO_FRAGMENT_CAP,
                _ => MANY_FRAGMENT_CAP,
            }
        };
        self.weight.clamp(-cap, cap)
    }
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
        for support in solved_supports(problem, &fragment.supports) {
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
        Cost, FragmentId, GtAllele, SiteKind,
        fixtures::ALT_1,
        problem::{Fragment, Guard, PhaseSite},
    };
    use SiteAllele::{First, Second};
    use seqair_types::Pos0;

    /// Sites are interchangeable to the solver — only their count and their
    /// positions matter — so the fixtures name them by index.
    fn sites(count: u32) -> Vec<PhaseSite> {
        (0..count)
            .map(|index| PhaseSite {
                pileup: index as usize,
                pos: Pos0::new(1000 + index * 10).expect("a small position"),
                kind: SiteKind::Snv,
                alleles: [GtAllele::Ref, GtAllele::Alt(ALT_1)],
                guard: None,
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
        solve(
            &PhaseProblem { sites: sites(count).into(), fragments },
            DEFAULT_MIN_SWITCH_QUALITY.into(),
        )
    }

    /// The solver without the cut, for tests of what comes before it.
    fn solve_uncut(count: u32, fragments: Vec<Fragment>) -> Solution {
        solve(&PhaseProblem { sites: sites(count).into(), fragments }, i32::MIN)
    }

    /// Each block oriented so its first site is `First`: the solver's answer
    /// is only defined up to a flip of the whole block.
    fn assignment(solution: &Solution) -> Vec<Vec<(usize, SiteAllele)>> {
        solution
            .blocks
            .iter()
            .map(|block| {
                let flip = block.first().is_some_and(|phased| phased.allele == Second);
                block
                    .iter()
                    .map(|&PhasedSite { site, allele, .. }| {
                        (site.get(), if flip { allele.flipped() } else { allele })
                    })
                    .collect()
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
    /// flip pass is what notices. Here the strongest edge (0–1, three
    /// fragments against two for the others) is the wrong one, outweighed by
    /// the two weaker edges that agree with each other.
    #[test]
    fn refinement_corrects_a_spanning_tree_mistake() {
        let solution = solve_uncut(
            3,
            vec![
                fragment(1, &[(0, First, 15), (1, Second, 15)]),
                fragment(2, &[(0, First, 15), (1, Second, 15)]),
                fragment(3, &[(0, First, 15), (1, Second, 15)]),
                fragment(4, &[(1, First, 30), (2, First, 30)]),
                fragment(5, &[(0, First, 30), (2, First, 30)]),
                fragment(6, &[(0, First, 30), (2, First, 30)]),
                fragment(7, &[(1, First, 30), (2, First, 30)]),
            ],
        );

        assert_eq!(
            assignment(&solution),
            [[(0, First), (1, First), (2, First)]],
            "four agreeing fragments beat the strongest link"
        );
    }

    #[test]
    fn a_link_is_capped_by_its_fragment_count() {
        let weights = |fragments: Vec<Fragment>| -> Vec<i32> {
            let problem = PhaseProblem { sites: sites(2).into(), fragments };
            build_edges(&problem).iter().map(|edge| edge.weight).collect()
        };
        let cis = |id| fragment(id, &[(0, First, 26), (1, First, 26)]);

        assert_eq!(weights(vec![cis(1)]), [26]);
        assert_eq!(weights(vec![cis(1), cis(2)]), [36], "not 52");
        assert_eq!(weights(vec![cis(1), cis(2), cis(3)]), [45], "not 78");
        assert_eq!(
            weights(vec![cis(1), fragment(2, &[(0, First, 5), (1, Second, 5)])]),
            [20],
            "a disagreeing fragment caps the link at 20, not 26 - 5"
        );
        assert_eq!(
            weights(vec![
                fragment(1, &[(0, First, 5), (1, First, 5)]),
                fragment(2, &[(0, First, 26), (1, Second, 26)]),
                fragment(3, &[(0, First, 26), (1, Second, 26)]),
            ]),
            [-20]
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

        let solution = solve_uncut(3 + LEAVES, fragments);

        let expected: Vec<_> = (0..3 + LEAVES as usize).map(|site| (site, First)).collect();
        assert_eq!(assignment(&solution), [expected]);
        assert_eq!(solution.stats.mec_cost, 59 * LEAVES, "only the trans fragments are wrong");
    }

    /// Two pairs of sites, each linked by three fragments, and joined by
    /// `bridges` fragments over sites 1 and 2.
    fn two_pairs_bridged(bridges: u64) -> Vec<Fragment> {
        let mut fragments: Vec<_> = (1..=3)
            .flat_map(|id| {
                [
                    fragment(id, &[(0, First, 26), (1, First, 26)]),
                    fragment(10 + id, &[(2, First, 26), (3, Second, 26)]),
                ]
            })
            .collect();
        fragments
            .extend((1..=bridges).map(|id| fragment(20 + id, &[(1, First, 26), (2, First, 26)])));
        fragments
    }

    #[test]
    fn a_single_fragment_bridge_is_cut() {
        let solution = solve_sites(4, two_pairs_bridged(1));
        assert_eq!(
            assignment(&solution),
            [[(0, First), (1, First)], [(2, First), (3, Second)]],
            "one fragment is a switch quality of 26, below 30"
        );
        assert_eq!(solution.stats.blocks, 2);

        let solution = solve_sites(4, two_pairs_bridged(2));
        assert_eq!(
            assignment(&solution),
            [[(0, First), (1, First), (2, First), (3, Second)]],
            "two agreeing fragments are 52"
        );

        let uncut =
            solve(&PhaseProblem { sites: sites(4).into(), fragments: two_pairs_bridged(1) }, 20);
        assert_eq!(uncut.stats.blocks, 1, "at Q20 one fragment joins the pairs");
    }

    #[test]
    fn a_switch_costs_what_the_crossing_fragments_lose() {
        let qualities = |fragments: Vec<Fragment>| {
            let problem = PhaseProblem { sites: sites(2).into(), fragments };
            let Forest { blocks, haplotype, component } =
                spanning_forest(&problem.sites, &build_edges(&problem));
            let restrictions = restrictions(&problem, &component, &haplotype);
            let by_site = restrictions_by_site(&problem.sites, &restrictions);
            blocks
                .iter()
                .map(|block| switch_qualities(block, &restrictions, &by_site, &haplotype))
                .collect::<Vec<_>>()
        };
        let cis = |id| fragment(id, &[(0, First, 26), (1, First, 26)]);
        let trans = |id| fragment(id, &[(0, First, 26), (1, Second, 26)]);

        assert_eq!(qualities(vec![cis(1), cis(2)]), [[52]]);
        assert_eq!(
            qualities(vec![cis(1), cis(2), trans(3), trans(4), cis(5)]),
            [[26]],
            "the two trans fragments cancel two of the cis ones"
        );
        // A cis and a trans fragment cancel to no edge at all, so a weak third
        // one is needed to make a block.
        let tie_breaker = fragment(3, &[(0, First, 1), (1, First, 1)]);
        assert_eq!(qualities(vec![cis(1), trans(2), tie_breaker]), [[1]]);
    }

    #[test]
    fn a_site_quality_is_the_calibrated_weaker_of_its_cuts() {
        // Cuts of 52 (two fragments) and 260 (ten) between three sites.
        let mut fragments: Vec<_> =
            (1..=2).map(|id| fragment(id, &[(0, First, 26), (1, First, 26)])).collect();
        fragments.extend((3..=12).map(|id| fragment(id, &[(1, First, 26), (2, First, 26)])));
        let solution = solve_sites(3, fragments);

        let qualities: Vec<_> =
            solution.blocks.iter().flatten().map(|phased| phased.quality).collect();
        assert_eq!(qualities, [Cost(31), Cost(31), Cost(37)]);

        assert_eq!(
            [9, 10, 19, 20, 29, 30, 119, 120].map(Cost::calibrated),
            [3, 8, 8, 22, 22, 31, 31, 37].map(Cost)
        );
    }

    /// Two blocks, `[0, 1]` and `[3, 4]`, each of four cis fragments, and a
    /// guarded site 2 between them: its second allele travels with site 1's
    /// first `to_left` times, its first with site 3's first `to_right` times.
    fn guarded_between(to_left: &[u8], to_right: &[u8]) -> Solution {
        let mut fragments: Vec<_> = (1..=4)
            .flat_map(|id| {
                [
                    fragment(id, &[(0, First, 26), (1, First, 26)]),
                    fragment(10 + id, &[(3, First, 26), (4, First, 26)]),
                ]
            })
            .collect();
        fragments.extend(
            (20..)
                .zip(to_left)
                .map(|(id, &cost)| fragment(id, &[(1, First, 26), (2, Second, cost)])),
        );
        fragments.extend(
            (40..)
                .zip(to_right)
                .map(|(id, &cost)| fragment(id, &[(2, First, cost), (3, First, 26)])),
        );
        let mut sites = sites(5);
        if let Some(site) = sites.get_mut(2) {
            site.guard = Some(Guard::HetDensity);
        }
        solve(&PhaseProblem { sites: sites.into(), fragments }, DEFAULT_MIN_SWITCH_QUALITY.into())
    }

    #[test]
    fn a_guarded_site_attaches_to_the_block_it_is_tied_to() {
        let solution = guarded_between(&[26, 26], &[6, 6]);
        assert_eq!(
            assignment(&solution),
            [vec![(0, First), (1, First), (2, Second)], vec![(3, First), (4, First)]],
            "two fragments are 36 to the left, 12 to the right"
        );
        let attached = solution.blocks.iter().flatten().find(|phased| phased.site.get() == 2);
        assert_eq!(attached.map(|phased| phased.quality), Some(Cost(31)));

        let solution = guarded_between(&[9], &[9]);
        assert_eq!(
            assignment(&solution),
            [[(0, First), (1, First)], [(3, First), (4, First)]],
            "9 is too weak to attach to either"
        );
    }

    #[test]
    fn a_guarded_site_attaches_on_two_fragments_not_one() {
        let solution = guarded_between(&[26], &[]);
        assert_eq!(
            assignment(&solution),
            [[(0, First), (1, First)], [(3, First), (4, First)]],
            "one fragment is 26, below the two-fragment bar"
        );

        let solution = guarded_between(&[26, 26], &[]);
        assert_eq!(
            assignment(&solution),
            [vec![(0, First), (1, First), (2, Second)], vec![(3, First), (4, First)]],
        );
    }

    /// A block `[0, 1, 2, 3]` of cis fragments and a site 4 tied cis to 0, 1
    /// and 2 by two fragments each, and to 3 by `against` trans and `with` cis
    /// fragments: the transitive votes put 4 cis with 3 whatever its direct
    /// reads say.
    fn contradicted_by_direct_reads(against: u64, with: u64, guarded: bool) -> Solution {
        let mut fragments: Vec<_> = (0..3)
            .flat_map(|site| {
                (1..=4).map(move |id| {
                    fragment(10 * u64::from(site) + id, &[(site, First, 26), (site + 1, First, 26)])
                })
            })
            .collect();
        fragments.extend((0..3).flat_map(|site| {
            (1..=2).map(move |id| {
                fragment(100 + 10 * u64::from(site) + id, &[(site, First, 26), (4, First, 26)])
            })
        }));
        fragments
            .extend((1..=against).map(|id| fragment(200 + id, &[(3, First, 26), (4, Second, 26)])));
        fragments
            .extend((1..=with).map(|id| fragment(300 + id, &[(3, First, 26), (4, First, 26)])));
        let mut sites = sites(5);
        if guarded && let Some(site) = sites.get_mut(4) {
            site.guard = Some(Guard::HetDensity);
        }
        solve(&PhaseProblem { sites: sites.into(), fragments }, DEFAULT_MIN_SWITCH_QUALITY.into())
    }

    #[test]
    fn direct_reads_against_an_attachment_unphase_the_site() {
        let block = [(0, First), (1, First), (2, First), (3, First)];
        let with_four = [(0, First), (1, First), (2, First), (3, First), (4, First)];

        assert_eq!(assignment(&contradicted_by_direct_reads(3, 0, true)), [block], "3:0 vetoes");
        assert_eq!(
            assignment(&contradicted_by_direct_reads(3, 2, true)),
            [with_four],
            "3:2 is not two to one"
        );
        assert_eq!(
            assignment(&contradicted_by_direct_reads(2, 0, true)),
            [with_four],
            "two fragments are too few"
        );
        assert_eq!(
            assignment(&contradicted_by_direct_reads(3, 0, false)),
            [with_four],
            "a solved site is left to the solve"
        );
    }

    /// Unguarded, two fragments of 15 would link site 2 to the right block at
    /// 30, a switch quality the cut keeps, and so join the blocks.
    #[test]
    fn a_guarded_site_tied_to_two_blocks_joins_only_the_stronger() {
        let solution = guarded_between(&[26, 26, 26, 26], &[15, 15]);
        assert_eq!(
            assignment(&solution),
            [vec![(0, First), (1, First), (2, Second)], vec![(3, First), (4, First)]],
        );
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
                    phased.iter().map(move |p| (p.site.get(), (block, p.allele)))
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
                let solution = solve_uncut(count, fragments.clone());
                let cost = mec_from_scratch(&fragments, &solution.blocks);
                prop_assert_eq!(solution.stats.mec_cost, cost);

                for (block, phased) in solution.blocks.iter().enumerate() {
                    for flip in 0..phased.len() {
                        let mut flipped = solution.blocks.clone();
                        if let Some(phased) = flipped.get_mut(block).and_then(|b| b.get_mut(flip)) {
                            phased.allele = phased.allele.flipped();
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
                    .is_some_and(|phased| truth[phased.site.get()] != phased.allele);
                for &PhasedSite { site, allele, .. } in block {
                    let expected =
                        if flipped { truth[site.get()].flipped() } else { truth[site.get()] };
                    prop_assert_eq!(allele, expected, "site {} is on the wrong haplotype", site.get());
                }
            }
        }
    }
}
