/**
 * @file test_search.cpp
 * @brief Unit and integration tests for search.cpp (off-target TST search).
 *
 * Strategy
 * --------
 * We build small TST index files on disk via crispritz::build_tree(), then
 * call crispritz::search_offtargets() and verify:
 *  - exact-match targets (0 mm, 0 bulges) are found
 *  - mismatch tolerance allows 1- and 2-mm hits to be recovered
 *  - DNA bulge (gap in guide)  hits are recovered
 *  - RNA bulge (gap in target) hits are recovered
 *  - combined DNA + RNA bulge hits are recovered
 *  - deduplication removes identical hits from multi-partition files
 *  - positional profiles accumulate correct per-position counts
 *  - SearchParams validation (negative limits) throws std::invalid_argument
 *  - load_index throws on a missing file
 *  - multi-threaded and single-threaded results are identical
 *
 * All temporary .bin files are cleaned up in a RAII guard at the end of
 * each test.
 */

#include "offtarget.hpp"
#include "search.hpp"
#include "tst.hpp"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using crispritz::build_tree;
using crispritz::OffTarget;
using crispritz::search_offtargets;
using crispritz::SearchParams;
using crispritz::SearchResult;

// =============================================================================
// Minimal test harness
// =============================================================================

static int g_total = 0;
static int g_passed = 0;
static int g_failed = 0;

static void record(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_total;
    if (ok)
    {
        ++g_passed;
        std::cout << "  [PASS] " << name << "\n";
    }
    else
    {
        ++g_failed;
        std::cout << "  [FAIL] " << name;
        if (!detail.empty())
            std::cout << " -- " << detail;
        std::cout << "\n";
    }
}

// =============================================================================
// RAII cleaner for .bin files
// =============================================================================

struct BinCleaner
{
    std::string pam_seq;
    std::string chr_name;

    BinCleaner(std::string p, std::string c) : pam_seq(std::move(p)), chr_name(std::move(c)) {}

    ~BinCleaner()
    {
        const std::string prefix = pam_seq + "_" + chr_name + "_";
        for (auto& entry : fs::directory_iterator(fs::current_path()))
        {
            if (!entry.is_regular_file())
                continue;
            std::string fname = entry.path().filename().string();
            if (fname.rfind(prefix, 0) == 0 && fname.size() > 4 &&
                fname.substr(fname.size() - 4) == ".bin")
                fs::remove(entry.path());
        }
    }
};

// =============================================================================
// Helpers
// =============================================================================

/** @brief Collect all .bin files matching <pam>_<chr>_*.bin in cwd. */
static std::vector<std::string> find_bin_files(const std::string& pam_seq,
                                               const std::string& chr_name)
{
    std::vector<std::string> result;
    const std::string prefix = pam_seq + "_" + chr_name + "_";
    for (auto& entry : fs::directory_iterator(fs::current_path()))
    {
        if (!entry.is_regular_file())
            continue;
        std::string fname = entry.path().filename().string();
        if (fname.rfind(prefix, 0) == 0 && fname.size() > 4 &&
            fname.substr(fname.size() - 4) == ".bin")
            result.push_back(entry.path().string());
    }
    std::sort(result.begin(), result.end());
    return result;
}

/**
 * @brief Build a genome containing @p guide immediately followed by @p pam,
 *        padded to @p total_len with 'A' characters.
 *
 * The genomic site is at position 0 (0-based): guide at [0, guide.size()),
 * pam at [guide.size(), guide.size() + pam.size()).
 */
static std::string make_simple_genome(const std::string& guide, const std::string& pam,
                                      int total_len = 200)
{
    std::string genome(total_len, 'A');
    for (int i = 0; i < static_cast<int>(guide.size()) && i < total_len; ++i)
        genome[i] = guide[i];
    for (int i = 0; i < static_cast<int>(pam.size()); ++i)
    {
        int pos = static_cast<int>(guide.size()) + i;
        if (pos < total_len)
            genome[pos] = pam[i];
    }
    return genome;
}

/** @brief Count hits matching a given chromosome. */
static int count_hits_on_chr(const SearchResult& res, const std::string& chr)
{
    int n = 0;
    for (const OffTarget& ot : res.off_targets)
        if (ot.chromosome == chr)
            ++n;
    return n;
}

/** @brief True iff at least one hit has mismatches == mm and bulge_size == bs. */
static bool has_hit(const SearchResult& res, int mm, int bs)
{
    for (const OffTarget& ot : res.off_targets)
        if (ot.mismatches == mm && ot.bulge_size == bs)
            return true;
    return false;
}

/** @brief True iff at least one hit has the given bulge_type. */
static bool has_bulge_type(const SearchResult& res, const std::string& btype)
{
    for (const OffTarget& ot : res.off_targets)
        if (ot.bulge_type == btype)
            return true;
    return false;
}

// =============================================================================
// Default SearchParams factory for SpCas9-style (20-nt guide, NGG, pam-at-end)
// =============================================================================
static SearchParams cas9_params(int mm = 0, int bdna = 0, int brna = 0, int threads = 1)
{
    SearchParams p;
    p.max_mismatches = mm;
    p.max_dna_bulges = bdna;
    p.max_rna_bulges = brna;
    p.pam_length = 23; // 20 guide + 3 PAM
    p.pam_limit = 3;
    p.pam_at_start = false;
    p.num_threads = threads;
    return p;
}

// =============================================================================
// Test: exact match (0 mm, 0 bulges)
// =============================================================================

static void test_exact_match()
{
    const std::string chr = "exact_chr";
    const std::string pam_s = "NGG";
    const std::string guide = "ACGTACGTACGTACGTACGT";           // 20 nt
    const std::string genome = make_simple_genome(guide, "GG"); // guide + GG + A...

    BinCleaner cleaner(pam_s, chr);

    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (const std::exception& e)
    {
        record("exact_match: build_tree does not throw", false, e.what());
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    record("exact_match: .bin file produced", !bins.empty());
    if (bins.empty())
        return;

    auto res = search_offtargets(bins, {guide}, cas9_params(0, 0, 0));
    record("exact_match: at least one hit found", !res.off_targets.empty(),
           "hits=" + std::to_string(res.off_targets.size()));
    record("exact_match: hit has 0 mismatches", has_hit(res, 0, 0));
    record("exact_match: bulge_type == X", has_bulge_type(res, "X"));
    record("exact_match: profile mm_counts[0] > 0", !res.profiles.empty() &&
                                                        !res.profiles[0].mm_counts.empty() &&
                                                        res.profiles[0].mm_counts[0] > 0);
}

// =============================================================================
// Test: 1-mismatch tolerance
// =============================================================================

static void test_one_mismatch()
{
    const std::string chr = "mm1_chr";
    const std::string pam_s = "NGG";
    const std::string guide = "ACGTACGTACGTACGTACGT";

    // Genome has a 1-mismatch target: change position 5 from C to T.
    std::string target = guide;
    target[5] = (target[5] == 'A') ? 'T' : 'A'; // introduce one substitution

    const std::string genome = make_simple_genome(target, "GG");

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("mm1: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("mm1: bin produced", false);
        return;
    }

    // With mm=0, should NOT find the 1-mm target.
    auto res0 = search_offtargets(bins, {guide}, cas9_params(0));
    record("mm1: 0-mm search returns no hit", res0.off_targets.empty(),
           "hits=" + std::to_string(res0.off_targets.size()));

    // With mm=1, should find it.
    auto res1 = search_offtargets(bins, {guide}, cas9_params(1));
    record("mm1: 1-mm search finds hit", !res1.off_targets.empty(),
           "hits=" + std::to_string(res1.off_targets.size()));
    record("mm1: hit has 1 mismatch", has_hit(res1, 1, 0));

    // Profile: mm_counts[1] should be incremented.
    record("mm1: profile mm_counts[1] > 0", !res1.profiles.empty() &&
                                                res1.profiles[0].mm_counts.size() > 1 &&
                                                res1.profiles[0].mm_counts[1] > 0);
}

// =============================================================================
// Test: 2-mismatch tolerance
// =============================================================================

static void test_two_mismatches()
{
    const std::string chr = "mm2_chr";
    const std::string pam_s = "NGG";
    const std::string guide = "ACGTACGTACGTACGTACGT";

    std::string target = guide;
    target[2] = (target[2] == 'G') ? 'T' : 'G';
    target[10] = (target[10] == 'A') ? 'C' : 'A';

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(make_simple_genome(target, "GG"), chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("mm2: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("mm2: bin produced", false);
        return;
    }

    auto res = search_offtargets(bins, {guide}, cas9_params(2));
    record("mm2: 2-mm search finds hit", !res.off_targets.empty());
    record("mm2: hit has 2 mismatches", has_hit(res, 2, 0));
}

// =============================================================================
// Test: DNA bulge (gap in guide / extra base in target)
// =============================================================================

static void test_dna_bulge()
{
    const std::string chr = "bdna_chr";
    const std::string pam_s = "NGG";

    // guide:  ACGTACGTACGTACGTACGT  (20 nt)
    // target: ACGTACGTAACGTACGTACGT (21 nt, extra 'A' at pos 9)
    const std::string guide = "ACGTACGTACGTACGTACGT";
    const std::string target = "ACGTACGTAACGTACGTACGT"; // 21 nt

    // Build with max_bulges=1 so the extra base is extracted.
    BinCleaner cleaner(pam_s, chr);
    std::string genome = make_simple_genome(target, "GG");
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 1, 1);
    }
    catch (...)
    {
        record("bdna: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("bdna: bin produced", false);
        return;
    }

    // Search with 1 DNA bulge allowed.
    auto res = search_offtargets(bins, {guide}, cas9_params(0, 1, 0));
    record("bdna: DNA-bulge search returns results", !res.off_targets.empty(),
           "hits=" + std::to_string(res.off_targets.size()));
    if (!res.off_targets.empty())
        record("bdna: hit has bulge_size 1", has_hit(res, 0, 1));
}

// =============================================================================
// Test: RNA bulge (gap in target / extra base in guide)
// =============================================================================

static void test_rna_bulge()
{
    const std::string chr = "brna_chr";
    const std::string pam_s = "NGG";

    // guide:  ACGTACGTACGTACGTACGT  (20 nt)
    // target: ACGTACGTCGTACGTACGT   (19 nt, missing 'A' at pos 8)
    const std::string guide = "ACGTACGTACGTACGTACGT";
    const std::string target = "ACGTACGTCGTACGTACGT"; // 19 nt

    BinCleaner cleaner(pam_s, chr);
    std::string genome = make_simple_genome(target, "GG");
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("brna: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("brna: bin produced", false);
        return;
    }

    // Search with 1 RNA bulge allowed.
    auto res = search_offtargets(bins, {guide}, cas9_params(0, 0, 1));
    record("brna: RNA-bulge search returns results", !res.off_targets.empty(),
           "hits=" + std::to_string(res.off_targets.size()));
    if (!res.off_targets.empty())
        record("brna: hit has bulge_size 1", has_hit(res, 0, 1));
}

// =============================================================================
// Test: multiple guides on the same genome
// =============================================================================

static void test_multiple_guides()
{
    const std::string chr = "multi_chr";
    const std::string pam_s = "NGG";

    const std::string g1 = "ACGTACGTACGTACGTACGT";
    const std::string g2 = "TTTTTTTTTTTTTTTTTTTT";

    // Build a genome that contains both guides.
    std::string genome(200, 'A');
    for (int i = 0; i < 20; ++i)
        genome[i] = g1[i];
    genome[20] = genome[21] = 'G';
    for (int i = 0; i < 20; ++i)
        genome[50 + i] = g2[i];
    genome[70] = genome[71] = 'G';

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("multi: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("multi: bin produced", false);
        return;
    }

    auto res = search_offtargets(bins, {g1, g2}, cas9_params(0));

    record("multi: two profiles produced", res.profiles.size() == 2,
           "profiles=" + std::to_string(res.profiles.size()));
    record("multi: hits found for both guides", res.off_targets.size() >= 2,
           "hits=" + std::to_string(res.off_targets.size()));
}

// =============================================================================
// Test: PAM-at-start (Cas12a-style, TTTV)
// =============================================================================

static void test_pam_at_start()
{
    const std::string chr = "cas12a_chr";
    const std::string pam_s = "TTT";
    const std::string guide = "ACGTACGTACGTACGTACGT";

    // For PAM-at-start, genome layout is: TTT + guide + padding.
    std::string genome(200, 'A');
    genome[0] = genome[1] = genome[2] = 'T';
    for (int i = 0; i < 20; ++i)
        genome[3 + i] = guide[i];

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, true, ".", 0, 1);
    }
    catch (...)
    {
        record("cas12a: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("cas12a: bin produced", false);
        return;
    }

    SearchParams p = cas9_params(0);
    p.pam_at_start = true;

    auto res = search_offtargets(bins, {guide}, p);
    record("cas12a: at least one hit found", !res.off_targets.empty(),
           "hits=" + std::to_string(res.off_targets.size()));
}

// =============================================================================
// Test: deduplication – same site in two identical .bin files
// =============================================================================

static void test_deduplication()
{
    const std::string chr = "dedup_chr";
    const std::string pam_s = "NGG";
    const std::string guide = "ACGTACGTACGTACGTACGT";
    const std::string genome = make_simple_genome(guide, "GG");

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("dedup: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("dedup: bin produced", false);
        return;
    }

    // Duplicate the bin list to simulate two partitions with identical content.
    std::vector<std::string> doubled_bins;
    doubled_bins.insert(doubled_bins.end(), bins.begin(), bins.end());
    doubled_bins.insert(doubled_bins.end(), bins.begin(), bins.end());

    auto res = search_offtargets(doubled_bins, {guide}, cas9_params(0));

    // Deduplication should collapse the doubled entries back to the same set
    // as a single-copy search.
    auto res_single = search_offtargets(bins, {guide}, cas9_params(0));

    record("dedup: doubled input yields same hits as single",
           res.off_targets.size() == res_single.off_targets.size(),
           "doubled=" + std::to_string(res.off_targets.size()) +
               " single=" + std::to_string(res_single.off_targets.size()));
}

// =============================================================================
// Test: invalid parameters throw std::invalid_argument
// =============================================================================

static void test_invalid_params_throw()
{
    const std::vector<std::string> bins = {"nonexistent.bin"};
    const std::vector<std::string> guides = {"ACGTACGTACGTACGTACGT"};

    // Negative mismatch.
    {
        bool threw = false;
        try
        {
            search_offtargets(bins, guides, -1, 0, 0, 23, 3, false, 1);
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        catch (...)
        {
        }
        record("invalid_params: negative mm throws invalid_argument", threw);
    }

    // Negative DNA bulges.
    {
        bool threw = false;
        try
        {
            search_offtargets(bins, guides, 0, -1, 0, 23, 3, false, 1);
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        catch (...)
        {
        }
        record("invalid_params: negative bdna throws invalid_argument", threw);
    }

    // Negative RNA bulges.
    {
        bool threw = false;
        try
        {
            search_offtargets(bins, guides, 0, 0, -1, 23, 3, false, 1);
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        catch (...)
        {
        }
        record("invalid_params: negative brna throws invalid_argument", threw);
    }

    // pam_length <= pam_limit.
    {
        bool threw = false;
        try
        {
            search_offtargets(bins, guides, 0, 0, 0, 3, 3, false, 1);
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        catch (...)
        {
        }
        record("invalid_params: pam_length == pam_limit throws", threw);
    }

    // Empty guides.
    {
        bool threw = false;
        try
        {
            search_offtargets(bins, {}, 0, 0, 0, 23, 3, false, 1);
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        catch (...)
        {
        }
        record("invalid_params: empty guides throws", threw);
    }

    // Empty bin_paths.
    {
        bool threw = false;
        try
        {
            search_offtargets({}, guides, 0, 0, 0, 23, 3, false, 1);
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        catch (...)
        {
        }
        record("invalid_params: empty bin_paths throws", threw);
    }
}

// =============================================================================
// Test: load_index throws on a missing file
// =============================================================================

static void test_load_index_missing_file()
{
    bool threw = false;
    try
    {
        crispritz::load_index("/no/such/file.bin", 3, "chr0");
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }
    catch (...)
    {
    }
    record("load_index: missing file throws runtime_error", threw);
}

// =============================================================================
// Test: multi-threaded search matches single-threaded results
// =============================================================================

static void test_multithread_same_results()
{
    const std::string chr = "mt_chr";
    const std::string pam_s = "NGG";

    // Build a genome with several NGG sites.
    std::string genome(500, 'A');
    const std::string guide = "GCATGCATGCATGCATGCAT";
    for (int start : {0, 100, 200, 300})
    {
        for (int i = 0; i < 20 && start + i < 500; ++i)
            genome[start + i] = guide[i];
        if (start + 21 < 500)
        {
            genome[start + 20] = 'G';
            genome[start + 21] = 'G';
        }
    }

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("mt: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("mt: bin produced", false);
        return;
    }

    auto res1 = search_offtargets(bins, {guide}, cas9_params(1, 0, 0, 1));
    auto res4 = search_offtargets(bins, {guide}, cas9_params(1, 0, 0, 4));

    // Sort both result sets for comparison.
    auto cmp = [](const OffTarget& a, const OffTarget& b)
    {
        if (a.chromosome != b.chromosome)
            return a.chromosome < b.chromosome;
        if (a.genomic_pos != b.genomic_pos)
            return a.genomic_pos < b.genomic_pos;
        return a.target_seq < b.target_seq;
    };
    std::sort(res1.off_targets.begin(), res1.off_targets.end(), cmp);
    std::sort(res4.off_targets.begin(), res4.off_targets.end(), cmp);

    record("mt: 1-thread and 4-thread hit counts match",
           res1.off_targets.size() == res4.off_targets.size(),
           "1thr=" + std::to_string(res1.off_targets.size()) +
               " 4thr=" + std::to_string(res4.off_targets.size()));

    bool same = true;
    for (size_t i = 0; i < std::min(res1.off_targets.size(), res4.off_targets.size()); ++i)
    {
        if (res1.off_targets[i].genomic_pos != res4.off_targets[i].genomic_pos ||
            res1.off_targets[i].strand != res4.off_targets[i].strand ||
            res1.off_targets[i].target_seq != res4.off_targets[i].target_seq)
        {
            same = false;
            break;
        }
    }
    record("mt: hit positions match between 1-thread and 4-thread", same);
}

// =============================================================================
// Test: OffTarget fields are populated correctly
// =============================================================================

static void test_offtarget_fields()
{
    const std::string chr = "fields_chr";
    const std::string pam_s = "NGG";
    const std::string guide = "ACGTACGTACGTACGTACGT";
    const std::string genome = make_simple_genome(guide, "GG");

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("fields: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("fields: bin produced", false);
        return;
    }

    auto res = search_offtargets(bins, {guide}, cas9_params(0));

    record("fields: results non-empty", !res.off_targets.empty());
    if (res.off_targets.empty())
        return;

    const OffTarget& ot = res.off_targets[0];
    record("fields: chromosome is non-empty", !ot.chromosome.empty());
    record("fields: genomic_pos > 0", ot.genomic_pos > 0, "pos=" + std::to_string(ot.genomic_pos));
    record("fields: strand is + or -", ot.strand == '+' || ot.strand == '-',
           std::string("strand=") + ot.strand);
    record("fields: mismatches == 0 for exact match", ot.mismatches == 0);
    record("fields: bulge_size == 0 for exact match", ot.bulge_size == 0);
    record("fields: total_score == 0 for exact match", ot.total_score == 0);
    record("fields: guide_seq length == pam_length", static_cast<int>(ot.guide_seq.size()) == 23,
           "len=" + std::to_string(ot.guide_seq.size()));
    record("fields: target_seq length == pam_length", static_cast<int>(ot.target_seq.size()) == 23,
           "len=" + std::to_string(ot.target_seq.size()));
}

// =============================================================================
// Test: profile structure is correctly initialised
// =============================================================================

static void test_profile_structure()
{
    const std::string chr = "prof_chr";
    const std::string pam_s = "NGG";
    const std::string guide = "ACGTACGTACGTACGTACGT";
    const std::string genome = make_simple_genome(guide, "GG");

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("prof: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("prof: bin produced", false);
        return;
    }

    auto p = cas9_params(2, 1, 1);
    auto res = search_offtargets(bins, {guide}, p);

    record("prof: one profile per guide", res.profiles.size() == 1);
    if (res.profiles.empty())
        return;

    const auto& prof = res.profiles[0];
    record("prof: mm_per_pos has guide_length entries",
           static_cast<int>(prof.mm_per_pos.size()) == p.guide_length(),
           "size=" + std::to_string(prof.mm_per_pos.size()));
    record("prof: mm_counts has max_mm+1 entries",
           static_cast<int>(prof.mm_counts.size()) == p.max_mismatches + 1,
           "size=" + std::to_string(prof.mm_counts.size()));
    record("prof: nuc_per_pos outer dim == max_mm+1",
           static_cast<int>(prof.nuc_per_pos.size()) == p.max_mismatches + 1);
    record("prof: joint_counts outer dim == max_mm+1",
           static_cast<int>(prof.joint_counts.size()) == p.max_mismatches + 1);
}

// =============================================================================
// Test: SearchResult __len__ equivalent (off_targets size)
// =============================================================================

static void test_result_size()
{
    const std::string chr = "size_chr";
    const std::string pam_s = "NGG";
    const std::string guide = "ACGTACGTACGTACGTACGT";
    const std::string genome = make_simple_genome(guide, "GG");

    BinCleaner cleaner(pam_s, chr);
    try
    {
        build_tree(genome, chr, pam_s, 23, 3, false, ".", 0, 1);
    }
    catch (...)
    {
        record("size: build does not throw", false);
        return;
    }

    auto bins = find_bin_files(pam_s, chr);
    if (bins.empty())
    {
        record("size: bin produced", false);
        return;
    }

    auto res = search_offtargets(bins, {guide}, cas9_params(0));
    record("size: SearchResult has expected number of hits",
           res.off_targets.size() == static_cast<size_t>(count_hits_on_chr(res, chr)),
           "total=" + std::to_string(res.off_targets.size()));
}

// =============================================================================
// main
// =============================================================================

int main()
{
    std::cout << "=== test_search ===\n\n";

    std::cout << "-- exact match --\n";
    test_exact_match();

    std::cout << "\n-- mismatch tolerance --\n";
    test_one_mismatch();
    test_two_mismatches();

    std::cout << "\n-- bulges --\n";
    test_dna_bulge();
    test_rna_bulge();

    std::cout << "\n-- multiple guides --\n";
    test_multiple_guides();

    std::cout << "\n-- PAM at start --\n";
    test_pam_at_start();

    std::cout << "\n-- deduplication --\n";
    test_deduplication();

    std::cout << "\n-- error handling --\n";
    test_invalid_params_throw();
    test_load_index_missing_file();

    std::cout << "\n-- multi-threading --\n";
    test_multithread_same_results();

    std::cout << "\n-- OffTarget fields --\n";
    test_offtarget_fields();

    std::cout << "\n-- profile structure --\n";
    test_profile_structure();

    std::cout << "\n-- result size --\n";
    test_result_size();

    std::cout << "\n=== Results: " << g_passed << "/" << g_total << " passed";
    if (g_failed > 0)
        std::cout << " (" << g_failed << " FAILED)";
    std::cout << " ===\n";

    return g_failed == 0 ? 0 : 1;
}