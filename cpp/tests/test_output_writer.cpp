/**
 * @file test_output_writer.cpp
 * @brief Unit tests for the output generation layer (output_writer.hpp/.cpp).
 *
 * Coverage:
 *   - valid output      : header + rows for both formats; column correctness
 *   - empty output      : header-only when there are no records; header
 * suppressed
 *   - malformed path    : write_to_file to a bad path throws
 *   - multiple records  : flat vector and SearchResult (multi-guide) paths
 *   - formatter factory : both formats; bad enum value rejected
 *   - construction      : null formatter rejected; format_name reported
 *
 * Most assertions run against an in-memory std::ostringstream so they need no
 * filesystem; the malformed-path and round-trip-file cases use real files.
 *
 * Shared lightweight record()/g_* harness, consistent with the other suites.
 */

#include "output_writer.hpp"

#include "offtarget.hpp"
#include "search_configuration.hpp"
#include "tst_search.hpp"

#include <cstdio> // std::remove
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace crispritz;

// =============================================================================
// Harness
// =============================================================================

static int g_total = 0, g_passed = 0, g_failed = 0;

static void record(const std::string &name, bool ok,
                   const std::string &detail = "") {
  ++g_total;
  if (ok) {
    ++g_passed;
    std::cout << "  [PASS] " << name << '\n';
  } else {
    ++g_failed;
    std::cout << "  [FAIL] " << name;
    if (!detail.empty())
      std::cout << " -- " << detail;
    std::cout << '\n';
  }
}

// =============================================================================
// Helpers
// =============================================================================

static OffTarget hit(const std::string &chrom, int pos, Strand strand,
                     int mm = 1, int bdna = 0, int brna = 0) {
  return OffTarget{
      chrom, pos,  strand, "ACGTACGTACGTACGTACGTNGG", "ACGTACGTACGTACGTACGTaGG",
      mm,    bdna, brna};
}

// Split into lines (no trailing empty element if text ends with '\n').
static std::vector<std::string> lines(const std::string &s) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == '\n') {
      out.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  if (!cur.empty())
    out.push_back(cur);
  return out;
}

static std::vector<std::string> fields(const std::string &line,
                                       char sep = '\t') {
  std::vector<std::string> out;
  std::string cur;
  for (char c : line) {
    if (c == sep) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  out.push_back(cur);
  return out;
}

// =============================================================================
// Formatter factory & construction
// =============================================================================

static void test_factory_tsv() {
  auto f = make_formatter(OutputFormat::Tsv);
  record("make_formatter(Tsv) not null", f != nullptr);
  record("make_formatter(Tsv) name == tsv", f && f->name() == "tsv");
}

static void test_factory_targets() {
  auto f = make_formatter(OutputFormat::Targets);
  record("make_formatter(Targets) not null", f != nullptr);
  record("make_formatter(Targets) name == targets",
         f && f->name() == "targets");
}

static void test_factory_bad_enum() {
  bool threw = false;
  try {
    (void)make_formatter(static_cast<OutputFormat>(99));
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("make_formatter(bad enum) throws invalid_argument", threw);
}

static void test_writer_null_formatter_rejected() {
  bool threw = false;
  try {
    OutputWriter w{std::unique_ptr<OffTargetFormatter>(nullptr)};
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("OutputWriter(null formatter) throws invalid_argument", threw);
}

static void test_writer_format_name() {
  OutputWriter w_tsv{OutputFormat::Tsv};
  OutputWriter w_tg{OutputFormat::Targets};
  record("writer format_name tsv", w_tsv.format_name() == "tsv");
  record("writer format_name targets", w_tg.format_name() == "targets");
}

// =============================================================================
// Valid output — TSV
// =============================================================================

static void test_tsv_header_and_row() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  std::size_t n = w.write({hit("chr1", 100, Strand::Forward)}, os);

  auto ls = lines(os.str());
  record("tsv: returns 1 written", n == 1u);
  record("tsv: 2 lines (header + 1 row)", ls.size() == 2u,
         "lines=" + std::to_string(ls.size()));

  if (ls.size() == 2u) {
    auto hdr = fields(ls[0]);
    record("tsv: header has 9 columns", hdr.size() == 9u);
    record("tsv: header[0] == chrom", !hdr.empty() && hdr[0] == "chrom");
    record("tsv: header[8] == bulge_type",
           hdr.size() == 9u && hdr[8] == "bulge_type");

    auto row = fields(ls[1]);
    record("tsv: row has 9 columns", row.size() == 9u);
    record("tsv: row chrom == chr1", row.size() > 0 && row[0] == "chr1");
    record("tsv: row pos == 100", row.size() > 1 && row[1] == "100");
    record("tsv: row strand == +", row.size() > 2 && row[2] == "+");
  }
}

// =============================================================================
// Valid output — Targets
// =============================================================================

static void test_targets_header_and_row() {
  OutputWriter w{OutputFormat::Targets};
  std::ostringstream os;
  // mm=2, dna=1, rna=1 → bulge_size=2, total=4, bulge_type=DNA,RNA
  w.write({hit("chr7", 42, Strand::Reverse, 2, 1, 1)}, os);

  auto ls = lines(os.str());
  record("targets: 2 lines (header + 1 row)", ls.size() == 2u);

  if (ls.size() == 2u) {
    auto hdr = fields(ls[0]);
    record("targets: header has 9 columns", hdr.size() == 9u);
    record("targets: header[0] == bulge_type",
           !hdr.empty() && hdr[0] == "bulge_type");

    auto row = fields(ls[1]);
    record("targets: row[0] bulge_type == DNA,RNA",
           row.size() > 0 && row[0] == "DNA,RNA");
    record("targets: row[3] chrom == chr7", row.size() > 3 && row[3] == "chr7");
    record("targets: row[4] pos == 42", row.size() > 4 && row[4] == "42");
    record("targets: row[5] strand == -", row.size() > 5 && row[5] == "-");
    record("targets: row[6] mismatches == 2", row.size() > 6 && row[6] == "2");
    record("targets: row[7] bulge_size == 2", row.size() > 7 && row[7] == "2");
    record("targets: row[8] total == 4", row.size() > 8 && row[8] == "4");
  }
}

// =============================================================================
// Empty output
// =============================================================================

static void test_empty_with_header() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  std::size_t n = w.write(std::vector<OffTarget>{}, os, /*write_header=*/true);
  auto ls = lines(os.str());
  record("empty+header: 0 rows written", n == 0u);
  record("empty+header: exactly 1 line (header only)", ls.size() == 1u,
         "lines=" + std::to_string(ls.size()));
}

static void test_empty_without_header() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  std::size_t n = w.write(std::vector<OffTarget>{}, os, /*write_header=*/false);
  record("empty-no-header: 0 rows written", n == 0u);
  record("empty-no-header: output is empty", os.str().empty());
}

static void test_no_header_flag_suppresses_header() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  w.write({hit("chr1", 1, Strand::Forward)}, os, /*write_header=*/false);
  auto ls = lines(os.str());
  record("no-header flag: 1 line (row only)", ls.size() == 1u);
  // first field of the single line should be the data, not "chrom"
  record("no-header flag: line is data not header",
         !ls.empty() && fields(ls[0])[0] == "chr1");
}

// =============================================================================
// Multiple records
// =============================================================================

static void test_multiple_records_vector() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  std::vector<OffTarget> recs = {
      hit("chr1", 100, Strand::Forward),
      hit("chr1", 200, Strand::Reverse),
      hit("chr2", 50, Strand::Forward),
  };
  std::size_t n = w.write(recs, os);
  auto ls = lines(os.str());
  record("multi: 3 rows written", n == 3u);
  record("multi: 4 lines (header + 3)", ls.size() == 4u,
         "lines=" + std::to_string(ls.size()));
}

static void test_multiple_records_search_result() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;

  SearchResult result;
  result.source_path = "chrX_NGG_1.bin";
  result.hits_by_guide.resize(3);
  result.hits_by_guide[0] = {hit("chr1", 100, Strand::Forward),
                             hit("chr1", 150, Strand::Reverse)};
  result.hits_by_guide[1] = {hit("chr2", 300, Strand::Forward)};
  // guide[2] intentionally empty

  std::size_t n = w.write(result, os);
  auto ls = lines(os.str());
  record("SearchResult: 3 rows written (2 + 1 + 0)", n == 3u,
         "n=" + std::to_string(n));
  record("SearchResult: 4 lines (header + 3)", ls.size() == 4u);
}

static void test_search_result_all_empty() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  SearchResult result;
  result.hits_by_guide.resize(5); // five guides, no hits
  std::size_t n = w.write(result, os);
  auto ls = lines(os.str());
  record("SearchResult all-empty: 0 rows", n == 0u);
  record("SearchResult all-empty: header only", ls.size() == 1u);
}

// =============================================================================
// File API + malformed path
// =============================================================================

static void test_file_round_trip() {
  OutputWriter w{OutputFormat::Tsv};
  const std::string path = "/tmp/crispritz_ow_roundtrip.tsv";
  std::vector<OffTarget> recs = {
      hit("chr1", 100, Strand::Forward),
      hit("chr2", 200, Strand::Reverse),
  };
  std::size_t n = w.write_to_file(recs, path);
  record("file: 2 rows written", n == 2u);

  std::ifstream in(path);
  record("file: exists and opens", in.is_open());
  std::string content((std::istreambuf_iterator<char>(in)),
                      std::istreambuf_iterator<char>());
  in.close();
  auto ls = lines(content);
  record("file: 3 lines (header + 2)", ls.size() == 3u,
         "lines=" + std::to_string(ls.size()));
  std::remove(path.c_str());
}

static void test_malformed_path_throws() {
  OutputWriter w{OutputFormat::Tsv};
  // A path inside a non-existent directory cannot be opened for writing.
  const std::string bad = "/nonexistent_dir_xyz/sub/output.tsv";
  bool threw = false;
  try {
    (void)w.write_to_file({hit("chr1", 1, Strand::Forward)}, bad);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  record("malformed path: write_to_file throws runtime_error", threw);
}

static void test_malformed_path_searchresult_throws() {
  OutputWriter w{OutputFormat::Targets};
  const std::string bad = "/nonexistent_dir_xyz/sub/result.targets";
  SearchResult result;
  result.hits_by_guide.resize(1);
  result.hits_by_guide[0] = {hit("chr1", 1, Strand::Forward)};
  bool threw = false;
  try {
    (void)w.write_to_file(result, bad);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  record("malformed path (SearchResult): throws runtime_error", threw);
}

// =============================================================================
// Separation-of-concerns sanity: row values identical across stream and file
// =============================================================================

static void test_stream_file_consistency() {
  OutputWriter w{OutputFormat::Targets};
  std::vector<OffTarget> recs = {hit("chr3", 999, Strand::Reverse, 1, 1, 0)};

  std::ostringstream os;
  w.write(recs, os);

  const std::string path = "/tmp/crispritz_ow_consistency.targets";
  w.write_to_file(recs, path);
  std::ifstream in(path);
  std::string file_content((std::istreambuf_iterator<char>(in)),
                           std::istreambuf_iterator<char>());
  in.close();
  std::remove(path.c_str());

  record("stream and file produce identical bytes", os.str() == file_content);
}

// =============================================================================
// Streaming Session — threshold flush
// =============================================================================

static void test_session_basic_stream() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  {
    auto sess = w.open_session(os, /*threshold=*/1000);
    sess.add(hit("chr1", 100, Strand::Forward));
    sess.add(hit("chr1", 200, Strand::Reverse));
    // not yet flushed (below threshold), but close()/dtor will flush
    std::size_t n = sess.close();
    record("session: close() reports 2 written", n == 2u);
  }
  auto ls = lines(os.str());
  record("session: 3 lines (header + 2)", ls.size() == 3u,
         "lines=" + std::to_string(ls.size()));
}

static void test_session_auto_flush_at_threshold() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  auto sess = w.open_session(os, /*threshold=*/3);

  sess.add(hit("chr1", 1, Strand::Forward));
  sess.add(hit("chr1", 2, Strand::Forward));
  record("session: below threshold, 0 flushed", sess.total_written() == 0u);
  record("session: 2 buffered", sess.buffered() == 2u);

  sess.add(hit("chr1", 3, Strand::Forward)); // hits threshold → auto-flush
  record("session: at threshold, 3 flushed", sess.total_written() == 3u);
  record("session: buffer cleared after flush", sess.buffered() == 0u);

  sess.add(hit("chr1", 4, Strand::Forward));
  record("session: 1 buffered after next add", sess.buffered() == 1u);
  sess.close();
  record("session: 4 total after close", sess.total_written() == 4u);
}

static void test_session_memory_is_bounded() {
  // The key property: with threshold N, the buffer never exceeds N regardless
  // of how many records are written.
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  auto sess = w.open_session(os, /*threshold=*/10);

  std::size_t max_buffered = 0;
  for (int i = 1; i <= 95; ++i) {
    sess.add(hit("chr1", i, Strand::Forward));
    if (sess.buffered() > max_buffered)
      max_buffered = sess.buffered();
  }
  record("session: buffer never exceeds threshold", max_buffered <= 10u,
         "max_buffered=" + std::to_string(max_buffered));
  sess.close();
  record("session: all 95 records written", sess.total_written() == 95u);

  auto ls = lines(os.str());
  record("session: 96 lines (header + 95)", ls.size() == 96u,
         "lines=" + std::to_string(ls.size()));
}

static void test_session_add_all() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  auto sess = w.open_session(os, /*threshold=*/2);
  std::vector<OffTarget> recs = {
      hit("chr1", 1, Strand::Forward), hit("chr1", 2, Strand::Forward),
      hit("chr1", 3, Strand::Forward), hit("chr1", 4, Strand::Forward),
      hit("chr1", 5, Strand::Forward),
  };
  sess.add_all(recs);
  sess.close();
  record("session add_all: 5 written", sess.total_written() == 5u);
  record("session add_all: header + 5 lines", lines(os.str()).size() == 6u);
}

static void test_session_dtor_flushes() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  {
    auto sess = w.open_session(os, /*threshold=*/1000);
    sess.add(hit("chr9", 42, Strand::Reverse));
    // no explicit close — destructor must flush the single buffered record
  }
  auto ls = lines(os.str());
  record("session dtor: flushes buffered record (header + 1)", ls.size() == 2u,
         "lines=" + std::to_string(ls.size()));
}

static void test_session_threshold_zero_rejected() {
  OutputWriter w{OutputFormat::Tsv};
  std::ostringstream os;
  bool threw = false;
  try {
    auto sess = w.open_session(os, /*threshold=*/0);
    (void)sess;
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("session: threshold 0 throws invalid_argument", threw);
}

static void test_session_to_file_threshold() {
  OutputWriter w{OutputFormat::Targets};
  const std::string path = "/tmp/crispritz_ow_session.targets";
  {
    auto sess = w.open_session_to_file(path, /*threshold=*/4);
    for (int i = 1; i <= 10; ++i)
      sess.add(hit("chr1", i, Strand::Forward));
    std::size_t n = sess.close();
    record("session file: 10 written", n == 10u);
  } // file closed here via owned stream

  std::ifstream in(path);
  record("session file: exists", in.is_open());
  std::string content((std::istreambuf_iterator<char>(in)),
                      std::istreambuf_iterator<char>());
  in.close();
  record("session file: header + 10 lines", lines(content).size() == 11u,
         "lines=" + std::to_string(lines(content).size()));
  std::remove(path.c_str());
}

static void test_session_to_file_bad_path() {
  OutputWriter w{OutputFormat::Tsv};
  bool threw = false;
  try {
    auto sess = w.open_session_to_file("/nonexistent_dir_xyz/s.tsv");
    (void)sess;
  } catch (const std::runtime_error &) {
    threw = true;
  }
  record("session file: bad path throws runtime_error", threw);
}

static void test_session_equivalent_to_batch_write() {
  // Streaming output must be byte-identical to the one-shot write() for the
  // same records and format — chunking is an implementation detail.
  std::vector<OffTarget> recs;
  for (int i = 1; i <= 50; ++i)
    recs.push_back(
        hit("chr2", i, (i % 2) ? Strand::Forward : Strand::Reverse, i % 3));

  OutputWriter w{OutputFormat::Tsv};

  std::ostringstream batch_os;
  w.write(recs, batch_os);

  std::ostringstream stream_os;
  {
    auto sess = w.open_session(stream_os, /*threshold=*/7);
    sess.add_all(recs);
    sess.close();
  }

  record("streaming output == batch output (byte-identical)",
         batch_os.str() == stream_os.str());
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_output_writer ===\n\n";

  std::cout << "-- factory & construction --\n";
  test_factory_tsv();
  test_factory_targets();
  test_factory_bad_enum();
  test_writer_null_formatter_rejected();
  test_writer_format_name();

  std::cout << "\n-- valid output (TSV) --\n";
  test_tsv_header_and_row();

  std::cout << "\n-- valid output (Targets) --\n";
  test_targets_header_and_row();

  std::cout << "\n-- empty output --\n";
  test_empty_with_header();
  test_empty_without_header();
  test_no_header_flag_suppresses_header();

  std::cout << "\n-- multiple records --\n";
  test_multiple_records_vector();
  test_multiple_records_search_result();
  test_search_result_all_empty();

  std::cout << "\n-- file API & malformed path --\n";
  test_file_round_trip();
  test_malformed_path_throws();
  test_malformed_path_searchresult_throws();
  test_stream_file_consistency();

  std::cout << "\n-- streaming Session (threshold flush) --\n";
  test_session_basic_stream();
  test_session_auto_flush_at_threshold();
  test_session_memory_is_bounded();
  test_session_add_all();
  test_session_dtor_flushes();
  test_session_threshold_zero_rejected();
  test_session_to_file_threshold();
  test_session_to_file_bad_path();
  test_session_equivalent_to_batch_write();

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}