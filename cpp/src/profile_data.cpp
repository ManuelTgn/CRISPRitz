#include "profile_data.hpp"

#include <cctype>      // std::toupper
#include <stdexcept>
#include <string>
#include <vector>

namespace crispritz {

// =============================================================================
// ProfileAccumulator — constructor
// =============================================================================

ProfileAccumulator::ProfileAccumulator(std::string guide,
                                       int         guide_len,
                                       int         pam_len,
                                       int         max_mm,
                                       int         max_bulge_dna,
                                       int         max_bulge_rna,
                                       bool        pam_at_start)
    : guide_(std::move(guide))
    , guide_len_(guide_len)
    , pam_len_(pam_len)
    , max_mm_(max_mm)
    , max_bulge_dna_(max_bulge_dna)
    , max_bulge_rna_(max_bulge_rna)
    , pam_at_start_(pam_at_start)
{
    if (guide_len_ <= 0)
        throw std::invalid_argument(
            "ProfileAccumulator: guide_len must be > 0, got " +
            std::to_string(guide_len_));
    if (pam_len_ < 0)
        throw std::invalid_argument(
            "ProfileAccumulator: pam_len must be >= 0, got " +
            std::to_string(pam_len_));
    if (max_mm_ < 0)
        throw std::invalid_argument(
            "ProfileAccumulator: max_mm must be >= 0, got " +
            std::to_string(max_mm_));
    if (max_bulge_dna_ < 0)
        throw std::invalid_argument(
            "ProfileAccumulator: max_bulge_dna must be >= 0, got " +
            std::to_string(max_bulge_dna_));
    if (max_bulge_rna_ < 0)
        throw std::invalid_argument(
            "ProfileAccumulator: max_bulge_rna must be >= 0, got " +
            std::to_string(max_bulge_rna_));

    const int mm_slots  = max_mm_ + 1;
    const int dna_slots = (max_bulge_dna_ > 0) ? max_bulge_dna_ : 1;
    const int rna_slots = (max_bulge_rna_ > 0) ? max_bulge_rna_ : 1;

    // File 1
    pos_mm_count_.assign(guide_len_, 0);
    offt_by_mm_.assign(mm_slots, 0);

    // File 2
    ext_mm_nuc_pos_.assign(mm_slots,
        std::vector<std::vector<int>>(4,
            std::vector<int>(guide_len_, 0)));
    ext_total_by_mm_.assign(mm_slots, std::vector<int>(guide_len_, 0));
    ext_dna_by_mm_pos_.assign(mm_slots, std::vector<int>(guide_len_, 0));
    ext_rna_by_mm_pos_.assign(mm_slots, std::vector<int>(guide_len_, 0));

    // File 3
    pos_bulge_dna_.assign(guide_len_, 0);
    pos_mm_in_dna_.assign(guide_len_, 0);
    offt_dna_.assign(mm_slots, std::vector<int>(dna_slots, 0));

    // File 4
    pos_bulge_rna_.assign(guide_len_, 0);
    pos_mm_in_rna_.assign(guide_len_, 0);
    offt_rna_.assign(mm_slots, std::vector<int>(rna_slots, 0));

    // File 5
    offt_complete_by_mm_.assign(mm_slots, 0);
    pos_mm_complete_.assign(guide_len_, 0);
}

// =============================================================================
// ProfileAccumulator::nuc_index  (static helper)
// =============================================================================

int ProfileAccumulator::nuc_index(char c) noexcept
{
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(c)))) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:  return -1; // 'N' or any ambiguous base
    }
}

// =============================================================================
// ProfileAccumulator::push
// =============================================================================

void ProfileAccumulator::push(const OffTarget& ot)
{
    const std::string& aln_grna   = ot.grna();
    const std::string& aln_target = ot.target();
    const std::size_t  aln_len    = aln_grna.size();

    const int mm_in_hit  = ot.mismatches();
    const int bd_in_hit  = ot.bulge_dna();
    const int br_in_hit  = ot.bulge_rna();
    const bool has_dna   = bd_in_hit > 0;
    const bool has_rna   = br_in_hit > 0;
    const bool has_bulge = has_dna || has_rna;

    // Clamp mm_in_hit to the configured budget when indexing into counters.
    // Hits beyond the budget should not appear (the searcher already enforces
    // this), but defensive clamping avoids undefined behaviour if they do.
    const int mm_idx = (mm_in_hit <= max_mm_) ? mm_in_hit : max_mm_;

    // -----------------------------------------------------------------------
    // Walk the alignment string column by column.
    //
    // body_pos tracks the current *guide body* position (0-based, PAM
    // columns excluded).  It advances on every non-gap guide column that
    // is *not* a PAM column.
    //
    // PAM columns are those at the leading (pam_at_start_) or trailing
    // (!pam_at_start_) end of the alignment.  Because the PAM length is
    // fixed and the alignment has no gaps in the PAM region (the PAM was
    // the filter that caused this site to be indexed), we can identify PAM
    // columns by their position in the raw alignment string — not by the
    // body_pos counter.
    //
    // The alignment string layout is:
    //   pam_at_start_:   [pam_len_ PAM cols] [guide_len_ + bulge cols]
    //   !pam_at_start_:  [guide_len_ + bulge cols] [pam_len_ PAM cols]
    //
    // We simply skip (advance past) PAM columns without touching any counter.
    // -----------------------------------------------------------------------

    // Precompute the half-open interval [body_start, body_end) in the raw
    // alignment that covers the guide body (including any bulge gaps).
    const std::size_t body_start = pam_at_start_
                                       ? static_cast<std::size_t>(pam_len_)
                                       : 0u;
    const std::size_t body_end   = pam_at_start_
                                       ? aln_len
                                       : aln_len - static_cast<std::size_t>(pam_len_);

    int body_pos = 0; // current guide body position (PAM-stripped)

    for (std::size_t i = body_start; i < body_end; ++i) {
        const char g = aln_grna[i];
        const char t = aln_target[i];

        if (g == '-') {
            // ── DNA bulge ─────────────────────────────────────────────────
            // The target has an extra base; the guide has a gap.  We do NOT
            // advance body_pos (this column does not correspond to a guide
            // position), but we record the bulge at the *current* body_pos
            // (the position the guide is about to match or has just matched).
            //
            // Guard body_pos within range: when a DNA bulge appears at the
            // very start of the body section, body_pos is 0 — still valid.
            if (body_pos < guide_len_) {
                pos_bulge_dna_[body_pos]++;
                ext_dna_by_mm_pos_[mm_idx][body_pos]++;
            }
            // body_pos does NOT advance for DNA bulge columns.

        } else if (t == '-') {
            // ── RNA bulge ─────────────────────────────────────────────────
            // The guide has an extra base; the target has a gap.
            if (body_pos < guide_len_) {
                pos_bulge_rna_[body_pos]++;
                ext_rna_by_mm_pos_[mm_idx][body_pos]++;
            }
            // body_pos advances: this alignment column consumed a guide base.
            ++body_pos;

        } else {
            // ── Match or substitution mismatch ────────────────────────────
            const char gu = static_cast<char>(std::toupper(static_cast<unsigned char>(g)));
            const char tu = static_cast<char>(std::toupper(static_cast<unsigned char>(t)));

            const bool is_mismatch = (gu != tu) && (gu != 'N') && (tu != 'N');

            if (is_mismatch && body_pos < guide_len_) {
                // File 5 (complete channel): every mismatch, regardless of bulge.
                pos_mm_complete_[body_pos]++;

                // File 1 (MM-only channel): only mismatches from no-bulge hits.
                // A mismatch occurring inside a bulge alignment belongs to the
                // DNA/RNA channels, not the mismatch-only profile.
                if (!has_bulge)
                    pos_mm_count_[body_pos]++;

                // Files 3/4: positional mismatch within a bulge alignment.
                if (has_dna) pos_mm_in_dna_[body_pos]++;
                if (has_rna) pos_mm_in_rna_[body_pos]++;

                // File 2: per-(mm_threshold, nucleotide, position).
                const int nuc = nuc_index(tu); // off-target nucleotide
                if (nuc >= 0) {
                    ext_mm_nuc_pos_[mm_idx][nuc][body_pos]++;
                    ext_total_by_mm_[mm_idx][body_pos]++;
                }
            }
            // body_pos always advances for non-gap guide columns.
            ++body_pos;
        }
    }

    // -----------------------------------------------------------------------
    // Bucket the hit into per-mm-count / per-bulge-size summary cells.
    // -----------------------------------------------------------------------

    // File 5: complete channel — every hit contributes here.
    // "On-target" in the complete sense means 0 mismatches in ANY channel
    // (MM-only, DNA-bulge, or RNA-bulge), so the bulge state is irrelevant.
    offt_complete_by_mm_[mm_idx]++;
    if (mm_idx == 0) {
        ont_count_complete_++;
    }

    if (!has_bulge) {
        // ── MM-only channel (File 1) ─────────────────────────────────────
        offt_by_mm_[mm_idx]++;
        if (mm_in_hit == 0) {
            ont_count_++;
        }
    }

    if (has_dna && max_bulge_dna_ > 0) {
        // ── DNA-bulge channel (File 3) ───────────────────────────────────
        const int bd_idx = bd_in_hit - 1; // bulge size 1 → index 0
        if (bd_idx < max_bulge_dna_) {
            offt_dna_[mm_idx][bd_idx]++;
        }
        if (mm_in_hit == 0) {
            ont_count_dna_++;
        }
    }

    if (has_rna && max_bulge_rna_ > 0) {
        // ── RNA-bulge channel (File 4) ───────────────────────────────────
        const int br_idx = br_in_hit - 1;
        if (br_idx < max_bulge_rna_) {
            offt_rna_[mm_idx][br_idx]++;
        }
        if (mm_in_hit == 0) {
            ont_count_rna_++;
        }
    }
}

// =============================================================================
// ProfileAccumulator::build
// =============================================================================

GuideProfile ProfileAccumulator::build() const
{
    GuideProfile p;

    p.guide           = guide_;
    p.guide_len       = guide_len_;
    p.pam_len         = pam_len_;
    p.pam_at_start    = pam_at_start_;
    p.max_mm          = max_mm_;
    p.max_bulge_dna   = max_bulge_dna_;
    p.max_bulge_rna   = max_bulge_rna_;

    // File 1
    p.pos_mm_count    = pos_mm_count_;
    p.ont_count       = ont_count_;
    p.offt_by_mm      = offt_by_mm_;

    // File 2
    p.ext_mm_nuc_pos   = ext_mm_nuc_pos_;
    p.ext_total_by_mm  = ext_total_by_mm_;
    p.ext_dna_by_mm_pos = ext_dna_by_mm_pos_;
    p.ext_rna_by_mm_pos = ext_rna_by_mm_pos_;

    // File 3
    p.pos_bulge_dna   = pos_bulge_dna_;
    p.pos_mm_in_dna   = pos_mm_in_dna_;
    p.offt_dna        = offt_dna_;
    p.ont_count_dna   = ont_count_dna_;

    // File 4
    p.pos_bulge_rna   = pos_bulge_rna_;
    p.pos_mm_in_rna   = pos_mm_in_rna_;
    p.offt_rna        = offt_rna_;
    p.ont_count_rna   = ont_count_rna_;

    // File 5
    p.offt_complete_by_mm  = offt_complete_by_mm_;
    p.ont_count_complete   = ont_count_complete_;
    p.pos_mm_complete      = pos_mm_complete_;

    return p;
}

// =============================================================================
// ProfileAccumulator::reset
// =============================================================================

void ProfileAccumulator::reset()
{
    const int mm_slots  = max_mm_ + 1;
    const int dna_slots = (max_bulge_dna_ > 0) ? max_bulge_dna_ : 1;
    const int rna_slots = (max_bulge_rna_ > 0) ? max_bulge_rna_ : 1;

    // File 1
    std::fill(pos_mm_count_.begin(), pos_mm_count_.end(), 0);
    ont_count_ = 0;
    std::fill(offt_by_mm_.begin(), offt_by_mm_.end(), 0);

    // File 2
    for (int m = 0; m < mm_slots; ++m) {
        for (int n = 0; n < 4; ++n)
            std::fill(ext_mm_nuc_pos_[m][n].begin(), ext_mm_nuc_pos_[m][n].end(), 0);
        std::fill(ext_total_by_mm_[m].begin(),   ext_total_by_mm_[m].end(),   0);
        std::fill(ext_dna_by_mm_pos_[m].begin(), ext_dna_by_mm_pos_[m].end(), 0);
        std::fill(ext_rna_by_mm_pos_[m].begin(), ext_rna_by_mm_pos_[m].end(), 0);
    }

    // File 3
    std::fill(pos_bulge_dna_.begin(), pos_bulge_dna_.end(), 0);
    std::fill(pos_mm_in_dna_.begin(), pos_mm_in_dna_.end(), 0);
    for (int m = 0; m < mm_slots; ++m)
        std::fill(offt_dna_[m].begin(), offt_dna_[m].end(), 0);
    ont_count_dna_ = 0;

    // File 4
    std::fill(pos_bulge_rna_.begin(), pos_bulge_rna_.end(), 0);
    std::fill(pos_mm_in_rna_.begin(), pos_mm_in_rna_.end(), 0);
    for (int m = 0; m < mm_slots; ++m)
        std::fill(offt_rna_[m].begin(), offt_rna_[m].end(), 0);
    ont_count_rna_ = 0;

    // File 5
    std::fill(offt_complete_by_mm_.begin(), offt_complete_by_mm_.end(), 0);
    ont_count_complete_ = 0;
    std::fill(pos_mm_complete_.begin(), pos_mm_complete_.end(), 0);

    (void)dna_slots; (void)rna_slots; // used only at construction
}

} // namespace crispritz