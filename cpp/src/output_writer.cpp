#include "output_writer.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

namespace crispritz
{

    // =========================================================================
    // TsvFormatter
    // =========================================================================

    std::string TsvFormatter::header() const
    {
        // Reuse the canonical column list owned by OffTarget so the header and
        // the rows can never drift apart.
        const std::vector<std::string> cols = OffTarget::tsv_header();
        std::string line;
        for (std::size_t i = 0; i < cols.size(); ++i)
        {
            if (i != 0)
                line += '\t';
            line += cols[i];
        }
        return line;
    }

    std::string TsvFormatter::format_row(const OffTarget& ot) const
    {
        // OffTarget already knows how to serialize itself in the canonical
        // schema; the formatter must not re-implement it.
        return ot.to_tsv_row('\t');
    }

    // =========================================================================
    // TargetsFormatter
    // =========================================================================
    //
    // Legacy "targets" column order:
    //   bulge_type, grna, target, chrom, pos, strand,
    //   mismatches, bulge_size, total
    //
    // bulge_size = bulge_dna + bulge_rna
    // total      = mismatches + bulge_size  (== total_edit_distance)
    // =========================================================================

    std::string TargetsFormatter::header() const
    {
        return "bulge_type\tgrna\ttarget\tchrom\tpos\tstrand\t"
               "mismatches\tbulge_size\ttotal";
    }

    std::string TargetsFormatter::format_row(const OffTarget& ot) const
    {
        const int bulge_size = ot.bulge_dna() + ot.bulge_rna();
        const int total      = ot.total_edit_distance();

        std::string row;
        row.reserve(ot.grna().size() + ot.target().size() + 64u);

        row += ot.bulge_type();
        row += '\t';
        row += ot.grna();
        row += '\t';
        row += ot.target();
        row += '\t';
        row += ot.chrom();
        row += '\t';
        row += std::to_string(ot.pos());
        row += '\t';
        row += to_char(ot.strand());
        row += '\t';
        row += std::to_string(ot.mismatches());
        row += '\t';
        row += std::to_string(bulge_size);
        row += '\t';
        row += std::to_string(total);

        return row;
    }

    // =========================================================================
    // Formatter factory
    // =========================================================================

    std::unique_ptr<OffTargetFormatter> make_formatter(OutputFormat fmt)
    {
        switch (fmt)
        {
        case OutputFormat::Tsv:
            return std::make_unique<TsvFormatter>();
        case OutputFormat::Targets:
            return std::make_unique<TargetsFormatter>();
        }
        throw std::invalid_argument(
            "make_formatter: unrecognised OutputFormat value " +
            std::to_string(static_cast<int>(fmt)));
    }

    // =========================================================================
    // OutputWriter
    // =========================================================================

    OutputWriter::OutputWriter(OutputFormat fmt)
        : formatter_(make_formatter(fmt))
    {
    }

    OutputWriter::OutputWriter(std::unique_ptr<OffTargetFormatter> formatter)
        : formatter_(std::move(formatter))
    {
        if (!formatter_)
            throw std::invalid_argument("OutputWriter: formatter must not be null");
    }

    namespace
    {
        /**
         * @brief Throw if the stream is in a fail/bad state.
         * @param os    Stream to check.
         * @param where Context string for the error message.
         */
        void check_stream(const std::ostream& os, const char* where)
        {
            if (!os.good())
                throw std::runtime_error(std::string("OutputWriter: write failed (") +
                                         where + ')');
        }
    } // namespace

    std::size_t OutputWriter::write(const std::vector<OffTarget>& records,
                                    std::ostream&                 os,
                                    bool                          write_header) const
    {
        if (!os.good())
            throw std::runtime_error("OutputWriter: output stream is not writable");

        if (write_header)
        {
            os << formatter_->header() << '\n';
            check_stream(os, "header");
        }

        std::size_t written = 0;
        for (const OffTarget& ot : records)
        {
            os << formatter_->format_row(ot) << '\n';
            check_stream(os, "row");
            ++written;
        }
        return written;
    }

    std::size_t OutputWriter::write(const SearchResult& result,
                                    std::ostream&       os,
                                    bool                write_header) const
    {
        if (!os.good())
            throw std::runtime_error("OutputWriter: output stream is not writable");

        if (write_header)
        {
            os << formatter_->header() << '\n';
            check_stream(os, "header");
        }

        std::size_t written = 0;
        for (const std::vector<OffTarget>& guide_hits : result.hits_by_guide)
        {
            for (const OffTarget& ot : guide_hits)
            {
                os << formatter_->format_row(ot) << '\n';
                check_stream(os, "row");
                ++written;
            }
        }
        return written;
    }

    std::size_t OutputWriter::write_to_file(const std::vector<OffTarget>& records,
                                            const std::string&            path) const
    {
        std::ofstream out(path, std::ios::out | std::ios::trunc);
        if (!out.is_open())
            throw std::runtime_error("OutputWriter: cannot open output file: " + path);

        // RAII: `out` closes on scope exit, including on the exception paths
        // inside write().
        return write(records, out, /*write_header=*/true);
    }

    std::size_t OutputWriter::write_to_file(const SearchResult& result,
                                            const std::string&  path) const
    {
        std::ofstream out(path, std::ios::out | std::ios::trunc);
        if (!out.is_open())
            throw std::runtime_error("OutputWriter: cannot open output file: " + path);

        return write(result, out, /*write_header=*/true);
    }

} // namespace crispritz