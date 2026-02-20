#include "tst.hpp"

#include "pam_search.hpp"

#include <iostream>
#include <omp.h>

namespace crispritz
{

    void build_tree(const std::string& sequence, const std::string& pam, int pam_length,
                    int pam_limit, bool pam_at_start)
    {
        // set pam search parameters
        pam::SearchParams params(pam_length, pam_limit, pam_at_start);

        // convert genomic sequence in bits
        pam::CompactGenome genome_bits(sequence);

        // perform pam search on bit encoded genome sequence
        auto sites = pam::search_pam_sites_fast(pam, genome_bits, params);

        int counter = 0;
        for (int i; i < sites.size(); i++)
        {
            if (sites[i] < 0)
            {
                std::cout << sites[i] << "\n";
                counter += 1;
            }
        }

        std::cout << "Found " << counter << " sites on forward";
    }

} // namespace crispritz
