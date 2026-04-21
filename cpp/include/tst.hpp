#pragma once

#include <string>

namespace crispritz
{
    void build_tree(const std::string& sequence, const std::string& pam, int pam_length,
                    int pam_limit, bool pam_at_start, int num_threads = 1);
}
