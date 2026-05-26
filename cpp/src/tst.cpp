#include "nucleotide_encoding.hpp"
#include "pam_search.hpp"
#include "tst.hpp"
#include "tst_utils.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace crispritz {

    TernarySearchTree::TernarySearchTree(const std::string& sequence, const std::string& chr_name,
                                        const std::string& pam_seq, int pam_length, int pam_limit,
                                        bool pam_at_start, int max_bulges, int num_threads)
        : sequence_(sequence), chr_name_(chr_name), pam_seq_(pam_seq),
        pam_length_(pam_length), pam_limit_(pam_limit), pam_at_start_(pam_at_start),
        guide_length_(pam_length - pam_limit), max_bulges_(max_bulges), num_threads_(num_threads)
    {
        if (pam_seq_.length() == static_cast<size_t>(pam_length_))
        {
            if (!pam_at_start_)
            {
                pam_rna_ = pam_seq_.substr(pam_length_ - pam_limit_, pam_limit_);
            }
            else
            {
                pam_rna_ = pam_seq_.substr(0, pam_limit_);
            }
        }
        else
        {
            pam_rna_ = pam_seq_;
        }
        std::transform(pam_rna_.begin(), pam_rna_.end(), pam_rna_.begin(), ::toupper);
    }

    TLeaf create_leaf(int idx, std::string target, std::string pam_dna) {
        TLeaf leaf;
        leaf.guide_index = idx;
        leaf.guide_dna = target;
        leaf.pam_dna = pam_dna;
        leaf.next = 0;
        return leaf;
    }

    void TernarySearchTree::build()
    {
        std::string upper_seq = sequence_;
        std::transform(upper_seq.begin(), upper_seq.end(), upper_seq.begin(), ::toupper);

        // search for PAM occurrences on input sequence
        pam::SearchParams params(pam_length_, pam_limit_, pam_at_start_, num_threads_);
        std::vector<int> pam_indexes = search_pam_sites(pam_rna_, upper_seq, params);

        leaves_.clear();
        leaves_.reserve(pam_indexes.size());

        size_t pamlen = static_cast<size_t>(pam_length_);
        size_t pam_lim = static_cast<size_t>(pam_limit_);

        if (pam_at_start_)
        {
            for (int idx : pam_indexes)
            {
                std::string target;
                if (idx < 0)
                {
                    size_t pos = static_cast<size_t>(-idx);
                    if (pos + pamlen + max_bulges_ > upper_seq.length()) 
                        continue;
                    target = upper_seq.substr(pos, pamlen + max_bulges_);
                    if (is_valid_genomic_window(target))
                        continue;
                    
                    // retrieve PAM sequence
                    std::string pam_dna = target.substr(0, pam_lim);
                    std::reverse(pam_dna.begin(), pam_dna.end());
                    // retrieve target sequence to encode in leaf
                    std::string target_dna = target.substr(pam_lim, pamlen - pam_lim + max_bulges_);

                    // add current target leaf to TST
                    TLeaf leaf = create_leaf(idx, target_dna, pam_dna);
                    leaves_.push_back(leaf);
                }
                else
                {
                    size_t pos = static_cast<size_t>(idx);
                    if (pos + pamlen + max_bulges_ > upper_seq.length()) 
                        continue;
                    target = upper_seq.substr(pos, pamlen + max_bulges_);
                    if (is_valid_genomic_window(target))
                        continue;
                    
                    // reverse complement target 
                    std::string tmp = iupac::reverse_complement(target);

                    // retrieve PAM sequence
                    std::string pam_dna = tmp.substr(0, pam_lim);
                    std::reverse(pam_dna.begin(), pam_dna.end());

                    // retrieve target sequence to encode in leaf
                    std:: string target_dna = tmp.substr(pam_lim, pamlen - pam_lim + max_bulges_);
                    
                    // add current target leaf to TST
                    TLeaf leaf = create_leaf(idx, target_dna, pam_dna);
                    leaves_.push_back(leaf);
                }

                
            }
        }
        else
        {
            for (int idx : pam_indexes)
            {
                std::string target;
                if (idx > 0)
                {
                    size_t pos = static_cast<size_t>(idx);
                    if (pos + pamlen + max_bulges_ > upper_seq.length()) 
                        continue;
                    target = upper_seq.substr(pos, pamlen + max_bulges_);
                    if (is_valid_genomic_window(target))
                        continue;

                    std::reverse(target.begin(), target.end());
                    
                    // retrieve PAM sequence
                    std::string pam_dna = target.substr(0, pam_lim);
                    // retrieve target sequence to encode in leaf
                    std::string target_dna = target.substr(pam_lim);

                    // add current target leaf to TST
                    TLeaf leaf = create_leaf(idx, target_dna, pam_dna);
                    leaves_.push_back(leaf);
                }
                else
                {
                    size_t pos = static_cast<size_t>(-idx);
                    if (pos + pamlen + max_bulges_ > upper_seq.length()) 
                        continue;
                    target = upper_seq.substr(pos, pamlen + max_bulges_);
                    if (is_valid_genomic_window(target))
                        continue;

                    // retrieve target sequence to encode in leaf
                    std::string tmp;
                    tmp.reserve(target.length());
                    for (char c : target)
                        tmp += iupac::complement(c);
                    std::string target_dna = tmp.substr(pam_lim);

                    // retrieve PAM sequence
                    std::string pam_dna = tmp.substr(0, pam_lim);

                    // add current target leaf to TST
                    TLeaf leaf = create_leaf(idx, target_dna, pam_dna);
                    leaves_.push_back(leaf);
                }
            }
        }

        leaves_.shrink_to_fit();  // remove unused cells

        std::sort(leaves_.begin(), leaves_.end(), [](const TLeaf& a, const TLeaf& b) {
            return a.guide_dna < b.guide_dna;
        });
    }

    void TernarySearchTree::insert_target(const std::string& target, int global_idx, int relative_idx,
                                         std::vector<TNode>& tree, int& node_used)
    {
        if (target.empty()) 
            return;
        const char* s = target.c_str();
        int current_node = 0;

        while (node_used > 0)
        {
            int d = *s - tree[current_node].splitchar;
            if (d == 0)  // null character or root
            {
                s++;
                if (*s == 0)
                {
                    leaves_[global_idx].next = tree[current_node].eqkid;
                    tree[current_node].eqkid = (relative_idx + 1) * -1;
                    return;
                }
                if (!tree[current_node].eqkid)
                {
                    tree[current_node].eqkid = node_used;
                    break;
                }
                current_node = tree[current_node].eqkid;
            }
            else if (d < 0)
            {
                if (!tree[current_node].lokid)
                {
                    tree[current_node].lokid = node_used;
                    break;
                }
                current_node = tree[current_node].lokid;
            }
            else
            {
                if (!tree[current_node].hikid)
                {
                    tree[current_node].hikid = node_used;
                    break;
                }
                current_node = tree[current_node].hikid;
            }
        }

        for (;;)
        {
            tree[node_used].splitchar = *s;
            tree[node_used].lokid = 0;
            tree[node_used].eqkid = 0;
            tree[node_used].hikid = 0;

            int prev_node = node_used;
            node_used++;

            s++;
            if (*s == 0)
            {
                tree[prev_node].eqkid = (relative_idx + 1) * -1;
                return;
            }
            tree[prev_node].eqkid = node_used;
        }
    }

    void TernarySearchTree::build_sub_tree(int l, int r, std::vector<TNode>& tree, int& node_used, int start)
    {
        if (r < l) return;
        int m = l + (r - l) / 2;
        int relative_idx = m - start;
        insert_target(leaves_[m].guide_dna, m, relative_idx, tree, node_used);
        build_sub_tree(l, m - 1, tree, node_used, start);
        build_sub_tree(m + 1, r, tree, node_used, start);
    }

    void TernarySearchTree::write_pair_nuc(std::ostream& os, char pair_nuc[2], uint8_t& bit_nuc) const
    {
        // serialize left and right nibble and write it to file
        bit_nuc = Serializer::serialize_left_nuc(pair_nuc[0]);
        bit_nuc += Serializer::serialize_right_nuc(pair_nuc[1]);
        os.put(bit_nuc);
    }

    bool TernarySearchTree::update_pair_nuc(char pair_nuc[2], bool switch_node, char nt, std::ostream& os, uint8_t& bit_nuc) const {
        if (switch_node) {
            pair_nuc[0] = nt;
            return false;
        }

        pair_nuc[1] = nt;
        write_pair_nuc(os, pair_nuc, bit_nuc);
        return true;
    }

    void TernarySearchTree::serialize_node(int node_idx, const std::vector<TNode>& tree, std::ostream& os,
                                           bool& switch_node, char pair_nuc[2], uint8_t& bit_nuc) const
    {
        const TNode& node = tree[node_idx];

        switch_node = update_pair_nuc(pair_nuc, switch_node, node.splitchar, os, bit_nuc);

        if (node.lokid > 0)  // lower characters lexicographically
            serialize_node(node.lokid, tree, os, switch_node, pair_nuc, bit_nuc);
        else 
            switch_node = update_pair_nuc(pair_nuc, switch_node, '0', os, bit_nuc);

        if (node.hikid > 0)  // higher characters lexicographically
            serialize_node(node.hikid, tree, os, switch_node, pair_nuc, bit_nuc);
        else
            switch_node = update_pair_nuc(pair_nuc, switch_node, '0', os, bit_nuc);

        if (node.eqkid > 0)  // equal characters lexicographically
            serialize_node(node.eqkid, tree, os, switch_node, pair_nuc, bit_nuc);
        else {
            pair_nuc[switch_node ? 0 : 1] = '_';
            write_pair_nuc(os, pair_nuc, bit_nuc);
            switch_node = true;
            os.write(reinterpret_cast<const char*>(&node.eqkid), sizeof(int));
        }
    }

    void TernarySearchTree::save(const std::string& outdir)
    {
        int counter_index = static_cast<int>(leaves_.size());
        if (counter_index == 0) return;

        // compute how many partitions must be built for current TST
        int group_tst = static_cast<int>(std::ceil(counter_index / (double)LEAVES_PER_GROUP));

        for (int jk = 0; jk < group_tst; jk++)
        {
            // get array for current partition
            int start = jk * LEAVES_PER_GROUP;
            int end = std::min(((jk + 1) * LEAVES_PER_GROUP), counter_index);
            int array_dim = end - start;

            std::vector<TNode> tree(array_dim * pam_length_);
            int nodes_used = 0;
            build_sub_tree(start, end - 1, tree, nodes_used, start);

            // open and write TST partition
            std::filesystem::path tst_bin = format_partition_filename(outdir, pam_seq_, chr_name_, jk + 1);
            std::ofstream file_tree(tst_bin, std::ios::out | std::ios::binary);
            if (!file_tree)
            {
                throw std::runtime_error("Failed to open output file for TST serialization: " + tst_bin.string());
            }

            file_tree.write(reinterpret_cast<const char*>(&array_dim), sizeof(int));
            file_tree.write(reinterpret_cast<const char*>(&guide_length_), sizeof(int));

            for (int i = start; i < end; i++)
            {
                file_tree.write(reinterpret_cast<const char*>(&leaves_[i].guide_index), sizeof(int));
                
                // serialize PAM sequence of each target
                uint8_t bit_nuc = 0;
                int k = 0;
                for (size_t char_idx = 0; char_idx < leaves_[i].pam_dna.length(); ++char_idx)
                {
                    char c = leaves_[i].pam_dna[char_idx];
                    bit_nuc += Serializer::serialize_pam(c);
                    k++;
                    bool is_last = (char_idx + 1 == leaves_[i].pam_dna.length());
                    if (is_last || k == 2)
                    {
                        file_tree.put(bit_nuc);
                        bit_nuc = 0;
                        k = 0;
                    }
                    if (bit_nuc)
                        bit_nuc <<= 4;
                }

                if (leaves_[i].next)
                {
                    file_tree.put('_');
                    file_tree.write(reinterpret_cast<const char*>(&leaves_[i].next), sizeof(int));
                }
                else
                {
                    file_tree.put('0');
                }
            }

            file_tree.write(reinterpret_cast<const char*>(&nodes_used), sizeof(int));

            if (nodes_used > 0)
            {
                bool flag = true;
                char pair_nuc[2] = {0, 0};
                uint8_t bit_nuc = 0;
                serialize_node(0, tree, file_tree, flag, pair_nuc, bit_nuc);
            }

            file_tree.close();  // close file
        }
    }
}





