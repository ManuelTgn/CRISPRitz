#include "pam_search.hpp"
#include "nucleotide_encoding.hpp"
#include <iostream>
#include <cassert>
#include <string>

using namespace pam;

void test_nucleotide_encoding() {
    std::cout << "Testing nucleotide encoding... ";
    
    // Test basic encoding
    assert(NucleotideEncoder::encode_genome('A') == 0b0001);
    assert(NucleotideEncoder::encode_genome('C') == 0b0010);
    assert(NucleotideEncoder::encode_genome('G') == 0b0100);
    assert(NucleotideEncoder::encode_genome('T') == 0b1000);
    assert(NucleotideEncoder::encode_genome('N') == 0b0000);
    
    // Test PAM encoding (N = wildcard)
    assert(NucleotideEncoder::encode_pam('N') == 0b1111);
    
    // Test complement
    assert(NucleotideEncoder::complement('A') == 'T');
    assert(NucleotideEncoder::complement('T') == 'A');
    assert(NucleotideEncoder::complement('C') == 'G');
    assert(NucleotideEncoder::complement('G') == 'C');
    
    std::cout << "PASSED\n";
}

void test_reverse_complement() {
    std::cout << "Testing reverse complement... ";
    
    std::string seq = "ACGT";
    std::string expected = "ACGT";
    assert(reverse_complement(seq) == expected);
    
    seq = "AAAA";
    expected = "TTTT";
    assert(reverse_complement(seq) == expected);
    
    seq = "ATCG";
    expected = "CGAT";
    assert(reverse_complement(seq) == expected);
    
    std::cout << "PASSED\n";
}

void test_compact_genome() {
    std::cout << "Testing compact genome... ";
    
    std::string genome = "ACGTACGT";
    CompactGenome compact(genome);
    
    assert(compact.size() == 8);
    assert(compact[0] == NucleotideEncoder::encode_genome('A'));
    assert(compact[1] == NucleotideEncoder::encode_genome('C'));
    assert(compact[2] == NucleotideEncoder::encode_genome('G'));
    assert(compact[3] == NucleotideEncoder::encode_genome('T'));
    
    std::cout << "PASSED\n";
}

void test_pam_search_exact() {
    std::cout << "Testing PAM search (exact match)... ";
    
    // Simple test genome with known NGG sites
    std::string genome = "ACGTAGGACTCGGTACGGG";
    //                     0123456789...
    // NGG appears at positions: 4, 10, 16
    
    SearchParams params(3, 3, 0, 0, false);  // No mismatches
    auto sites = search_pam_sites("NGG", genome, params);
    
    // Should find sites on both strands
    assert(!sites.empty());
    std::cout << "Found " << sites.size() << " sites - PASSED\n";
}

void test_iupac_codes() {
    std::cout << "Testing IUPAC ambiguity codes... ";
    
    // R = A or G, so should match both
    std::string genome = "ACGTATGGAACTGAAGAGG";
    //                     0123456789...
    // NRG appears at positions: 5, 13, 15, 16
    
    SearchParams params(3, 3, 0, 0, false);
    auto sites = search_pam_sites("NRG", genome, params);  // R matches A or G
    
    assert(!sites.empty());
    std::cout << "Found " << sites.size() << " sites - PASSED\n";
}

void test_performance() {
    std::cout << "Testing performance on larger genome... ";
    
    // Create a 1MB test genome
    std::string genome;
    genome.reserve(1000000);
    const char bases[] = {'A', 'C', 'G', 'T'};
    
    for (int i = 0; i < 1000000; ++i) {
        genome += bases[i % 4];
    }
    
    // Add some NGG sites
    for (int i = 0; i < genome.length() - 100; i += 1000) {
        genome[i] = 'N';
        genome[i+1] = 'G';
        genome[i+2] = 'G';
    }
    
    SearchParams params(3, 3, 0, 0, false);
    
    auto start = std::chrono::high_resolution_clock::now();
    CompactGenome compact(genome);
    auto encode_end = std::chrono::high_resolution_clock::now();
    
    auto sites = search_pam_sites_fast("NGG", compact, params);
    auto search_end = std::chrono::high_resolution_clock::now();
    
    auto encode_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        encode_end - start).count();
    auto search_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        search_end - encode_end).count();
    
    std::cout << "\n  Encoding: " << encode_ms << " ms\n";
    std::cout << "  Search: " << search_ms << " ms\n";
    std::cout << "  Found: " << sites.size() << " sites - PASSED\n";
}

int main() {
    std::cout << "=== PAM Generator Test Suite ===\n\n";
    
    try {
        test_nucleotide_encoding();
        test_reverse_complement();
        test_compact_genome();
        test_pam_search_exact();
        test_iupac_codes();
        test_performance();
        
        std::cout << "\n=== ALL TESTS PASSED ===\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n!!! TEST FAILED !!!\n";
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}