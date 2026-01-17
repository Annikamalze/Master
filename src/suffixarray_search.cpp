#include <divsufsort.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

int compare_suffix(std::vector<seqan3::dna5> const & text,
                   size_t sa_pos,
                   std::vector<seqan3::dna5> const & pattern)
{
    size_t i = 0;
    while (i < pattern.size() && sa_pos + i < text.size())
    {
        if (text[sa_pos + i] < pattern[i]) return -1;
        if (text[sa_pos + i] > pattern[i]) return 1;
        ++i;
    }

    if (i == pattern.size())
        return 0; // full match

    return -1; // suffix shorter than pattern
}

size_t find_LP(std::vector<seqan3::dna5> const & text,
               std::vector<saidx_t> const & sa,
               std::vector<seqan3::dna5> const & pattern)
{
    size_t L = 0, R = sa.size();
    while (L < R)
    {
        size_t M = (L + R) / 2;
        if (compare_suffix(text, sa[M], pattern) >= 0)
            R = M;
        else
            L = M + 1;
    }
    return L;
}

size_t find_RP(std::vector<seqan3::dna5> const & text,
               std::vector<saidx_t> const & sa,
               std::vector<seqan3::dna5> const & pattern)
{
    size_t L = 0, R = sa.size();
    while (L < R)
    {
        size_t M = (L + R) / 2;
        if (compare_suffix(text, sa[M], pattern) > 0)
            R = M;
        else
            L = M + 1;
    }
    return L;
}

int main(int argc, char const * const * argv)
{
    seqan3::argument_parser parser{
        "suffixarray_search", argc, argv, seqan3::update_notifications::off};

    std::filesystem::path reference_file{};
    std::filesystem::path query_file{};
    size_t number_of_queries{100};

    parser.add_option(reference_file, '\0', "reference", "path to the reference file");
    parser.add_option(query_file, '\0', "query", "path to the query file");
    parser.add_option(number_of_queries, '\0', "query_ct",
                      "number of queries (duplicated if needed)");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & e)
    {
        std::cerr << "Parsing error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    seqan3::sequence_file_input reference_stream{reference_file};
    std::vector<seqan3::dna5> reference;

    for (auto & record : reference_stream)
    {
        auto const & seq = record.sequence();
        reference.insert(reference.end(), seq.begin(), seq.end());
    }

    seqan3::sequence_file_input query_stream{query_file};
    std::vector<std::vector<seqan3::dna5>> queries;

    for (auto & record : query_stream)
        queries.push_back(record.sequence());

    while (queries.size() < number_of_queries)
    {
        size_t old_size = queries.size();
        queries.resize(2 * old_size);
        std::copy_n(queries.begin(), old_size, queries.begin() + old_size);
    }
    queries.resize(number_of_queries);

 
    std::vector<saidx_t> suffixarray(reference.size());
    sauchar_t const * str =
        reinterpret_cast<sauchar_t const *>(reference.data());

    if (divsufsort(str, suffixarray.data(), reference.size()) != 0)
    {
        std::cerr << "Suffix array construction failed\n";
        return EXIT_FAILURE;
    }


    for (auto const & q : queries)
    {
        size_t LP = find_LP(reference, suffixarray, q);
        size_t RP = find_RP(reference, suffixarray, q);

        std::vector<size_t> hits;
        for (size_t i = LP; i < RP; ++i)
            hits.push_back(suffixarray[i]);

        // optional output
        // std::cout << hits.size() << "\n";
    }

    return 0;
}
