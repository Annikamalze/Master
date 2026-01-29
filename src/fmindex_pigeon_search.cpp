#include <algorithm>
#include <filesystem>
#include <fstream>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <cereal/archives/binary.hpp>


std::vector<std::vector<seqan3::dna5>>
cut_query(size_t pieces,
          std::vector<seqan3::dna5> const & query)
{
    std::vector<std::vector<seqan3::dna5>> parts;
    size_t l = query.size() / pieces;

    for (size_t i = 0; i < pieces; ++i) {
        if (i == pieces - 1) {
            parts.emplace_back(query.begin() + l * i, query.end());
        } else {
            parts.emplace_back(query.begin() + l * i,
                               query.begin() + l * (i + 1));
        }
    }

    return parts;
}

bool verify(std::vector<seqan3::dna5> const & ref,
            std::vector<seqan3::dna5> const & query,
            size_t pos,
            size_t max_errors)
{
    if (pos + query.size() > ref.size())
        return false;

    size_t errors = 0;

    for (size_t i = 0; i < query.size(); ++i) {
        if (ref[pos + i] != query[i]) {
            ++errors;
            if (errors > max_errors)
                return false;
        }
    }

    return true;
}

int main(int argc, char const * const * argv)
{
    seqan3::argument_parser parser{
        "fmindex_pigeon_search",
        argc,
        argv,
        seqan3::update_notifications::off
    };

    parser.info.author  = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path index_path;
    std::filesystem::path reference_file;
    std::filesystem::path query_file;

    size_t  number_of_queries = 100;
    uint8_t number_of_errors  = 0;

    parser.add_option(index_path,        '\0', "index",    "path to the fm-index");
    parser.add_option(reference_file,    '\0', "reference","path to the reference file");
    parser.add_option(query_file,        '\0', "query",    "path to the query file");
    parser.add_option(number_of_queries, '\0', "query_ct", "number of queries");
    parser.add_option(number_of_errors,  '\0', "errors",   "allowed hamming distance");

    try {
        parser.parse();
    } catch (seqan3::argument_parser_error const & e) {
        seqan3::debug_stream << "Parsing error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    seqan3::sequence_file_input reference_stream{reference_file};
    seqan3::sequence_file_input query_stream{query_file};

    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto & rec : reference_stream)
        reference.push_back(rec.sequence());

    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto & rec : query_stream)
        queries.push_back(rec.sequence());

    using Index = seqan3::fm_index<std::vector<std::vector<seqan3::dna5>>>;
    Index index;

    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive archive{is};
        archive(index);
        seqan3::debug_stream << "done\n";
    }

    while (queries.size() < number_of_queries) {
        size_t old = queries.size();
        queries.resize(2 * old);
        std::copy_n(queries.begin(), old, queries.begin() + old);
    }
    queries.resize(number_of_queries);

    seqan3::configuration const cfg =
        seqan3::search_cfg::max_error_total{
            seqan3::search_cfg::error_count{number_of_errors}
        };

    std::vector<std::vector<std::pair<size_t, size_t>>> results;

    for (size_t q = 0; q < queries.size(); ++q) {
        auto const & query = queries[q];
        auto parts = cut_query(3, query);

        std::vector<std::pair<size_t, size_t>> query_hits;

        for (size_t p = 0; p < parts.size(); ++p) {
            size_t offset = p * parts[0].size();

            for (auto const & hit : seqan3::search(parts[p], index, cfg)) {
                size_t ref_id  = hit.reference_id();
                size_t hit_pos = hit.reference_begin_position();

                if (hit_pos < offset)
                    continue;

                size_t candidate = hit_pos - offset;

                if (verify(reference[ref_id], query, candidate, number_of_errors)) {
                    query_hits.emplace_back(ref_id, candidate);
                }
            }
        }

        std::sort(query_hits.begin(), query_hits.end());
        query_hits.erase(std::unique(query_hits.begin(), query_hits.end()),
                         query_hits.end());

        results.push_back(std::move(query_hits));
    }

    return 0;
}
