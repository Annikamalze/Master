#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <vector>
#include <string>

std::vector<std::vector<seqan3::dna5>>
cut_query(int pieces, std::vector<seqan3::dna5> const & query)
{
    std::vector<std::vector<seqan3::dna5>> pieced_queries;

    size_t l = query.size() / pieces;

    for (int i = 0; i < pieces; i++)
    {
        if (i == pieces - 1)
        {
            pieced_queries.emplace_back(query.begin() + l * i, query.end());
        }
        else
        {
            pieced_queries.emplace_back(query.begin() + l * i,
                                        query.begin() + l * (i + 1));
        }
    }

    return pieced_queries;
}


bool verify(std::vector<seqan3::dna5> const & ref,
            std::vector<seqan3::dna5> const & query,
            size_t pos,
            size_t max_errors)
{
    // Bounds-Check
    if (pos + query.size() > ref.size())
        return false;

    size_t errors = 0;

    for (size_t i = 0; i < query.size(); i++) {   // <- i++
        if (ref[pos + i] != query[i]) {
            errors++;

            if (errors > max_errors)
                return false; // early exit
        }
    }

    return true;
}


int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed hamming distance errors");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read reference into memory
    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto& record : reference_stream) {
        reference.push_back(record.sequence());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // loading fm-index into memory
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Some hack
    Index index; // construct fm-index
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        seqan3::debug_stream << "done\n";
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    // split

    // create index
    // verify query in text (where is it in original/ not pieced text)
    seqan3::configuration const cfg =
        seqan3::search_cfg::max_error_total{
            seqan3::search_cfg::error_count{number_of_errors}
        };

    std::vector<std::vector<std::pair<size_t, size_t>>> results;

    for (size_t q = 0; q < queries.size(); q++)
    {   
        auto const & query = queries[q];

        // Teile die Query in 3 Stücke (Pigeonhole-Prinzip)
        auto parts = cut_query(3, query);

        std::vector<std::pair<size_t, size_t>> query_hits;

        size_t offset = 0; // Kumulativer Offset für die Teile
        for (size_t p = 0; p < parts.size(); ++p)
        {
            for (auto const & hit : seqan3::search(parts[p], index, cfg))
            {
                size_t ref_id  = hit.reference_id();
                size_t hit_pos = hit.reference_begin_position();

                if (hit_pos < offset)
                    continue;

                size_t candidate = hit_pos - offset;

                // verify prüft, ob candidate gültig ist
                if (verify(reference[ref_id], query, candidate, number_of_errors))
                {
                    query_hits.emplace_back(ref_id, candidate);
                }
            }

            // Offset für das nächste Stück aufsummieren
            offset += parts[p].size();
        }

        // Treffer sortieren & duplizierte entfernen
        std::sort(query_hits.begin(), query_hits.end());
        query_hits.erase(std::unique(query_hits.begin(), query_hits.end()),
                     query_hits.end());

        results.push_back(std::move(query_hits));
    }
}