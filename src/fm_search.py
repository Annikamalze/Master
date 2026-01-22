# Annika Malze (5565183), Clara Daßio (5568983), Janina Krentz (5566628)
import iv2py as iv
import time

def number_of_queries(queries_orig, number_wanted):
    if len(queries_orig) >= number_wanted:
        queries = queries_orig[:number_wanted] # auf gewünschte Anzahl kürzen
        return queries  
    else:
        # Duplizieren bis gewünschte Anzahl
        extended = queries_orig.copy()
        while len(extended) < number_wanted:
            needed = number_wanted - len(extended)
            extended += queries_orig[:needed]
        return extended


def run_query_search(reference, n_queries, n):
    """
    function constructs FMIndex and searches for n specific queries in a given query list
    input:
    reference: text to be searched
    n_queries: list of queries
    n: number of queries from the query list to be searched
    output:
    list of lists with occurences of the respective queries
    """
    start_index = time.perf_counter()
    index = iv.fmindex(reference=reference, samplingRate=16)
    start = time.perf_counter()

    res = []
    for i in range(n):
        exists = index.search(n_queries[i])
        res.append(exists[0][1] if exists else []) # index search returns list of one tuple

    end = time.perf_counter()

    runtime = end - start
    runtime_index = end - start_index
    print(f"runtime for {n} queries: {runtime}")
    print(f"runtime for {n} queries with fm construction: {runtime_index}")
    print(f"first 10 matches: {res[:10]}")
    return res


reference = []
for record in iv.fasta.reader(file="../data/hg38_partial.fasta.gz"):
    reference.append(str(record.seq))

queries_100 = []
for record in iv.fasta.reader(file="../data/illumina_reads_100.fasta.gz"):
    queries_100.append(str(record.seq))

# tests
# reference = ["PINEAPPLEPEN"]
# index = iv.fmindex(reference=reference, samplingRate=16)
# result = index.search(queries_100[0])

# result = run_query_search(reference, queries_100, 10)


############### exercise 4 #######################
# [10**3, 10**4, 10**5, 10**6]
n = 10**6
n_queries = number_of_queries(queries_100, n)
res = run_query_search(reference, n_queries, n)


############### exercise 5 #######################
# print(type(reference), type(reference[0])) # 455, list of strings
whole_genome = []
for record in iv.fasta.reader(file="../data/GRCh38_genomic.fna.gz"):
    whole_genome.append(str(record.seq))

def run_query_search_on_genome(reference, n_queries, n):
    """
    function constructs FMIndex and searches for n specific queries in a given query list
    input:
    reference: text to be searched
    n_queries: list of queries
    n: number of queries from the query list to be searched
    output:
    list of lists of tuples with occurences of the respective queries in each reference sequence
    """
    start_index = time.perf_counter()
    index = iv.fmindex(reference=reference, samplingRate=16)
    start = time.perf_counter()

    res = []
    for i in range(n):
        res.append(index.search(n_queries[i])) # return tuple, need to know from which reference it is

    end = time.perf_counter()

    runtime = end - start
    runtime_index = end - start_index
    print(f"runtime for {n} queries: {runtime}")
    print(f"runtime for {n} queries with fm construction: {runtime_index}")
    print(f"first 10 matches: {res[:10]}")
    return res


####### length 40 ############
queries_l40 = []
for record in iv.fasta.reader(file="../data/illumina_reads_40.fasta.gz"):
    queries_l40.append(str(record.seq))

# print(len(queries_l40)) # 100.000

####### length 60 ############
queries_l60 = []
for record in iv.fasta.reader(file="../data/illumina_reads_60.fasta.gz"):
    queries_l60.append(str(record.seq))

# print(len(queries_l60)) # 100.000

####### length 80 ############
queries_l80 = []
for record in iv.fasta.reader(file="../data/illumina_reads_80.fasta.gz"):
    queries_l80.append(str(record.seq))

# print(len(queries_l80)) # 100.000

n = 10**4 # so it doesn't take to long for the whole genome
# comment in / out the respective searches
res_40 = run_query_search_on_genome(whole_genome, queries_l40, n)
res_60 = run_query_search_on_genome(whole_genome, queries_l60, n)
res_80 = run_query_search_on_genome(whole_genome, queries_l80, n)
res_100 = run_query_search_on_genome(whole_genome, queries_100, n)