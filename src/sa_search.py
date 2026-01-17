import numpy as np
import iv2py as iv
import time 

# helper functions for suffix array search
def find_LP(pattern, SA, text):
    m = len(pattern)
    # falls außerhalb der bounds
    if pattern <= text[SA[0] : SA[0] + m]:
        lp = 0
    elif pattern > text[SA[-1] : SA[-1] + m]:
        lp = len(SA)
    else:
        L, R = 0, len(SA)-1
        while R - L > 1:
            M = int(np.ceil((L + R) / 2))
            if pattern <= text[SA[M] : SA[M] + m]: # [0:m]
                R = M
            else:
                L = M
        lp = R
    return lp

def find_RP(pattern, SA, text):
    m = len(pattern)
    # falls außerhalb der bounds
    if pattern < text[SA[0] : SA[0] + m]:
        rp = 0
    elif pattern >= text[SA[-1] : SA[-1] + m]:
        rp = len(SA)
    else:
        L, R = 0, len(SA)-1
        while R - L > 1:
            M = int(np.ceil((L + R) / 2))
            if text[SA[M] : SA[M]+m] > pattern: # [0:m]
                R = M
            else:
                L = M
        rp = R
    return rp


def SA_search(text, pattern, sa):
    LP = find_LP(pattern, sa, text)
    RP = find_RP(pattern, sa, text)
    return [sa[i] for i in range(LP, RP)]

# file contains 100 000 queries, function to return as many as needed
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


def run_query_search(reference, n_queries, sa, n):
    """
    function to search for n specific queries in a given query list
    input:
    reference: text to be searched
    n_queries: list of queries
    sa: the suffix array of the reference
    n: number of queries from the query list to be searched
    output:
    list of lists with occurences of the respective queries
    """
    start = time.perf_counter()

    res = []
    for i in range(n):
        res.append(SA_search(reference, n_queries[i], sa))

    end = time.perf_counter()

    runtime = end - start
    print(f"runtime for {n} queries: {runtime}")
    print(f"first 50 matches: {res[:10]}")
    return res

# human genome
reference = ""
for record in iv.fasta.reader(file="../data/hg38_partial.fasta.gz"):
    reference = str(record.seq)

# queries length 100
queries_100 = []
for record in iv.fasta.reader(file="../data/illumina_reads_100.fasta.gz"):
    queries_100.append(str(record.seq))

############### excercise 4 #######################

# 10**3, 10**4, 10**5, 10**6
# n = 10**6
# n_queries = number_of_queries(queries_100, n)
# sa = iv.create_suffixarray(reference)

# res_100 = run_query_search(reference, n_queries, sa, n)

############### excercise 5 #######################

n = 10**3

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

sa = iv.create_suffixarray(reference)

# res_40 = run_query_search(reference, queries_l40, n)
# res_60 = run_query_search(reference, queries_l60, n)
res_80 = run_query_search(reference, queries_l80, n)