import iv2py as iv
import time

reference = ""
for record in iv.fasta.reader(file="../data/hg38_partial.fasta.gz"):
    reference = str(record.seq)

queries_100 = []
for record in iv.fasta.reader(file="../data/illumina_reads_100.fasta.gz"):
    queries_100.append(str(record.seq))

def naiveSearch(text, pattern):
    if len(pattern) == 0 or len(text) < len(pattern):
        return
    res = []
    m = len(pattern)
    for i in range(len(text)): # -m+1
        if text[i:i+m] == pattern:
            res.append(i)
    return res

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

# preparation
# 10**3, 10**4, 10**5, 10**6
n = 10**6
n_queries = number_of_queries(queries_100, n)

# time the actual search
start = time.perf_counter()

res = []
for i in range(n):
    res.append(naiveSearch(reference, n_queries[i]))
# print(len(querys_100)) # 100 000 querys
# print(res)

end = time.perf_counter()

runtime = end - start
print(f"runtime for {n} queries: {runtime:.2f}")
print(f"first 50 matches: {res[:10]}")
