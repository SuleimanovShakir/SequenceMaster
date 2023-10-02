from source import fastq_red as fr

# Поработать с условиями в 56 строке - мне кажется, что можно как-то упростить
# Поставить значения по умолчанию - обязательно нужно сделать - непонятно как они стакаются с аннотациями
def fastq_filter(seqs: dict, gc_bound: tuple, length_bound: tuple, quality_threshold: float) -> dict:
    dict = seqs.copy()
    for value in seqs:
        fastq = seqs[value][0]
        if (fr.gc_content(fastq) < gc_bound[0] or fr.gc_content(fastq) > gc_bound[1]) or (fr.seq_length(fastq) < length_bound[0] or fr.seq_length(fastq) > length_bound[1]) or (fr.quality_score(fastq) < quality_threshold):
            del dict[value]      
    return dict