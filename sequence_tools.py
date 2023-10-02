from source import fastq_read as fr


# Function for filtering FASTQ file. 
def fastq_filter(seqs: dict, gc_bound: tuple = (0, 100), length_bound: tuple = (0, 2**32), quality_threshold: float = 0) -> dict:
    dict = seqs.copy()
    for value in seqs:
        fastq = seqs[value][0]
        if type(gc_bound) == int or type(gc_bound) == float:
            gc_check = fr.gc_content(fastq) >= gc_bound
        else: 
            gc_check = fr.gc_content(fastq) < gc_bound[0] or fr.gc_content(fastq) > gc_bound[1]
        len_check = fr.seq_length(fastq) < length_bound[0] or fr.seq_length(fastq) > length_bound[1]
        quality_check = fr.quality_score(fastq) < quality_threshold
        if gc_check or len_check or quality_check:
            del dict[value]      
    return dict
