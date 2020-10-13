from random import randint

def generate_barcoded_sequences(barcode_length, upstream_seq, downstream_seq, num_sequences):
    bases="ACGT"
    sequences=[]
    for _ in range(num_sequences):
        rand_barcode=""
        for _ in range(barcode_length):
            rand_base = bases[randint(0,3)]
            rand_barcode += rand_base
        sequences.append(upstream_seq + rand_barcode + downstream_seq)
    return sequences

def generate_random_sequences(sequence_length, num_sequences):
    bases="ACGT"
    sequences=[]
    for _ in range(num_sequences):
        rand_seq=""
        for _ in range(sequence_length):
            rand_base = bases[randint(0,3)]
            rand_seq += rand_base
        sequences.append(rand_seq)
    return sequences
