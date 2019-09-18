import numpy as np
from keras_preprocessing.sequence import pad_sequences


def seq_encoder(seq_data, char_encoder, max_len, unknown_char=True):
    if unknown_char:
        one_hot_encoder = {
            -1: [0.25 for i in range(len(char_encoder.keys()) - 1)],
            0: [0 for i in range(0)] + [1] + [0 for i in range(len(char_encoder.keys()) - 2 - 0)],
            1: [0 for i in range(1)] + [1] + [0 for i in range(len(char_encoder.keys()) - 2 - 1)],
            2: [0 for i in range(2)] + [1] + [0 for i in range(len(char_encoder.keys()) - 2 - 2)],
            3: [0 for i in range(3)] + [1] + [0 for i in range(len(char_encoder.keys()) - 2 - 3)],
        }
    else:
        one_hot_encoder = {
            -1: [0.25 for i in range(len(char_encoder.keys()))],
            0: [0 for i in range(0)] + [1] + [0 for i in range(len(char_encoder.keys()) - 1 - 0)],
            1: [0 for i in range(1)] + [1] + [0 for i in range(len(char_encoder.keys()) - 1 - 1)],
            2: [0 for i in range(2)] + [1] + [0 for i in range(len(char_encoder.keys()) - 1 - 2)],
            3: [0 for i in range(3)] + [1] + [0 for i in range(len(char_encoder.keys()) - 1 - 3)],
        }
    encoded_sequences = []
    for sequence in seq_data:
        encoded_sequence = [char_encoder[char] for char in sequence]
        encoded_sequences.append(encoded_sequence)

    encoded_sequences = pad_sequences(encoded_sequences, maxlen=max_len, padding='post', truncating='post', value=-1)
    onehot_sequences = []
    for encoded_sequence in encoded_sequences.tolist():
        onehot_sequence = [one_hot_encoder[enc] for enc in encoded_sequence]
        onehot_sequences.append(onehot_sequence)

    onehot_sequences = np.array(onehot_sequences)
    return onehot_sequences
