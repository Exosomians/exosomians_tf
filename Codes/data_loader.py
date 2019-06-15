import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from keras.preprocessing.sequence import pad_sequences

class DataLoader(object):
    def __init__(self, data_path=None, label_path=None, shuffle=True, skip_cols=None):
        self.data_path = data_path
        self.label_path = label_path
        self.data = pd.read_csv(data_path, index_col=0)
        self.labels = pd.read_csv(label_path, index_col=0)["label"]
        self.shuffle = shuffle
        self._one_hot_encode()
        self.train_test_split(test_size=0.25)

    def train_test_split(self, test_size=0.25):
        x_train, x_test, y_train, y_test = train_test_split(self.seq_data, self.labels, test_size=test_size,
                                                            shuffle=self.shuffle, stratify=True)
        self.x_train, self.y_train = x_train, y_train
        self.x_test, self.y_test = x_test, y_test

    def _one_hot_encode(self):
        label_encoder = LabelEncoder()
        self.labels = label_encoder.fit_transform(np.reshape(self.labels.values, (-1, 1)))

        sequences = self.data['seq'].values
        label_encoder = LabelEncoder()
        label_encoder = label_encoder.fit(list('ATCGN'))
        integer_encoded_seqs = [label_encoder.transform(list(sequences[i])).tolist() for i in
                                range(sequences.shape[0])]
        integer_encoded_seqs = pad_sequences(integer_encoded_seqs, maxlen=50)
        # integer_encoded_seqs = np.array(integer_encoded_seqs)
        print(integer_encoded_seqs.shape)

        onehot_encoder = OneHotEncoder(sparse=False)
        onehot_encoder = onehot_encoder.fit(np.reshape(integer_encoded_seqs[0], (-1, 1)))
        onehot_encoded_seq = [onehot_encoder.transform(np.reshape(integer_encoded_seqs[i], (-1, 1))).tolist() for i in
                              range(len(integer_encoded_seqs))]
        onehot_encoded_seq = np.array(onehot_encoded_seq)

        self.seq_data = onehot_encoded_seq
        # print(self.seq_data)
