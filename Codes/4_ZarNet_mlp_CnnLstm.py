## 4_ZarNet_mlp_CnnLstm
## all features are included

import argparse

import numpy as np
import pandas as pd
from keras import backend as K
from keras.callbacks import EarlyStopping, CSVLogger
from keras.layers import Conv1D, Dense, Flatten, Input, MaxPooling1D, BatchNormalization, Dropout, LeakyReLU, \
    concatenate, LSTM
from keras.models import Model
from keras.optimizers import Adam, Nadam
from keras.preprocessing.sequence import pad_sequences
from keras.regularizers import l2, l1
from keras.utils import multi_gpu_model
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, normalize
from matplotlib import pyplot as plt


model_name = '4_ZarNet_mlp_CnnLstm'

parser = argparse.ArgumentParser(description='The higher order of ZarNet!!')
arguments_group = parser.add_argument_group("Parameters")
arguments_group.add_argument('-r', '--learning_rate', type=float, default=0.002, required=False,
                             help='Learning Rate')
arguments_group.add_argument('-l', '--lambda', type=float, default=1.0, required=False,
                             help='L1 and L2 coefficient')
arguments_group.add_argument('-e', '--epochs', type=int, default=500, required=False,
                             help='Number of epochs')
arguments_group.add_argument('-b', '--batch_size', type=int, default=128, required=False,
                             help='Batch size')
arguments_group.add_argument('-p', '--patience', type=int, default=100, required=False,
                             help='Patience of EarlyStopping')
arguments_group.add_argument('-d', '--dropout_rate', type=float, default=0.4, required=False,
                             help='Dropout rate')
arguments_group.add_argument('-w', '--weight', type=float, default=0.1, required=False,
                             help='class weight for Label 1')


args = vars(parser.parse_args())

testsize = 0.25
learning_rate = args['learning_rate']
lambda_value = args['lambda']
dropout_rate = args['dropout_rate']
patience = args['patience']
batch_size = args['batch_size']
epochs = args['epochs']
weight = args['weight']



data = pd.read_csv("../Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_element_string.csv")
data.drop(['Unnamed: 0', 'id', 'annotation', 'rnaType', 'ic', 'ev','element_string_number'], axis=1, inplace=True)
deepbind = pd.read_csv('../Data/deepbind_features.txt',sep='\t')
print(data.shape)
data = pd.concat([data, deepbind], axis=1)
print(data.shape)
rowsToRemove = [i for i in range(len(data['seq'])) if set(data['seq'][i]) == {'C', 'A', 'G', 'U', 'N'} ]
data.drop(rowsToRemove, axis=0, inplace=True)
data = data.reset_index()
data.drop('index', axis=1, inplace=True)


## preprocess cnn input data
string_features = ['seq', 'DB', 'element_string']
sequences = data['seq']
dotbracket = data['DB']
element = data['element_string']
labels = data['label']
labels = np.array(labels)
labels[labels == 'NO'] = 0
labels[labels == 'YES'] = 1

max_len = sequences.str.len().max()
print('max_len: ',max_len)



def preprocess_mlp(data):
    count_features = data.columns.tolist()[data.columns.tolist().index('AAAA'):data.columns.tolist().index('G111')]
    for column in count_features:
        data[column] /= data['length']

    count_features_df = data.drop(['strand', 'chr'], axis=1)
    count_features_df = pd.DataFrame(normalize(count_features_df.values, norm='max'), columns=count_features_df.columns)

    strand = data['strand']
    le = LabelEncoder()
    data['strand'] = pd.Series(le.fit_transform(strand.values))

    chromosome = data['chr']
    le = LabelEncoder()
    data['chr'] = pd.Series(le.fit_transform(chromosome.values))

    return pd.concat([data['chr'], data['strand'], count_features_df], axis=1)

data.drop(string_features, axis=1, inplace=True)
mlp_data = preprocess_mlp(data.drop(['label'], axis=1))



## encoding and mergeing the sequance features for CNN
def encoder(VectorOfStrings, charsToFit, max_len):
    ## integer encoding
    le = LabelEncoder()
    le.fit(list(charsToFit))
    integer_encoded = [le.transform(list(VectorOfStrings[i])).tolist() for i in range(VectorOfStrings.shape[0])]

    ## one-hot encoding
    oe = OneHotEncoder(categories='auto')
    oe = oe.fit(np.reshape(list(range(len(charsToFit))), (-1, 1)))
    onehot_encoded = [oe.transform(np.reshape(integer_encoded[i], (-1, 1))).A for i in range(len(integer_encoded))]
    onehot_encoded = np.array(onehot_encoded)

    def _add_pad(max_len, a_onehot, charsToFit):
        pad = np.zeros((max_len - len(a_onehot), len(charsToFit)))
        a_padded_onehot = np.concatenate((pad, a_onehot))
        return (a_padded_onehot)

    onehot_encoded_padded = [_add_pad(max_len, onehot_encoded[i], charsToFit) for i in range(len(onehot_encoded))]
    onehot_encoded_padded = np.array(onehot_encoded_padded)
    print('shape after padding', onehot_encoded_padded.shape)

    return (onehot_encoded_padded)


encoded_sequences = encoder(sequences, 'AUCG', max_len)
dotbracket = dotbracket.str.replace(')', '(', regex=False)
encoded_dotbracket = encoder(dotbracket, '.(', max_len)
encoded_element = encoder(element, 'ftsimh', max_len)
merged_sequences = np.concatenate((encoded_sequences, encoded_dotbracket, encoded_element), axis=2)
print(merged_sequences.shape)




def sensitivity(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    return true_positives / (possible_positives + K.epsilon())

def specificity(y_true, y_pred):
    true_negatives = K.sum(K.round(K.clip((1 - y_true) * (1 - y_pred), 0, 1)))
    possible_negatives = K.sum(K.round(K.clip(1 - y_true, 0, 1)))
    return true_negatives / (possible_negatives + K.epsilon())




def zarnet(seq_len, n_mlp_features, onehot_len):
    sequences = Input(shape=(seq_len, onehot_len,))
    mlp_input = Input(shape=(n_mlp_features,))

    def cnn_rnn(sequences):
        conv = Conv1D(filters=8, kernel_size=8, activation='linear', padding='same',
                      kernel_regularizer=l2(lambda_value))(sequences)
        conv = LeakyReLU()(conv)
        max_pool = MaxPooling1D(pool_size=2)(conv)

        conv = Conv1D(filters=16, kernel_size=5, activation='linear', padding='same',
                      kernel_regularizer=l2(lambda_value))(max_pool)
        conv = LeakyReLU()(conv)
        max_pool = MaxPooling1D(pool_size=2)(conv)

        conv = Conv1D(filters=16, kernel_size=3, activation='linear', padding='same',
                      kernel_regularizer=l2(lambda_value))(max_pool)
        conv = LeakyReLU()(conv)
        max_pool = MaxPooling1D(pool_size=2)(conv)

        conv = Conv1D(filters=16, kernel_size=3, activation='linear', padding='same',
                      kernel_regularizer=l2(lambda_value))(max_pool)
        conv = LeakyReLU()(conv)
        max_pool = MaxPooling1D(pool_size=2)(conv)

        # (None, 31, 32)
        lstm = LSTM(64, activation='linear', return_sequences=False)(max_pool)  # (None, 16)
        lstm = LeakyReLU()(lstm)

        return lstm

    def mlp(mlp_input):
        dense = Dense(128, activation='linear', kernel_regularizer=l1(lambda_value))(mlp_input)
        dense = BatchNormalization()(dense)
        dense = LeakyReLU()(dense)
        dense = Dropout(dropout_rate)(dense)
        return dense

    dense = concatenate([cnn_rnn(sequences), mlp(mlp_input)], axis=1)

    dense = Dense(128, activation='linear', kernel_regularizer=l2(lambda_value))(dense)
    dense = BatchNormalization()(dense)
    dense = LeakyReLU()(dense)
    dense = Dropout(dropout_rate)(dense)

    output = Dense(1, activation='sigmoid')(dense)

    model = Model(inputs=[sequences, mlp_input], outputs=output)
    model.compile(optimizer=Nadam(lr=learning_rate), loss='binary_crossentropy',
                  metrics=['acc', sensitivity, specificity])
    model.summary()

    return model

model = zarnet(max_len, mlp_data.shape[1], merged_sequences.shape[2])
#gpu_model = multi_gpu_model(model, gpus=4)
#gpu_model.compile(optimizer=Nadam(lr=learning_rate), loss='binary_crossentropy',
#                  metrics=['acc', sensitivity, specificity])

# Train/Test split
x1_train, x1_test, x2_train, x2_test, y_train, y_test = train_test_split(merged_sequences, mlp_data, labels,
                                                                         test_size=testsize, stratify=data[["label"]])
print(x1_train.shape, x2_train.shape, y_train.shape)
print(x1_test.shape, x2_test.shape, y_test.shape)



early_stopping = EarlyStopping(patience=patience, monitor='val_loss', mode='min')
csv_logger = CSVLogger(filename="./"+model_name+"_train.log")


model.fit(x=[x1_train, x2_train],
          y=y_train,
          validation_data=([x1_test, x2_test], y_test),
          epochs=epochs,
          batch_size=batch_size,
          class_weight={0: 10.0 - weight, 1: weight},
          verbose=2,
          callbacks=[early_stopping, csv_logger],
          shuffle=True)

model.save("./"+model_name+".h5")



def visualize_results():
    history = pd.read_csv("./" + model_name + "_train.log")
    plt.close("all")
    plt.plot(history['acc'])
    plt.plot(history['val_acc'])
    plt.title('Model Accuracy')
    plt.ylabel('Accuracy')
    plt.ylim((0.0, 1.0))
    plt.xlabel('Epoch')
    plt.legend(['Training', 'Validation'], loc='best')
    plt.savefig("./" + model_name + "_Accuracy.pdf", dpi=300)
    plt.close()

    plt.plot(history['sensitivity'])
    plt.plot(history['val_sensitivity'])
    plt.title('Model Sensitivity')
    plt.ylabel('Sensitivity')
    plt.ylim((0.0, 1.0))
    plt.xlabel('Epoch')
    plt.legend(['Training', 'Validation'], loc='best')
    plt.savefig("./" + model_name + "_Sensitivity.pdf", dpi=300)
    plt.close()

    plt.plot(history['specificity'])
    plt.plot(history['val_specificity'])
    plt.title('model Specificity')
    plt.ylabel('Specificity')
    plt.ylim((0.0, 1.0))
    plt.xlabel('Epoch')
    plt.legend(['Training', 'Validation'], loc='best')
    plt.savefig("./" + model_name + "_Specificity.pdf", dpi=300)
    plt.close()

    plt.plot(history['loss'])
    plt.plot(history['val_loss'])
    plt.title('Model Loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['Training', 'Validation'], loc='best')
    plt.savefig("./" + model_name + "_Loss.pdf", dpi=300)
    plt.close()


visualize_results()
