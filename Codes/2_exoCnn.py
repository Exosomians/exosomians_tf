## 2_exoCnn
## includes only the independent features


import argparse

import numpy as np
import pandas as pd
from keras import backend as K
from keras.callbacks import EarlyStopping, CSVLogger
from keras.layers import Conv1D, Dense, Flatten, Input, MaxPooling1D, BatchNormalization, Dropout, LeakyReLU
from keras.models import Model
from keras.optimizers import Adam
from keras.regularizers import l2
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

model_name = '2_exoCnn'
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

## loading data and initial preprocessing
data = pd.read_csv("../Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_element_string.csv")
rowsToRemove = [i for i in range(len(data['seq'])) if set(data['seq'][i]) == {'C', 'A', 'G', 'U', 'N'}]
data.drop(rowsToRemove, axis=0, inplace=True)
data = data.reset_index()
data.drop('index', axis=1, inplace=True)

## preprocess cnn input data
sequences = data['seq']
dotbracket = data['DB']
element = data['element_string']
max_len = sequences.str.len().max()
labels = data['label']
labels = np.array(labels)
labels[labels == 'NO'] = 0
labels[labels == 'YES'] = 1

print('max_len: ', max_len)


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


## defining performance measures
def sensitivity(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    return true_positives / (possible_positives + K.epsilon())


def specificity(y_true, y_pred):
    true_negatives = K.sum(K.round(K.clip((1 - y_true) * (1 - y_pred), 0, 1)))
    possible_negatives = K.sum(K.round(K.clip(1 - y_true, 0, 1)))
    return true_negatives / (possible_negatives + K.epsilon())


def cnn(seq_len, onehot_len):
    cnn_input = Input(shape=(seq_len, onehot_len,))

    conv = Conv1D(filters=8, kernel_size=8, padding='same', kernel_regularizer=l2(lambda_value))(cnn_input)
    conv = BatchNormalization()(conv)
    conv = LeakyReLU()(conv)
    max_pool = MaxPooling1D(pool_size=2)(conv)

    conv = Conv1D(filters=16, kernel_size=5, padding='same', kernel_regularizer=l2(lambda_value))(max_pool)
    conv = BatchNormalization()(conv)
    conv = LeakyReLU()(conv)
    max_pool = MaxPooling1D(pool_size=2)(conv)

    conv = Conv1D(filters=16, kernel_size=3, padding='same', kernel_regularizer=l2(lambda_value))(max_pool)
    conv = BatchNormalization()(conv)
    conv = LeakyReLU()(conv)
    max_pool = MaxPooling1D(pool_size=2)(conv)

    conv = Conv1D(filters=16, kernel_size=3, padding='same', kernel_regularizer=l2(lambda_value))(max_pool)
    conv = BatchNormalization()(conv)
    conv = LeakyReLU()(conv)
    max_pool = MaxPooling1D(pool_size=2)(conv)

    flat = Flatten()(max_pool)

    dense = Dense(128, kernel_regularizer=l2(lambda_value))(flat)
    dense = BatchNormalization()(dense)
    dense = LeakyReLU()(dense)
    dense = Dropout(dropout_rate)(dense)

    dense = Dense(32, kernel_regularizer=l2(lambda_value))(dense)
    dense = BatchNormalization()(dense)
    dense = LeakyReLU()(dense)
    dense = Dropout(dropout_rate)(dense)

    output = Dense(1, activation='sigmoid')(dense)

    model = Model(inputs=cnn_input, outputs=output)
    model.compile(optimizer=Adam(lr=learning_rate), loss='binary_crossentropy',
                  metrics=['acc', sensitivity, specificity])
    model.summary()
    return model


model = cnn(max_len, merged_sequences.shape[2])

# Train/Test Split
x_train, x_test, y_train, y_test = train_test_split(merged_sequences, labels, test_size=testsize,
                                                    stratify=data[["label"]])
print(x_train.shape, y_train.shape)
print(x_test.shape, y_test.shape)

early_stopping = EarlyStopping(patience=patience, monitor='val_loss', mode='min')
csv_logger = CSVLogger(filename="./" + model_name + "_train.log")

model.fit(x=x_train,
          y=y_train,
          validation_data=(x_test, y_test),
          epochs=epochs,
          batch_size=batch_size,
          class_weight={0: 10.0 - weight, 1: weight},
          verbose=2,
          callbacks=[early_stopping, csv_logger],
          shuffle=True)

model.save("./" + model_name + ".h5")


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
