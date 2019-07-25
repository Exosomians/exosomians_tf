## 6_exoMLP
## includes only the independent features



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
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, normalize
from matplotlib import pyplot as plt


model_name = '6_exoMlp_subset'
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
arguments_group.add_argument('-w', '--weight', type=float, default=9.0, required=False,
                             help='class weight for Label 1')


args = vars(parser.parse_args())

testsize = 0.25
subsample_size = 5000
learning_rate = args['learning_rate']
lambda_value = args['lambda']
dropout_rate = args['dropout_rate']
patience = args['patience']
batch_size = args['batch_size']
epochs = args['epochs']
weight = args['weight']



## loading data and initial preprocessing 
data = pd.read_csv("../Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_element_string.csv")
deepbind = pd.read_csv('../Data/deepbind_features.txt',sep='\t')
print(data.shape)

data = pd.concat([data, deepbind], axis=1)
rowsToRemove = [i for i in range(len(data['seq'])) if set(data['seq'][i]) == {'C', 'A', 'G', 'U', 'N'} ]
data.drop(rowsToRemove, axis=0, inplace=True)
data.drop(['Unnamed: 0', 'id', 'annotation', 'rnaType', 'ic', 'ev', 'seq', 'DB', 'element_string','element_string_number'], axis=1, inplace=True)
data = data.reset_index()
data.drop('index', axis=1, inplace=True)
print(data.shape)



def preprocess_mlp(data):
    count_features = data.columns.tolist()[data.columns.tolist().index('AAAA'):]
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

mlp_data = preprocess_mlp(data.drop(['label'], axis=1))



## prepare the labels

labels = data['label']
labels = np.array(labels)
labels[labels == 'NO'] = 0
labels[labels == 'YES'] = 1


## defining performance measures

def sensitivity(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    return true_positives / (possible_positives + K.epsilon())

def specificity(y_true, y_pred):
    true_negatives = K.sum(K.round(K.clip((1 - y_true) * (1 - y_pred), 0, 1)))
    possible_negatives = K.sum(K.round(K.clip(1 - y_true, 0, 1)))
    return true_negatives / (possible_negatives + K.epsilon())




def exo_mlp(num_features):
    mlp_input = Input(shape=(num_features,))

    dense = Dense(128, activation='linear', kernel_regularizer=l2(lambda_value))(mlp_input)
    dense = BatchNormalization()(dense)
    dense = LeakyReLU()(dense)
    dense = Dropout(dropout_rate)(dense)

    dense = Dense(32, activation='linear', kernel_regularizer=l2(lambda_value))(dense)
    dense = BatchNormalization()(dense)
    dense = LeakyReLU()(dense)
    dense = Dropout(dropout_rate)(dense)

    output = Dense(1, activation='sigmoid')(dense)

    model = Model(inputs=mlp_input, outputs=output)
    model.compile(optimizer=Nadam(lr=learning_rate), loss='binary_crossentropy',
                  metrics=['acc', sensitivity, specificity])
    model.summary()

    return model


model = exo_mlp(mlp_data.shape[1])

# Subsampling
mlp_data_no = mlp_data.values[labels == 0]
mlp_data_yes = mlp_data.values[labels == 1]

labels_no = labels[labels == 0]
labels_yes = labels[labels == 1]

random_indices = np.random.choice(mlp_data_no.shape[0], subsample_size, replace=True)
mlp_data_no = mlp_data_no[random_indices]
labels_no = labels_no[random_indices]

mlp_data = np.concatenate([mlp_data_no, mlp_data_yes], axis=0)
labels = np.concatenate([labels_no, labels_yes], axis=0)

# Train/Test Split
x_train, x_test, y_train, y_test = train_test_split(mlp_data, labels, test_size=testsize, stratify=labels)
print(x_train.shape, y_train.shape)
print(x_test.shape, y_test.shape)

early_stopping = EarlyStopping(patience=patience, monitor='val_loss', mode='min')
csv_logger = CSVLogger(filename="./"+model_name+"_train.log")

model.fit(x=x_train,
          y=y_train,
          validation_data=(x_test, y_test),
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
