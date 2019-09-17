import math
import os

import keras
import numpy as np
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
from keras.layers import Input, Dense, BatchNormalization, Dropout, Conv1D, MaxPooling1D, Flatten, LSTM
from keras.layers.advanced_activations import LeakyReLU
from keras.models import Model, load_model
from keras.optimizers import Nadam
from keras.utils import to_categorical

from exoNet.models._losses import METRICS, LOSSES
from exoNet.models._network import Network
from exoNet.utils import label_encoder, train_test_split_data


class ExoCNNLSTM(Network):
    def __init__(self, seq_len, n_channels, n_classes=2, **kwargs):
        super().__init__()
        self.seq_len = seq_len
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.model_name = ExoCNNLSTM.__name__

        self.lr = kwargs.get("learning_rate", 0.001)
        self.dr_rate = kwargs.get("dropout_rate", 0.2)
        self.model_path = kwargs.get("model_path", "./models/FCN/")
        self.loss_fn = kwargs.get("loss_fn", 'cce')
        self.lambda_l1 = kwargs.get('lambda_l1', 0.0)
        self.lambda_l2 = kwargs.get('lambda_l2', 0.0)
        self.use_batchnorm = kwargs.get("use_batchnorm", False)

        self.x = Input(shape=(self.seq_len, self.n_channels,), name="data")

        self.network_kwargs = {
            "seq_len": self.seq_len,
            "n_channels": self.n_channels,
            "n_classes": self.n_classes,
            "dropout_rate": self.dr_rate,
            "loss_fn": self.loss_fn,
            "use_batchnorm": self.use_batchnorm,
        }

        self.training_kwargs = {
            "learning_rate": self.lr,
            "model_path": self.model_path,
            "lambda_l1": self.lambda_l1,
            "lambda_l2": self.lambda_l2,
        }

        self.init_w = keras.initializers.glorot_normal()
        self.regularizer = keras.regularizers.l1_l2(self.lambda_l1, self.lambda_l2)
        self._create_network()
        self._compile_models()

        print_summary = kwargs.get("print_summary", True)

        if print_summary:
            self.model.summary()

    def _create_network(self):
        h = Conv1D(filters=32, kernel_size=8, padding='same', kernel_initializer=self.init_w, use_bias=False,
                   kernel_regularizer=self.regularizer)(self.x)
        if self.use_batchnorm:
            h = BatchNormalization(trainable=True)(h)
        h = LeakyReLU()(h)
        h = MaxPooling1D(pool_size=2)(h)

        h = LSTM(16, activation='linear', return_sequences=False)(h)  # (None, 16)
        h = LeakyReLU()(h)

        h = Flatten()(h)

        h = Dense(32, kernel_initializer=self.init_w, kernel_regularizer=self.regularizer, use_bias=False)(h)
        if self.use_batchnorm:
            h = BatchNormalization()(h)
        h = LeakyReLU()(h)
        if self.dr_rate > 0:
            h = Dropout(self.dr_rate)(h)

        output = Dense(self.n_classes, kernel_initializer=self.init_w, kernel_regularizer=self.regularizer,
                       bias_initializer=keras.initializers.Constant(-math.log((1 - 0.01) / 0.01)),
                       activation='softmax')(h)

        self.model = Model(inputs=self.x, outputs=output)

    def _compile_models(self):
        self.optimizer = Nadam(lr=self.lr)
        self.model.compile(optimizer=self.optimizer, loss=LOSSES[self.loss_fn],
                           metrics=['acc', METRICS['sensitivity'], METRICS['specificity']])

    def to_latent(self):
        pass

    def predict(self, adata):
        adata = remove_sparsity(adata)

        return self.label_encoder.inverse_transform(np.argmax(self.model.predict(adata.X), axis=1))

    def save_model(self):
        self.model.save(os.path.join(self.model_path, f"{self.model_name}.h5"), overwrite=True)

    def restore_model(self):
        self.model = load_model(os.path.join(self.model_path, f"{self.model_name}.h5"), compile=False)
        self._compile_models()

    def train(self, adata, label_key, le=None, n_epochs=500, batch_size=32, early_stopping_kwargs={},
              lr_reducer_kwargs={}, verbose=2):
        adata = remove_sparsity(adata)

        train_adata, valid_adata = train_test_split_data(adata, 0.80)

        train_labels, self.label_encoder = label_encoder(train_adata, label_key=label_key, label_encoder=le)
        train_labels = to_categorical(train_labels, num_classes=self.n_classes)

        valid_labels, self.label_encoder = label_encoder(valid_adata, label_key=label_key, label_encoder=le)
        valid_labels = to_categorical(valid_labels, num_classes=self.n_classes)

        callbacks = []

        if early_stopping_kwargs != {}:
            callbacks.append(EarlyStopping(**early_stopping_kwargs))

        if lr_reducer_kwargs != {}:
            callbacks.append(ReduceLROnPlateau(**lr_reducer_kwargs))

        x_train = train_adata.X.reshape(shape=(-1, self.seq_len, self.n_channels))
        y_train = train_labels

        x_valid = valid_adata.X.reshape(shape=(-1, self.seq_len, self.n_channels))
        y_valid = valid_labels

        self.model.fit(x=x_train,
                       y=y_train,
                       validation_data=(x_valid, y_valid),
                       epochs=n_epochs,
                       batch_size=batch_size,
                       verbose=verbose,
                       callbacks=callbacks,
                       )
