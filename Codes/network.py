import tensorflow as tf
from keras.layers import Input, Dense, BatchNormalization, Dropout, LeakyReLU
from keras.models import Model
from keras.optimizers import Adam
from utils import sensitivity, specificity, weighted_cross_entropy


# Dialating CNN, CNN -> Capture Structure of RNA -> Saliency Map, Attention
# Instead of K-mer -> Use DeepBind : Which RBP is attached

class Network(object):
    def __init__(self, input_shape, model="LSTM-Dense", **kwargs):
        tf.reset_default_graph()
        self.input_shape = input_shape
        self.dropout_rate = kwargs.get("dropout_rate", 0.25)
        self.learning_rate = kwargs.get("learning_rate", 0.001)
        self.alpha = kwargs.get("alpha", 0.2)
        self.beta = kwargs.get("beta", 100.0)

        self._create_model()
        self._compile_model()

    def _create_model(self):
        x = Input(shape=self.input_shape)

        dense = Dense(128, activation='relu')(x)
        dense = BatchNormalization()(dense)
        dense = LeakyReLU()(dense)
        dense = Dropout(self.dropout_rate)(dense)

        dense = Dense(16)(dense)
        dense = BatchNormalization()(dense)
        dense = LeakyReLU()(dense)
        dense = Dropout(self.dropout_rate)(dense)

        output = Dense(1, activation='sigmoid')(dense)

        self.model = Model(inputs=x, outputs=output)
        self.model.summary()

    def _compile_model(self):
        optimizer = Adam(self.learning_rate)
        class_weights = {"0": 1.0 - self.alpha, "1": self.alpha}
        loss = lambda y_true, y_pred: self.beta * weighted_cross_entropy(class_weights)(y_true, y_pred)
        metrics = ['accuracy', sensitivity, specificity]
        self.model.compile(optimizer=optimizer, loss=loss, metrics=metrics)
