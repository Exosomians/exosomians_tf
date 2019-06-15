from data_loader import DataLoader
from keras.callbacks import EarlyStopping, ReduceLROnPlateau, CSVLogger
from network import Network


def main():
    print("Loading Data ...")
    data_loader = DataLoader(data_path="../Data/TertiaryData.csv", label_path="../Data/TertiaryLabelsMat.csv")
    x_train, y_train = data_loader.x_train, data_loader.y_train
    x_test, y_test = data_loader.x_test, data_loader.y_test
    print("Data has been loaded!")

    net = Network(input_shape=(None, 4),
                  model="Dense",
                  dropout_rate=0.25,
                  alpha=0.2,
                  beta=100.0
                  )
    reduce_lr = ReduceLROnPlateau(monitor='loss', factor=0.5,
                                  patience=20, min_lr=0.000001)
    csv_logger = CSVLogger(filename="./csv_logger.log")
    early_stopping = EarlyStopping(patience=25, monitor='val_loss')
    callbacks = [early_stopping, csv_logger]
    net.model.fit(
        x=x_train,
        y=y_train,
        batch_size=128,
        epochs=1000,
        verbose=2,
        callbacks=callbacks,
        shuffle=True,
        validation_data=(x_test, y_test)
    )


if __name__ == '__main__':
    main()
