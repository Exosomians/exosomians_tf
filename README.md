# Exosomians

## Introduction 
A Keras (with tensorflow backend) implementation of exoNet. 
<div float="left">
  <img src="https://www.tensorflow.org/images/tf_logo_transp.png" height="80" >
  <img src="https://s3.amazonaws.com/keras.io/img/keras-logo-2018-large-1200.png" height="80">
</div>
<div float="right">
</div>

## Getting Started

## Installation

### Installation with pip
To install the latest version from PyPI, simply use the following bash script:
```bash
pip install exonet
```
or install the development version via pip: 
```bash
pip install git+https://github.com/saberi1/Exosomians.git
```

or you can clone this repository and install via setup.py file:
```bash
git clone https://github.com/saberi1/Exosomians
cd Exosomians
python setup.py -q
``` 

## Examples

### Preprocessing

#### Sequence

You can encoded raw sequences with a character encoder simply by doing the following:

```python
import exoNet
import numpy as np

# Read raw sequences 
rna_sequences = np.loadtxt("./data/rna_sequnces.csv", delimiter=",")

char_encoder = {
    'N': -1, # exoNet Can create a vector of [0.25, 0.25, 0.25, 0.25]
    'A': 0,
    'C': 1,
    'G': 2,
    'U': 3,
}

# One-hot encoding
encoded_sequences = exoNet.pp.seq_encoder(seq_data=rna_sequences,
                                          char_encoder=char_encoder,
                                          unknown_char=True)
                                          
```


### ExoFCN
```python
import exoNet
import numpy as np

# Reading pre-processed data
fcn_data = np.load("./data/features.npy")
labels = np.load("./data/labels.npy")

# Create network with arbitrary hyper-parameters
network = exoNet.models.ExoFCN(x_dimension=fcn_data.shape[1], 
                               n_classes=len(np.unique(labels).tolist()),
                               architecture=[128, 32],
                               use_batchnorm=True,
                               dropout_rate=0.2,
                               model_path="./models/ExoFCN/",
                               lr=0.0001,
)

# Training network
network.train(data=fcn_data,
              labels=labels,
              le={"NO": 0, "YES": 1},
              n_epochs=300,
              batch_size=32,
              early_stopping_kwargs={"patience": 20, "monitor": "val_loss"},
              verbose=2,
)
# Get last hidden layer activations output (For visualzation purposes)
latent = network.to_latent(fcn_data)

# Calculate accuracy on validation data
valid_fcn_data = np.load("./data/valid_features.npy")
valid_labels = np.load("./data/valid_labels.npy")
pred_labels = network.predict(valid_fcn_data)
print(np.sum(pred_labels == valid_labels) / (valid_labels.shape[0]))
```

### ExoCNN
```python
import exoNet
import numpy as np

# Reading pre-processed data
seq_data = np.load("./data/sequences.npy")
labels = np.load("./data/labels.npy")

# Create network with arbitrary hyper-parameters
network = exoNet.models.ExoCNN(seq_len=seq_data.shape[1],
                               n_channels=seq_data.shape[2], 
                               n_classes=len(np.unique(labels).tolist()),
                               use_batchnorm=True,
                               dropout_rate=0.2,
                               model_path="./models/ExoFCN/",
                               lr=0.0001,
)

# Training network
network.train(seq_data=seq_data,
              labels=labels,
              le={"NO": 0, "YES": 1},
              n_epochs=300,
              batch_size=32,
              early_stopping_kwargs={"patience": 20, "monitor": "val_loss"},
              verbose=2,
)
# Get last hidden layer activations output (For visualzation purposes)
latent = network.to_latent(seq_data)
exoNet.pl.plot_umap(latent, labels)

# Calculate accuracy on validation data
valid_seq_data = np.load("./data/valid_sequences.npy")
valid_labels = np.load("./data/valid_labels.npy")
pred_labels = network.predict(valid_seq_data)
print(np.sum(pred_labels == valid_labels) / (valid_labels.shape[0]))
```

### ExoCNNLSTM
```python
import exoNet
import numpy as np

# Reading pre-processed data
seq_data = np.load("./data/sequences.npy")
labels = np.load("./data/labels.npy")

# Create network with arbitrary hyper-parameters
network = exoNet.models.ExoCNNLSTM(seq_len=seq_data.shape[1],
                                   n_channels=seq_data.shape[2], 
                                   n_classes=len(np.unique(labels).tolist()),
                                   use_batchnorm=True,
                                   dropout_rate=0.2,
                                   model_path="./models/ExoFCN/",
                                   lr=0.0001,
)

# Training network
network.train(seq_data=seq_data,
              labels=labels,
              le={"NO": 0, "YES": 1},
              n_epochs=300,
              batch_size=32,
              early_stopping_kwargs={"patience": 20, "monitor": "val_loss"},
              verbose=2,
)
# Get last hidden layer activations output (For visualzation purposes)
latent = network.to_latent(seq_data)
exoNet.pl.plot_umap(latent, labels)

# Calculate accuracy on validation data
valid_seq_data = np.load("./data/valid_sequences.npy")
valid_labels = np.load("./data/valid_labels.npy")
pred_labels = network.predict(valid_seq_data)
print(np.sum(pred_labels == valid_labels) / (valid_labels.shape[0]))
```

### ExoLSTM
```python
import exoNet
import numpy as np

# Reading pre-processed data
seq_data = np.load("./data/sequences.npy")
labels = np.load("./data/labels.npy")

# Create network with arbitrary hyper-parameters
network = exoNet.models.ExoCNNLSTM(seq_len=seq_data.shape[1],
                                   n_channels=seq_data.shape[2], 
                                   n_classes=len(np.unique(labels).tolist()),
                                   use_batchnorm=True,
                                   dropout_rate=0.2,
                                   model_path="./models/ExoFCN/",
                                   lr=0.0001,
)

# Training network
network.train(seq_data=seq_data,
              labels=labels,
              le={"NO": 0, "YES": 1},
              n_epochs=300,
              batch_size=32,
              early_stopping_kwargs={"patience": 20, "monitor": "val_loss"},
              verbose=2,
)
# Get last hidden layer activations output (For visualzation purposes)
latent = network.to_latent(seq_data)
exoNet.pl.plot_umap(latent, labels)

# Calculate accuracy on validation data
valid_seq_data = np.load("./data/valid_sequences.npy")
valid_labels = np.load("./data/valid_labels.npy")
pred_labels = network.predict(valid_seq_data)
print(np.sum(pred_labels == valid_labels) / (valid_labels.shape[0]))
```

### ExoFCNN
```python
import exoNet
import numpy as np

# Reading pre-processed data
fcn_data = np.load("./data/features.npy")
seq_data = np.load("./data/sequences.npy")
labels = np.load("./data/labels.npy")

# Create network with arbitrary hyper-parameters
network = exoNet.models.ExoFCNN(n_features=fcn_data.shape[1],
                                seq_len=seq_data.shape[1],
                                n_channels=seq_data.shape[2], 
                                n_classes=len(np.unique(labels).tolist()),
                                use_batchnorm=True,
                                dropout_rate=0.2,
                                model_path="./models/ExoFCN/",
                                lr=0.0001,
)

# Training network
network.train(seq_data=seq_data,
              fcn_data=fcn_data,
              labels=labels,
              le={"NO": 0, "YES": 1},
              n_epochs=300,
              batch_size=32,
              early_stopping_kwargs={"patience": 20, "monitor": "val_loss"},
              verbose=2,
)
# Get last hidden layer activations output (For visualzation purposes)
latent = network.to_latent(seq_data)
exoNet.pl.plot_umap(latent, labels)

# Calculate accuracy on validation data
valid_seq_data = np.load("./data/valid_sequences.npy")
valid_fcn_data = np.load("./data/valid_features.npy")
valid_labels = np.load("./data/valid_labels.npy")
pred_labels = network.predict(valid_seq_data)
print(np.sum(pred_labels == valid_labels) / (valid_labels.shape[0]))
```


### ExoNet
```python
import exoNet
import numpy as np

# Reading pre-processed data
fcn_data = np.load("./data/features.npy")
seq_data = np.load("./data/sequences.npy")
labels = np.load("./data/labels.npy")

# Create network with arbitrary hyper-parameters
network = exoNet.models.ExoNet(n_features=fcn_data.shape[1],
                               seq_len=seq_data.shape[1],
                               n_channels=seq_data.shape[2], 
                               n_classes=len(np.unique(labels).tolist()),
                               use_batchnorm=True,
                               dropout_rate=0.2,
                               model_path="./models/ExoFCN/",
                               lr=0.0001,
)

# Training network
network.train(seq_data=seq_data,
              fcn_data=fcn_data,
              labels=labels,
              le={"NO": 0, "YES": 1},
              n_epochs=300,
              batch_size=32,
              early_stopping_kwargs={"patience": 20, "monitor": "val_loss"},
              verbose=2,
)
# Get last hidden layer activations output (For visualzation purposes)
latent = network.to_latent(seq_data)
exoNet.pl.plot_umap(latent, labels)

# Calculate accuracy on validation data
valid_seq_data = np.load("./data/valid_sequences.npy")
valid_fcn_data = np.load("./data/valid_features.npy")
valid_labels = np.load("./data/valid_labels.npy")
pred_labels = network.predict(valid_seq_data)
print(np.sum(pred_labels == valid_labels) / (valid_labels.shape[0]))
```
