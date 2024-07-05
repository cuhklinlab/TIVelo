## TIVelo
A Tool to be Published on Nature

## Usage && Installation
Please follow the [Tutorials](https://tivelo.readthedocs.io/en/latest/) for installation and Usage.



+ ImportError: keras.optimizers.legacy is not supported in Keras 3. When using tf.keras, to continue using a tf.keras.optimizers.legacy optimizer, you can install the tf_keras package (Keras 2) and set the environment variable TF_USE_LEGACY_KERAS=True to configure TensorFlow to use tf_keras when accessing tf.keras.
```bash
pip install tf_keras
```
```python
import os
os.environ['TF_USE_LEGACY_KERAS'] = 'True'
```