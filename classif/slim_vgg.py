from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import glob
import os

import tensorflow as tf
import tensorflow.contrib.slim.nets as nets

slim = tf.contrib.slim
inception = nets.inception

def parser(serialized_example):
  """Parses a single tf.Example into image and label tensors."""
  features = tf.parse_single_example(
      serialized_example,
      features={
        "data": tf.FixedLenFeature((), tf.string, default_value=""),
        "label": tf.FixedLenFeature((), tf.int64, default_value=0)
      })
  image = tf.decode_raw(features['data'], tf.uint8)
  # TF does not suppport operations on uint8
  image = tf.cast(image, tf.float32)
  image = tf.reshape(image, [150, 2014, 4])

  #label = tf.cast(features["label"], tf.int64)
  #print(label)
  label = features["label"]
  #print(label)

  return image, label

def train_input_fn():
  # Import MNIST data
  batch_size = 128
  fnames = glob.glob(os.path.join("/mnt/disk3/felix/tensors/tfrecords", "*.tfrecords"))
  #train_filenames = "../data/tensors/NA12878.tfrecords"
  dataset = tf.data.TFRecordDataset(fnames)

  # Map the parser over dataset, and batch results by up to batch_size
  dataset = dataset.map(parser, num_parallel_calls=32)
  dataset = dataset.batch(batch_size)
  dataset = dataset.repeat()
  dataset = dataset.prefetch(buffer_size=10*batch_size)

  return dataset.make_one_shot_iterator().get_next()

train_log_dir = "/mnt/disk4/felix/log"

if not tf.gfile.Exists(train_log_dir):
  tf.gfile.MakeDirs(train_log_dir)

with tf.Graph().as_default():
  # Set up the data loading:
  images, labels = train_input_fn()
  labels = slim.one_hot_encoding(labels, 6)

  # Define the model:
  predictions, endpoints = inception.inception_v3(images, is_training=True, num_classes=6,
          spatial_squeeze=True)

  print(predictions.get_shape())
  # Specify the loss function:
  slim.losses.softmax_cross_entropy(endpoints["Logits"], labels)

  total_loss = slim.losses.get_total_loss()
  tf.summary.scalar('losses/total_loss', total_loss)

  # Specify the optimization scheme:
  optimizer = tf.train.GradientDescentOptimizer(learning_rate=.001)

  # create_train_op that ensures that when we evaluate it to get the loss,
  # the update_ops are done and the gradient updates are computed.
  train_tensor = slim.learning.create_train_op(total_loss, optimizer)

  # Actually runs training.
  slim.learning.train(train_tensor, train_log_dir)
