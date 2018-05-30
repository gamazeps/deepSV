from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import glob
import os

import numpy as np
import tensorflow as tf

tf.logging.set_verbosity(tf.logging.INFO)

def sv_model_fn(features, labels, mode):
  """Model function for CNN."""
  # Input Layer
  # Reshape X to 4-D tensor: [batch_size, height, width, channels]
  #input_layer = tf.reshape(features, [-1, 150, 2014, 4])

  # Flatten tensor into a batch of vectors
  flat = tf.reshape(features, [-1, 150 * 2014 * 4])

  # Logits layer
  # Input Tensor Shape: [batch_size, 1024]
  # Output Tensor Shape: [batch_size, 10]
  logits = tf.layers.dense(inputs=flat, units=6)

  predictions = {
      # Generate predictions (for PREDICT and EVAL mode)
      "classes": tf.argmax(input=logits, axis=1),
      # Add `softmax_tensor` to the graph. It is used for PREDICT and by the
      # `logging_hook`.
      "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
  }
  if mode == tf.estimator.ModeKeys.PREDICT:
    return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)

  # Calculate Loss (for both TRAIN and EVAL modes)
  loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)

  # Configure the Training Op (for TRAIN mode)
  if mode == tf.estimator.ModeKeys.TRAIN:
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
    train_op = optimizer.minimize(
        loss=loss,
        global_step=tf.train.get_global_step())
    return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)

  # Add evaluation metrics (for EVAL mode)
  eval_metric_ops = {
      "accuracy": tf.metrics.accuracy(
          labels=labels, predictions=predictions["classes"])}
  return tf.estimator.EstimatorSpec(
      mode=mode, loss=loss, eval_metric_ops=eval_metric_ops)

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

  return image, features['label']

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

  return dataset

def main(unused_argv):
  # Load training and eval data

  # Create the Estimator
  sv_classifier = tf.estimator.Estimator(
      model_fn=sv_model_fn, model_dir="/tmp/sv_convnet_model")

  # Set up logging for predictions
  # Log the values in the "Softmax" tensor with label "probabilities"
  tensors_to_log = {"probabilities": "softmax_tensor"}
  logging_hook = tf.train.LoggingTensorHook(
      tensors=tensors_to_log, every_n_iter=50)

  # Train the model
  sv_classifier.train(
      input_fn=train_input_fn,
      steps=10000,
      hooks=[logging_hook, tf.train.ProfilerHook(show_memory=True, save_steps=100)])

  # Evaluate the model and print results
  eval_results = mnist_classifier.evaluate(input_fn=train_input_fn, steps=100)
  print(eval_results)


if __name__ == "__main__":
  tf.app.run()

