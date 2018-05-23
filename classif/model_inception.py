import tensorflow as tf
import tensorflow.contrib.slim.nets as nets
import h5py

slim = tf.contrib.slim


def main():
    train = h5py.File("/mnt/disk4/felix/h5/node13.h5")
    eval = h5py.File("/mnt/disk4/felix/h5/node13_p.h5")

    train_data = train["data"]
    train_labels = train["labels"]
    eval_data = eval["data"]
    eval_labels = eval["labels"]

    g = tf.Graph()
    with g.as_default():
        model = nets.inception.inception_v3(num_classes=6)

        images, labels, _ = data_providers.make_batches(
            dataset.get_slim_dataset(), model, FLAGS.batch_size, mode='TRAIN')
        endpoints = model.create(images, dataset.num_classes, is_training=True)
        labels = slim.one_hot_encoding(labels, dataset.num_classes)
        total_loss = loss(
            endpoints['Logits'], labels, label_smoothing=FLAGS.label_smoothing)

        # Setup the moving averages:
        moving_average_variables = slim.get_model_variables()
        moving_average_variables.extend(slim.losses.get_losses())
        moving_average_variables.append(total_loss)

        variable_averages = tf.train.ExponentialMovingAverage(
            FLAGS.moving_average_decay, slim.get_or_create_global_step())

        tf.add_to_collection(tf.GraphKeys.UPDATE_OPS,
                             variable_averages.apply(moving_average_variables))

        # Configure the learning rate using an exponetial decay.
        decay_steps = int(((1.0 * dataset.num_examples) / FLAGS.batch_size) *
                          FLAGS.num_epochs_per_decay)

        learning_rate = tf.train.exponential_decay(
            FLAGS.learning_rate,
            slim.get_or_create_global_step(),
            decay_steps,
            FLAGS.learning_rate_decay_factor,
            staircase=True)

        opt = tf.train.RMSPropOptimizer(learning_rate, FLAGS.rmsprop_decay,
                                        FLAGS.rmsprop_momentum,
                                        FLAGS.rmsprop_epsilon)

        # Create training op
        train_tensor = slim.learning.create_train_op(
            total_loss,
            optimizer=opt,
            update_ops=tf.get_collection(tf.GraphKeys.UPDATE_OPS))

        # Summaries:
        slim.summaries.add_histogram_summaries(slim.get_model_variables())
        slim.summaries.add_scalar_summaries(slim.losses.get_losses(), 'losses')
        slim.summaries.add_scalar_summary(total_loss, 'Total_Loss', 'losses')
        slim.summaries.add_scalar_summary(learning_rate, 'Learning_Rate',
                                          'training')
        slim.summaries.add_histogram_summaries(endpoints.values())
        slim.summaries.add_zero_fraction_summaries(endpoints.values())
        # redacted

        # Set start-up delay
        startup_delay_steps = FLAGS.task * FLAGS.startup_delay_steps

        init_fn = model_init_function(model, dataset.num_classes,
                                      FLAGS.start_from_checkpoint)

        saver = tf.train.Saver(
            max_to_keep=FLAGS.max_checkpoints_to_keep,
            keep_checkpoint_every_n_hours=FLAGS.keep_checkpoint_every_n_hours)

        # Train model
        slim.learning.train(
            train_tensor,
            number_of_steps=FLAGS.number_of_steps,
            logdir=FLAGS.train_dir,
            master=target,
            init_fn=init_fn,
            is_chief=is_chief,
            saver=saver,
            startup_delay_steps=startup_delay_steps,
            save_summaries_secs=FLAGS.save_summaries_secs,
            save_interval_secs=FLAGS.save_interval_secs)
