
def build_model(num_classes):
    inputs = layers.Input(shape=model_input_shape)
    model = keras.applications.EfficientNetV2S(
        include_top=False,
        weights="imagenet",
        input_tensor=inputs,
        # classes=num_classes,
        # classifier_activation="sigmoid",
        include_preprocessing=True)

    # Freeze the pretrained weights
    model.trainable = False

    # Rebuild top
    x = layers.GlobalAveragePooling2D(name="avg_pool")(model.output)
    x = layers.BatchNormalization()(x)

    top_dropout_rate = 0.2
    x = layers.Dropout(top_dropout_rate, name="top_dropout")(x)
    outputs = layers.Dense(num_classes, activation="sigmoid", name="pred")(x)

    # Compile
    model = keras.Model(inputs, outputs, name="EfficientNet")
    optimizer = keras.optimizers.Adam(learning_rate=model_learning_rate)
    model.compile(
        optimizer=optimizer, loss="binary_crossentropy", metrics=model_metrics
    )
    return model


def unfreeze_model(model):
    # We unfreeze the top 20 layers while leaving BatchNorm layers frozen
    for layer in model.layers[-20:]:
        if not isinstance(layer, layers.BatchNormalization):
            layer.trainable = True

    optimizer = keras.optimizers.Adam(
        learning_rate=model_finetune_learning_rate)

    model.compile(
        optimizer=optimizer,
        loss="binary_crossentropy",
        metrics=model_metrics
    )
