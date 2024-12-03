
model = timm.create_model(model_name=model_name,
                          pretrained=True,
                          in_chans=n_input_channel,
                          num_classes=1)

# count_parameters(model)
net = NeuralNetBinaryClassifier(
    model,
    max_epochs=model_epochs,
    lr=model_learning_rate,
    batch_size=model_batch_size,
    criterion=nn.BCEWithLogitsLoss,
    optimizer=optim.Adam,
    device='cpu' if not torch.cuda.is_available() else 'cuda',
    # freeze bottom layers
    callbacks=[
        Freezer(lambda x: not x.startswith('classifier')),
        # EpochScoring("f1"), SacredLogger(,log_on_batch_end=True),
        LRScheduler(policy='StepLR', step_size=7, gamma=0.1)
    ],
    # shuffle training data
    iterator_train__shuffle=True,
    # class weights
    # criterion__weight=torch.Tensor(list(model_class_weight.values())),
    # validation data
    train_split=predefined_split(
        Dataset(torch.Tensor(x_val.swapaxes(1, 3).swapaxes(2, 3)), torch.Tensor(y_val)))
)

# change the dimension of y_train to 2D


net.fit(torch.Tensor(x_train.swapaxes(1, 3).swapaxes(2, 3)),
        torch.Tensor(y_train))

net.predict_proba(torch.Tensor(
    x_indtest.swapaxes(1, 3).swapaxes(2, 3)))

net.predict(torch.Tensor(
    x_indtest.swapaxes(1, 3).swapaxes(2, 3)))
