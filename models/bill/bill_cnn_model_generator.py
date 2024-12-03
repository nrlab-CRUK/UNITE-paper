from torchvision import models
import torch
import torch.nn as nn


def return_efficientnet(size='small', dev='cpu', in_channels=1, out_channels=3, use_pretrained=True):
    if size == 'small':
        model_generator = models.efficientnet_v2_s
        weights = models.EfficientNet_V2_S_Weights.IMAGENET1K_V1
    elif size == 'medium':
        model_generator = models.efficientnet_v2_m
        weights = models.EfficientNet_V2_M_Weights.IMAGENET1K_V1
    elif size == 'large':
        model_generator = models.efficientnet_v2_l
        weights = models.EfficientNet_V2_L_Weights.IMAGENET1K_V1
    else:
        assert (1 == 2)
    if use_pretrained:
        axial_tile_model = model_generator(weights=weights).to(dev)
    else:
        axial_tile_model = model_generator().to(dev)

    if size == 'large':
        axial_tile_model.features[0][0] = nn.Conv2d(in_channels, 32, kernel_size=(
            3, 3), stride=(2, 2), padding=(1, 1), bias=False).to(dev)
    else:
        axial_tile_model.features[0][0] = nn.Conv2d(in_channels, 24, kernel_size=(
            3, 3), stride=(2, 2), padding=(1, 1), bias=False).to(dev)

    axial_tile_model.classifier[-1] = nn.Linear(
        1280, out_channels, bias=True).to(dev)

    return axial_tile_model.to(dev)


def return_swin(size='small', dev='cpu', in_channels=1, out_channels=3, use_pretrained=True):
    if size == 'small':
        model_generator = models.swin_v2_t
        weights = models.Swin_V2_T_Weights.IMAGENET1K_V1
    elif size == 'medium':
        model_generator = models.swin_v2_s
        weights = models.Swin_V2_S_Weights.IMAGENET1K_V1
    elif size == 'large':
        model_generator = models.swin_v2_b
        weights = models.Swin_V2_B_Weights.IMAGENET1K_V1
    else:
        assert (1 == 2)
    if use_pretrained:
        axial_tile_model = model_generator(weights=weights).to(dev)
    else:
        axial_tile_model = model_generator().to(dev)

    if size == 'large':
        axial_tile_model.features[0][0] = nn.Conv2d(
            in_channels, 128, kernel_size=(4, 4), stride=(4, 4)).to(dev)
        axial_tile_model.head = nn.Linear(
            1024, out_channels, bias=True).to(dev)
    else:
        axial_tile_model.features[0][0] = nn.Conv2d(
            in_channels, 96, kernel_size=(4, 4), stride=(4, 4)).to(dev)
        axial_tile_model.head = nn.Linear(768, out_channels, bias=True).to(dev)

    return axial_tile_model.to(dev)


def return_resnet(size='small', dev='cpu', in_channels=1, out_channels=3, use_pretrained=True):
    if size == 'small':
        model_generator = models.resnet50
        weights = models.ResNet50_Weights.IMAGENET1K_V1
    elif size == 'medium':
        model_generator = models.resnet101
        weights = models.ResNet101_Weights.IMAGENET1K_V1
    elif size == 'large':
        model_generator = models.resnet152
        weights = models.ResNet152_Weights.IMAGENET1K_V1
    else:
        assert (1 == 2)

    if use_pretrained:
        axial_tile_model = model_generator(weights=weights).to(dev)
    else:
        axial_tile_model = model_generator().to(dev)

    axial_tile_model.conv1 = nn.Conv2d(
        in_channels, 64, kernel_size=(4, 4), stride=(4, 4)).to(dev)
    axial_tile_model.fc = nn.Linear(2048, out_channels, bias=True).to(dev)

    return axial_tile_model.to(dev)


def return_resnext(size='small', dev='cpu', in_channels=1, out_channels=3, use_pretrained=True):
    if size == 'small':
        model_generator = models.resnext50_32x4d
        weights = models.ResNeXt50_32X4D_Weights.IMAGENET1K_V1
    elif size == 'medium':
        model_generator = models.resnext101_64x4d
        weights = models.ResNeXt101_64X4D_Weights.IMAGENET1K_V1
    elif size == 'large':
        model_generator = models.resnext101_32x8d
        weights = models.ResNeXt101_32X8D_Weights.IMAGENET1K_V1
    else:
        assert (1 == 2)

    if use_pretrained:
        axial_tile_model = model_generator(weights=weights).to(dev)
    else:
        axial_tile_model = model_generator().to(dev)

    axial_tile_model.conv1 = nn.Conv2d(in_channels, 64, kernel_size=(
        7, 7), stride=(2, 2), padding=(3, 3), bias=False).to(dev)
    axial_tile_model.fc = nn.Linear(2048, out_channels, bias=True).to(dev)

    return axial_tile_model.to(dev)


def return_convnext(size='small', dev='cpu', in_channels=1, out_channels=3, use_pretrained=True):
    if size == 'small':
        model_generator = models.convnext_tiny
        weights = models.ConvNeXt_Tiny_Weights.IMAGENET1K_V1
    elif size == 'medium':
        model_generator = models.convnext_small
        weights = models.ConvNeXt_Small_Weights.IMAGENET1K_V1
    elif size == 'large':
        model_generator = models.convnext_base
        weights = models.ConvNeXt_Base_Weights.IMAGENET1K_V1
    else:
        assert (1 == 2)

    if use_pretrained:
        axial_tile_model = model_generator(weights=weights).to(dev)
    else:
        axial_tile_model = model_generator().to(dev)

    if size == 'large':
        axial_tile_model.features[0][0] = nn.Conv2d(
            in_channels, 128, kernel_size=(4, 4), stride=(4, 4)).to(dev)
        axial_tile_model.classifier[-1] = nn.Linear(
            1024, out_channels, bias=True).to(dev)
    elif size == 'medium':
        axial_tile_model.features[0][0] = nn.Conv2d(
            in_channels, 96, kernel_size=(4, 4), stride=(4, 4)).to(dev)
        axial_tile_model.classifier[-1] = nn.Linear(
            768, out_channels, bias=True).to(dev)
    else:
        axial_tile_model.features[0][0] = nn.Conv2d(
            in_channels, 64, kernel_size=(4, 4), stride=(4, 4)).to(dev)
        axial_tile_model.classifier[-1] = nn.Linear(
            512, out_channels, bias=True).to(dev)

    return axial_tile_model.to(dev)


def get_model(model_name='resnet', size='small', use_pretrained=True, dev='cpu', in_channels=1, out_channels=3):
    acceptable_models = ['resnet', 'resnext',
                         'convnext', 'swin', 'efficientnet']
    if model_name.lower() not in acceptable_models:
        raise ValueError(
            'Model name not recognized, please choose from: {}'.format(acceptable_models))

    if model_name.lower() == 'resnet':
        return return_resnet(size=size, dev=dev, in_channels=in_channels, out_channels=out_channels, use_pretrained=use_pretrained)
    elif model_name.lower() == 'resnext':
        return return_resnext(size=size, dev=dev, in_channels=in_channels, out_channels=out_channels, use_pretrained=use_pretrained)
    elif model_name.lower() == 'convnext':
        return return_convnext(size=size, dev=dev, in_channels=in_channels, out_channels=out_channels, use_pretrained=use_pretrained)
    elif model_name.lower() == 'swin':
        return return_swin(size=size, dev=dev, in_channels=in_channels, out_channels=out_channels, use_pretrained=use_pretrained)
    elif model_name.lower() == 'efficientnet':
        return return_efficientnet(size=size, dev=dev, in_channels=in_channels, out_channels=out_channels, use_pretrained=use_pretrained)


if __name__ == '__main__':
    import numpy as np
    import os
    im = torch.Tensor(np.load(os.path.join(os.getcwd(), 'input.npy'))[
                      :100]).swapaxes(1, 3).swapaxes(2, 3)
    targ = torch.Tensor(np.load(os.path.join(os.getcwd(), 'targ.npy'))[:100])

    # EXPECT INPUT TO HAVE SHAPE BATCH x CHANNELS x HEIGHT x WIDTH, which in your case is BATCH x 3 x 465 x 201
    # POSSIBLE MODEL NAMES ARE: 'resnet','resnext','convnext','swin','efficientnet'
    # POSSIBLE SIZES ARE: 'small','medium','large' - the actual number of parameters will vary depending on the model -
    # but 'small' in one model will have a similar number of parameters as 'small' in all other models. This holds
    # true for all 'medium' and 'large' models as well.

    # this script will automatically download pretrained weights from the internet if use_pretrained_weights is set to
    # True (thank you pytorch). If you want to train the model from scratch, set use_pretrained_weights to False

    in_channels = 3
    out_channels = 1
    use_pretrained_weights = True
    dev = 'cpu'
    model = get_model(model_name='efficientnet', size='small', use_pretrained=use_pretrained_weights,
                      dev=dev, in_channels=in_channels, out_channels=out_channels)
    print(model)

    # scale x to be between 0 and 255
    im = im*255

    model.eval()
    pred = model(im)
    print(pred)
    print(targ)
