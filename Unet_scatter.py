# -*- coding: utf-8 -*-
# This code was witten by Taisei Shimomura (U of Tokushima) for scatter correction on CBCT projection data

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cv2
import os
from scipy import ndimage
from PIL import Image
from numpy import inf
import time
import csv
from skimage.measure import compare_psnr as psnr
from skimage.measure import compare_ssim as ssim
from keras.models import Model
from keras.layers import Input, PReLU,LeakyReLU, BatchNormalization, Activation, Dropout, Concatenate
from keras.layers.convolutional import Conv2D, UpSampling2D, Conv3D, UpSampling3D, ZeroPadding3D
import h5py
from keras.models import model_from_json
from keras.models import load_model

from keras.optimizers import Adam
from keras.callbacks import CSVLogger, ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
import keras.backend as K
import tensorflow as tf
import load_image    #.py
import augmentation  #.py




def Unet2D(input_count, output_count, USE_indiv=False):
    if INPUT_CHANNELS==1: ENERGY = spectral[0]
    elif INPUT_CHANNELS==2: ENERGY = spectral[0] + spectral[1]
    elif INPUT_CHANNELS==3: ENERGY = spectral[0] + spectral[1] + spectral[2]
    else: ENERGY = '4E'
    drop_rate_list = [0.4, 0.3, 0.2, 0.1]

    if USE_indiv==False:
        unet_inputs = Input(shape=(INPUT_HEIGHT, INPUT_WIDTH, input_count), name="unet_input")
    if USE_indiv==True:
        ENERGY = spectral[0]
        unet_inputs = Input(shape=(INPUT_HEIGHT, INPUT_WIDTH, 1), name="unet_input")

    #num_downsample = int(np.floor(np.log(INPUT_HEIGHT)/np.log(2)))#Fujiwara
    #list_filter_count = [FILTER_COUNT*min(8, (2**i)) for i in range(num_downsample)]
    num_downsample = int(np.floor(np.log(INPUT_HEIGHT)/np.log(2)))-1
    list_filter_count = [FILTER_COUNT*min(64, (2**i)) for i in range(num_downsample)]
    
    
    
    
    ##scatter 3convolutions in 1 layer
    name_bn="unet__batchnorm_1st"###add
    enc1 = Conv2D(list_filter_count[0], CONV_FILTER_SIZE, strides=1, padding='same', kernel_initializer='glorot_uniform', name="unet_"+ENERGY+"_conv2D_1_a")(unet_inputs)#scatter
    enc1 = BatchNormalization(name=name_bn+"a")(enc1)###add
    #enc1 = LeakyReLU(0.2)(enc1)
    enc1 =Activation('relu')(enc1)
    #enc1 = PReLU()(enc1)#shimo
    enc2 = Conv2D(list_filter_count[0], CONV_FILTER_SIZE, strides=1, padding='same', kernel_initializer='glorot_uniform', name="unet_"+ENERGY+"_conv2D_1_b")(enc1)#scatter
    enc2 = BatchNormalization(name=name_bn+"b")(enc2)###add
    #enc2 = LeakyReLU(0.2)(enc2)
    enc2 =Activation('relu')(enc2)
    #enc2 = PReLU()(enc2)#shimo
    enc3 = Conv2D(list_filter_count[0], CONV_FILTER_SIZE, strides=1, padding='same', kernel_initializer='glorot_uniform', name="unet_"+ENERGY+"_conv2D_1_c")(enc2)#scatter
    enc3 = BatchNormalization(name=name_bn+"c")(enc3)###add
    #enc3= LeakyReLU(0.2)(enc3)
    enc3 =Activation('relu')(enc3)
    #enc3 = PReLU()(enc3)#shimo
    list_encoder = [enc3]#scatter
    
    
    for i, f in enumerate(list_filter_count[1:]):
        name_conv = "unet_" + ENERGY + "_conv2D_" + str(i+2)
        name_bn = "unet_" + ENERGY + "_batchnorm_" + str(i+1)
        enc = encoding_layer(list_encoder[-1], f, name_conv, name_bn)#def encoding_layer
        list_encoder.append(enc)
        
    
    # Decoder
    list_filter_count = list_filter_count[:-1][::-1]
    if len(list_filter_count) < num_downsample-1:
        list_filter_count.append(FILTER_COUNT)

    dec1 = decoding_layer(list_encoder[-1], list_encoder[-2], list_filter_count[0], "unet_"+ENERGY+"_upconv2D_1", "unet_"+ENERGY+"_upbatchnorm_1", "unet_"+ENERGY+"_upconv2D2_1", dropout=True, num_drop=drop_rate_list[0])
    list_decoder = [dec1]
    for i, f in enumerate(list_filter_count[1:]):
        name_conv = "unet_" + ENERGY + "_upconv2D_" + str(i+2)#i=0,1,2,,,
        name_bn = "unet_" + ENERGY + "_upbatchnorm_" + str(i+2)
        name_conv2 = "unet_" + ENERGY + "_upconv2D2_" + str(i+2)
        if i<3:
            d = True
            num_drop = drop_rate_list[i+1]
        else:
            d = False
        dec = decoding_layer(list_decoder[-1], list_encoder[-(i+3)], f, name_conv, name_bn, name_conv2, dropout=d, num_drop=num_drop)#i=0,1,2,,,  Fujiwara koko
        list_decoder.append(dec)

        
    #scatter 
    #x = Activation('relu')(list_decoder[-1])#scatter
    #x_mse = Conv2D(output_count, 1, padding='same', strides=1,name="unet_"+ENERGY+"_upconv2D_last")(x)#scatter
    x_mse = Conv2D(output_count, 1, padding='same', strides=1,name="unet_"+ENERGY+"_upconv2D_last")(list_decoder[-1])#scatter
    x_mse = Activation('relu')(x_mse)#scatter
    
   
    model = Model(inputs=unet_inputs, outputs=x_mse)
    return model

def encoding_layer(x, filter_count, name_conv, name_bn, USE_3D=False):


    
    ##### Scatter  # 3convolutions in 1 layer encoder
    x = Conv2D(filter_count, CONV_FILTER_SIZE, padding='same', strides=CONV_STRIDE, name=name_conv+"down")(x)
    #x = LeakyReLU(0.2)(x)
    x =Activation('relu')(x)
    #x = PReLU()(x)#shimo
    x = Conv2D(filter_count, CONV_FILTER_SIZE, padding='same', strides=1, kernel_initializer='glorot_uniform', name=name_conv+"a")(x)
    x = BatchNormalization(name=name_bn+"a")(x)###add
    #x = LeakyReLU(0.2)(x)
    x =Activation('relu')(x)
    #x = PReLU()(x)#shimo
    x = Conv2D(filter_count, CONV_FILTER_SIZE, padding='same', strides=1, kernel_initializer='glorot_uniform', name=name_conv+"b")(x)
    x = BatchNormalization(name=name_bn+"b")(x)###add
    #x = LeakyReLU(0.2)(x)
    x =Activation('relu')(x)
    #x = PReLU()(x)#shimo
    x = Conv2D(filter_count, CONV_FILTER_SIZE, padding='same', strides=1, kernel_initializer='glorot_uniform', name=name_conv+"c")(x)
    x = BatchNormalization(name=name_bn+"c")(x)###add
    #x = LeakyReLU(0.2)(x)
    x =Activation('relu')(x)
    #x = PReLU()(x)#shimo
    
    
    return x


    
def decoding_layer(x, encoded_x, filter_count, name_conv, name_bn, name_conv2, dropout=False, num_drop=0.4, USE_3D=False):
    #x = Activation('relu')(x)#Fujiwara
    #x = PReLU()(x)#shimo

    
    if USE_3D==False:

        
        
        x = UpSampling2D(UPSAMPLE_SIZE,interpolation='bilinear')(x)
        x = Concatenate(axis=CONCATENATE_AXIS)([x, encoded_x])
        x = Conv2D(filter_count, DECONV_FILTER_SIZE, strides=1, padding='same', kernel_initializer='glorot_uniform',name=name_conv+"a")(x)
        x = BatchNormalization(name=name_bn+"a")(x)###add
        x = Activation('relu')(x)#Shimo
        #x = PReLU()(x)#shimo
        x = Conv2D(filter_count, DECONV_FILTER_SIZE, strides=1, padding='same', kernel_initializer='glorot_uniform', name=name_conv+"b")(x)
        x = BatchNormalization(name=name_bn+"b")(x)###add
        x = Activation('relu')(x)#Shimo
        #x = PReLU()(x)#shimo
        x = Conv2D(filter_count, DECONV_FILTER_SIZE, strides=1, padding='same', kernel_initializer='glorot_uniform', name=name_conv+"c")(x)
        x = BatchNormalization(name=name_bn+"c")(x)###add
        x = Activation('relu')(x)#Shimo
        #x = PReLU()(x)#shimo
        
        

    if USE_3D==True:
        x = UpSampling3D(UPSAMPLE3D_SIZE)(x)
        x = ZeroPadding3D(CONV3D_PADDING)(x)
        x = Conv3D(filter_count, DECONV3D_FILTER_SIZE, name=name_conv)(x)
    #x = BatchNormalization(name=name_bn)(x)#batch normalize is not used in scatteing
    if dropout: x = Dropout(num_drop)(x)
    if (x.shape[1] != encoded_x.shape[1]):
        x = Conv2D(filter_count, (2,2), padding='valid', strides=1, name=name_conv2)(x)
        #x = BatchNormalization()(x)#batch normalize is not used in scatteing
    #x = Concatenate(axis=CONCATENATE_AXIS)([x, encoded_x])#Fujiwara skipconnection
    return x


def Unet3D(input_count, output_count):
    if INPUT_CHANNELS==1: ENERGY = spectral[0]
    elif INPUT_CHANNELS==2: ENERGY = spectral[0] + spectral[1]
    elif INPUT_CHANNELS==3: ENERGY = spectral[0] + spectral[1] + spectral[2]
    else: ENERGY = '4E'
    drop_rate_list = [0.4, 0.3, 0.2, 0.1]

    unet_input = Input(shape=(INPUT_HEIGHT, INPUT_WIDTH, input_count, 1), name="unet_input")
    num_downsample = int(np.floor(np.log(INPUT_HEIGHT)/np.log(2)))
    list_filter_count = [FILTER_COUNT*min(8, (2**i)) for i in range(num_downsample)]
    
    # Encoder
    enc1 = ZeroPadding3D(CONV3D_PADDING)(unet_input)
    enc1 = Conv3D(list_filter_count[0], CONV3D_FILTER_SIZE, strides=CONV3D_STRIDE, padding='same', name="unet_"+ENERGY+"_conv3D_1")(unet_input)
    list_encoder = [enc1]
    for i, f in enumerate(list_filter_count[1:]):
        name_conv = "unet_" + ENERGY + "_conv3D_" + str(i+2)
        name_bn = "unet_" + ENERGY + "_batchnorm_" + str(i+1)
        enc = encoding_layer(list_encoder[-1], f, name_conv, name_bn, USE_3D=True)
        list_encoder.append(enc)
    
    # Decoder
    list_filter_count = list_filter_count[:-2][::-1]
    if len(list_filter_count) < num_downsample-1:
        list_filter_count.append(FILTER_COUNT)

    dec1 = decoding_layer(list_encoder[-1], list_encoder[-2], list_filter_count[0], "unet_"+ENERGY+"_upconv3D_1", "unet_"+ENERGY+"_upbatchnorm_1", dropout=True, num_drop=drop_rate_list[0], USE_3D=True)
    list_decoder = [dec1]
    for i, f in enumerate(list_filter_count[1:]):
        name_conv = "unet_" + ENERGY + "_upconv3D_" + str(i+2)
        name_bn = "unet_" + ENERGY + "_upbatchnorm_" + str(i+2)
        if i<3:
            d = True
            num_drop = drop_rate_list[i+1]
        else:
            d = False
        dec = decoding_layer(list_decoder[-1], list_encoder[-(i+3)], f, name_conv, name_bn, dropout=d, num_drop=num_drop, USE_3D=True)#d=true or false
        list_decoder.append(dec)

    x = Activation('relu')(list_decoder[-1])
    x = UpSampling3D(UPSAMPLE3D_SIZE)(x)
    x = ZeroPadding3D(CONV3D_PADDING)(x)
    x_mse = Conv3D(output_count, (3,3,input_count), name="unet_"+ENERGY+"_upconv3D_last")(x)### should be changed

    model = Model(inputs=unet_input, outputs=x_mse)
    return model


def ssim_acc(y_true, y_pred):
    return tf.reduce_mean(tf.image.ssim(y_true, y_pred, max_val=1.0))

def ssim_loss(y_true, y_pred):
    return 1.0 - ssim_acc(y_true, y_pred)

def dice_coef(y_true, y_pred):
    y_true = K.flatten(y_true)
    y_pred = K.flatten(y_pred)
    intersection = K.sum(y_true * y_pred)    
    return (2.0 * intersection + 1) / (K.sum(y_true) + K.sum(y_pred) + 1)

def dice_coef_loss(y_true, y_pred):
    return 1.0 - dice_coef(y_true, y_pred)

def total_acc(y_true, y_pred):
    pred = K.cast(K.greater_equal(y_pred, 0.5), "float")
    flag = K.cast(K.equal(y_true, pred), "float")
    return K.prod(flag, axis=-1)


def train_step(X, Y, batch_size, LR, file_weight, file_csv, aug_flip, aug_RICAP, aug_scale, aug_rotate, USE_3D=False):
    X_train_set, X_val_set, Y_train_set, Y_val_set = load_image.split_train_val(X, Y, num_validation=10)
    X_train, Y_train = augmentation.generate(X_train_set, Y_train_set, aug_flip, aug_RICAP, aug_scale, aug_rotate)
    X_val, Y_val = augmentation.generate(X_val_set, Y_val_set, aug_flip, aug_RICAP, aug_scale, aug_rotate)
    if USE_3D==True:
        X_train = np.reshape(X_train, [X_train.shape[0], X_train.shape[1], X_train.shape[2], X_train.shape[3], 1])
        Y_train = np.reshape(Y_train, [Y_train.shape[0], Y_train.shape[1], Y_train.shape[2], 1, Y_train.shape[3]])
        X_val = np.reshape(X_val, [X_val.shape[0], X_val.shape[1], X_val.shape[2], X_val.shape[3], 1])
        Y_val = np.reshape(Y_val, [Y_val.shape[0], Y_val.shape[1], Y_val.shape[2], 1, Y_val.shape[3]])
    cp_cb = ModelCheckpoint(filepath = file_weight, monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=True, mode='auto')
    cl_cb = CSVLogger(file_csv, separator=',', append=True)
    lr_cb = ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=0.00001)
    train = model.fit(X_train, Y_train, batch_size=batch_size, epochs=n_epoch*30, validation_data=(X_val, Y_val), callbacks=[cp_cb, cl_cb, lr_cb])

def train_step_cv(X, Y, batch_size,  LR, file_weight, file_csv, aug_flip, aug_RICAP, aug_scale, aug_rotate, USE_indiv=False, USE_3D=False):
    epochcounter = 0
    last_loss, min_loss = 100000, 100000
    patience = 5 #if counter is higher than patience, LR would be (x0.5) 
    counter = 0
    for k in range(30):
        #X_train_set, X_val_set, Y_train_set, Y_val_set = load_image.split_train_val(X, Y, num_validation=10*4)#train & validation
        #X_train_set, X_val_set, Y_train_set, Y_val_set = load_image.split_train_val(X, Y, num_validation=10*4*2)#train & validation
        #X_train_set, X_val_set, Y_train_set, Y_val_set = load_image.split_train_val(X, Y, num_validation=90)#train & validation
        #X_train_set, X_val_set, Y_train_set, Y_val_set = load_image.split_train_val(X, Y, num_validation=600)#train & validation
        
        ########paper 30% validation!
        X_train_set, X_val_set, Y_train_set, Y_val_set = load_image.split_train_val(X, Y, num_validation=810)#train & validation
        
        X_train, Y_train = augmentation.generate(X_train_set, Y_train_set, aug_flip, aug_RICAP, aug_scale, aug_rotate)
        X_val, Y_val = augmentation.generate(X_val_set, Y_val_set, aug_flip, aug_RICAP, aug_scale, aug_rotate)

        if USE_indiv==True:
            X_train_list = [X_train[:,:,:,0,np.newaxis], X_train[:,:,:,1,np.newaxis]]
            X_val_list = [X_val[:,:,:,0,np.newaxis], X_val[:,:,:,1,np.newaxis]]
            if INPUT_CHANNELS > 2:
                X_train_list.append(X_train[:,:,:,2,np.newaxis])
                X_val_list.append(X_val[:,:,:,2,np.newaxis])
            if INPUT_CHANNELS > 3:
                X_train_list.append(X_train[:,:,:,3,np.newaxis])
                X_val_list.append(X_val[:,:,:,3,np.newaxis])

        if USE_3D==True:
            X_train = np.reshape(X_train, [X_train.shape[0], X_train.shape[1], X_train.shape[2], X_train.shape[3], 1])
            Y_train = np.reshape(Y_train, [Y_train.shape[0], Y_train.shape[1], Y_train.shape[2], 1, Y_train.shape[3]])
            X_val = np.reshape(X_val, [X_val.shape[0], X_val.shape[1], X_val.shape[2], X_val.shape[3], 1])
            Y_val = np.reshape(Y_val, [Y_val.shape[0], Y_val.shape[1], Y_val.shape[2], 1, Y_val.shape[3]])

        for o in range(n_epoch):
            if counter >= patience:
                LR = LR * 0.5
                #model.compile(loss='mse', optimizer=Adam(lr=LR), metrics=['accuracy'])
                model.compile(loss='mae', optimizer=Adam(lr=LR), metrics=['accuracy'])
                #model.compile(loss='mape', optimizer=Adam(lr=LR), metrics=['accuracy'])#shimo for Unet scattering
                counter = 0
                patience += 1
                last_loss = 100000
            print(epochcounter, LR, last_loss)
            #cp_cb = ModelCheckpoint(filepath = file_weight, monitor='val_loss', verbose=0, save_best_only=True, mode='auto')
            cl_cb = CSVLogger(file_csv, separator=',', append=True)
            #lr_cb = ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=0.00001)
            if USE_indiv==False: train = model.fit(X_train, Y_train, batch_size=batch_size, epochs=1, validation_data=(X_val, Y_val), callbacks=[cl_cb])
            if USE_indiv==True: train = model.fit(X_train_list, Y_train, batch_size=batch_size, epochs=1, validation_data=(X_val_list, Y_val), callbacks=[cl_cb])
            val_loss = min(train.history['val_loss'])
            if val_loss > last_loss:
                last_loss = val_loss
                counter += 1
            else:
                last_loss = val_loss
                counter = 0
                #if epochcounter > 50 and min_loss >= last_loss:#Fujiwara
                if epochcounter > 20 and min_loss >= last_loss:#shimo
                    min_loss = last_loss
                    model.save_weights(file_weight)
            epochcounter += 1

def plot_history(file_csv, plot_loss, plot_acc, plot_range):
    df = pd.read_csv(file_csv, index_col=[0])
    plt.plot(range(len(df['loss'])), df['loss'], label='train_loss')
    plt.plot(range(len(df['val_loss'])), df['val_loss'], label='val_loss')
    val = np.array(df['val_loss'])
    print(np.argmin(val)+1, len(val))
    plt.scatter(x=np.argmin(val), y=min(val), marker='o', c='k')
    plt.ylim(0.0, plot_range)
    plt.title('model loss')
    plt.xlabel('epoch')
    plt.ylabel('loss')
    plt.legend(loc='upper left')
    plt.savefig(plot_loss)
    plt.close('all')

    plt.plot(range(len(df['acc'])), df['acc'], label='train_acc')
    plt.plot(range(len(df['val_acc'])), df['val_acc'], label='val_acc')
    plt.title('model accuracy')
    plt.xlabel('epoch')
    plt.ylabel('accuracy')
    plt.legend(loc='upper left')
    plt.savefig(plot_acc)

#def predict(X_test, Y_test, INPUT_ELEMENTS, file_model, file_weight, folder_name, model_name,validate, save_img=False, USE_WholeBody=False, USE_indiv=False, USE_3D=False):
    
def predict(X_test, Y_test, INPUT_ELEMENTS, file_model, file_weight, folder_name, model_name,validate,gender,age,height,weight,save_img=False, USE_WholeBody=False, USE_indiv=False, USE_3D=False):# for test 2
    
    RMSE = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    RMSPE = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    MAE = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    PSNR = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    SSIM = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    RMSE_bbox = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    RMSPE_bbox = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    MAE_bbox = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    PSNR_bbox = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    SSIM_bbox = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    RMSE_bbox_in = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    RMSPE_bbox_in = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    MAE_bbox_in = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    PSNR_bbox_in = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    SSIM_bbox_in = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    RMSE_mid = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    RMSPE_mid = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    MAE_mid = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    PSNR_mid = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    SSIM_mid = np.zeros([Y_test.shape[0], INPUT_ELEMENTS])
    if USE_indiv==True:
        X_test_list = [X_test[:,:,:,0,np.newaxis], X_test[:,:,:,1,np.newaxis]]
        if INPUT_CHANNELS > 2:
            X_test_list.append(X_test[:,:,:,2,np.newaxis])
        if INPUT_CHANNELS > 3:
            X_test_list.append(X_test[:,:,:,3,np.newaxis])
    if USE_3D==True: X_test = np.reshape(X_test, [X_test.shape[0], X_test.shape[1], X_test.shape[2], X_test.shape[3], 1])

    model = model_from_json(open(file_model).read())
    model.load_weights(file_weight)
    start_time = time.time()
    if USE_indiv==False: predicted = model.predict(X_test)
    if USE_indiv==True: predicted = model.predict(X_test_list)
    print("time:%lf" %(time.time() - start_time))
    

    for Np in range(Y_test.shape[0]):
        img_ref = np.sum(Y_test[Np,:,:,:], axis=-1)
        for mat in range(INPUT_ELEMENTS):
            img_pred = predicted[Np,:,:,mat]
            img_true = Y_test[Np,:,:,mat]
            RMSE[Np,mat] = rmse(img_true.astype('float64'), img_pred.astype('float64'))
            RMSPE[Np,mat] = rmspe(img_true.astype('float64'), img_pred.astype('float64'))
            MAE[Np,mat] = mae(img_true.astype('float64'), img_pred.astype('float64'))
            PSNR[Np,mat] = psnr(img_true.astype('float64'), img_pred.astype('float64'))
            SSIM[Np,mat] = ssim(img_true.astype('float64'), img_pred.astype('float64'))
            RMSE_bbox[Np,mat] = rmse_bbox(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            RMSPE_bbox[Np,mat] = rmspe_bbox(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            MAE_bbox[Np,mat] = mae_bbox(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            PSNR_bbox[Np,mat] = psnr_bbox(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            SSIM_bbox[Np,mat] = ssim_bbox(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            RMSE_bbox_in[Np,mat] = rmse_bbox_inner(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            RMSPE_bbox_in[Np,mat] = rmspe_bbox_inner(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            MAE_bbox_in[Np,mat] = mae_bbox_inner(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            PSNR_bbox_in[Np,mat] = psnr_bbox_inner(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            SSIM_bbox_in[Np,mat] = ssim_bbox_inner(img_ref, img_true.astype('float64'), img_pred.astype('float64'))
            RMSE_mid[Np,mat] = rmse_mid(img_ref, img_true.astype('float64'), img_pred.astype('float64'), 50)
            RMSPE_mid[Np,mat] = rmspe_mid(img_ref, img_true.astype('float64'), img_pred.astype('float64'), 50)
            MAE_mid[Np,mat] = mae_mid(img_ref, img_true.astype('float64'), img_pred.astype('float64'), 50)
            PSNR_mid[Np,mat] = psnr_mid(img_ref, img_true.astype('float64'), img_pred.astype('float64'), 50)
            SSIM_mid[Np,mat] = ssim_mid(img_ref, img_true.astype('float64'), img_pred.astype('float64'), 50)
            
    if save_img==True:
        
        ##### paper
        for iz in range(predicted.shape[0]):
            img_res = Y_test[0,:,:,:] - predicted[0,:,:,:]
            img_res = np.transpose(img_res, (2, 0, 1))
            #img_res.astype('float32').tofile('%s/res_image/%s%dy%dh%dw/%s_res_float_%03d.raw' % (folder_name, gender, age, height,weight,model_name,iz))
            imgout = predicted[iz,:,:,:]
            imgout = np.transpose(imgout, (2, 0, 1))
            #imgout.astype('float32').tofile('%s/predicted_scatter/%s_pred_float_%03d.raw' % (folder_name, model_name,iz))
            imgout.astype('float32').tofile('%s/predicted_scatter/%s%dy%dh%dw/%s_pred_float_%03d.raw' % (folder_name, gender, age, height,weight, model_name,iz))###### for predict2
        

    df1 = pd.DataFrame(RMSE, columns=['H'])
    df2 = pd.DataFrame(RMSPE, columns=['H'])
    df3 = pd.DataFrame(MAE, columns=['H'])
    df4 = pd.DataFrame(PSNR, columns=['H'])
    df5 = pd.DataFrame(SSIM, columns=['H'])
    df1_ = pd.DataFrame(RMSE_bbox, columns=['H'])
    #df11_ = pd.DataFrame(RMSE_mid, columns=['H'])
    df2_ = pd.DataFrame(RMSPE_bbox, columns=['H'])
    df3_ = pd.DataFrame(MAE_bbox, columns=['H'])
    df4_ = pd.DataFrame(PSNR_bbox, columns=['H'])
    df5_ = pd.DataFrame(SSIM_bbox, columns=['H'])
    #df44_ = pd.DataFrame(SSIM_mid, columns=['H'])
    df1__ = pd.DataFrame(RMSE_bbox_in, columns=['H'])
    df2__ = pd.DataFrame(RMSPE_bbox_in, columns=['H'])
    df3__ = pd.DataFrame(MAE_bbox_in, columns=['H'])
    df4__ = pd.DataFrame(PSNR_bbox_in, columns=['H'])
    df5__ = pd.DataFrame(SSIM_bbox_in, columns=['H'])
    df11 = pd.DataFrame(RMSE_mid, columns=['H'])
    df22 = pd.DataFrame(RMSPE_mid, columns=['H'])
    df33 = pd.DataFrame(MAE_mid, columns=['H'])
    df44 = pd.DataFrame(PSNR_mid, columns=['H'])
    df55 = pd.DataFrame(SSIM_mid, columns=['H'])

    
    df1.to_csv('%s/%s%s_rmse_table.csv' % (folder_name, model_name, validate))
    df2.to_csv('%s/%s%s_rmspe_table.csv' % (folder_name, model_name, validate))
    df3.to_csv('%s/%s%s_mae_table.csv' % (folder_name, model_name, validate))
    df4.to_csv('%s/%s%s_psnr_table.csv' % (folder_name, model_name, validate))
    df5.to_csv('%s/%s%s_ssim_table.csv' % (folder_name, model_name, validate))
    df1_.to_csv('%s/%s%s_rmse_bbox_table.csv' % (folder_name, model_name, validate))
    df2_.to_csv('%s/%s%s_rmspe_bbox_table.csv' % (folder_name, model_name, validate))
    df3_.to_csv('%s/%s%s_mae_bbox_table.csv' % (folder_name, model_name, validate))
    df4_.to_csv('%s/%s%s_psnr_bbox_table.csv' % (folder_name, model_name, validate))
    df5_.to_csv('%s/%s%s_ssim_bbox_table.csv' % (folder_name, model_name, validate))
    df1__.to_csv('%s/%s%s_rmse_bbox_inner_table.csv' % (folder_name, model_name, validate))
    df2__.to_csv('%s/%s%s_rmspe_bbox_inner_table.csv' % (folder_name, model_name, validate))
    df3__.to_csv('%s/%s%s_mae_bbox_inner_table.csv' % (folder_name, model_name, validate))
    df4__.to_csv('%s/%s%s_psnr_bbox_inner_table.csv' % (folder_name, model_name, validate))
    df5__.to_csv('%s/%s%s_ssim_bbox_inner_table.csv' % (folder_name, model_name, validate))
    df11.to_csv('%s/%s%s_rmse_mid_table.csv' % (folder_name, model_name, validate))
    df22.to_csv('%s/%s%s_rmspe_mid_table.csv' % (folder_name, model_name, validate))
    df33.to_csv('%s/%s%s_mae_mid_table.csv' % (folder_name, model_name, validate))
    df44.to_csv('%s/%s%s_psnr_mid_table.csv' % (folder_name, model_name, validate))
    df55.to_csv('%s/%s%s_ssim_mid_table.csv' % (folder_name, model_name, validate))

def rmse(img_a, img_b):
    err = (np.sum( (img_a - img_b)*(img_a - img_b)/(img_a.shape[0]*img_a.shape[1]) ))**0.5
    return err

def rmspe(img_a, img_b):
    #"""
    with np.errstate(divide='ignore', invalid='ignore'):
        err = (img_a - img_b)/img_a
        err[err==inf] = 0
        err[err==-inf] = 0
        err = np.nan_to_num(err)
        err[err==np.nan] = 0
    #"""
    #err = (img_a - img_b)/(img_a+0.00001)
    err = ( np.sum(err**2)/(img_a.shape[0]*img_a.shape[1]) )**0.5
    return err*100

def mae(img_a, img_b):
    err = np.sum( np.abs(img_a - img_b)/(img_a.shape[0]*img_a.shape[1]) )
    return err


def find_outer_rectangle(img):
    cols = np.any(img > 0, axis=1)
    rows = np.any(img > 0, axis=0)
    y_min, y_max = np.where(cols)[0][[0, -1]]
    x_min, x_max = np.where(rows)[0][[0, -1]]
    return y_min, y_max, x_min, x_max

def find_inner_rectangle(img):
    img[img>0] = 1.0
    y_min, y_max, x_min, x_max = find_outer_rectangle(img)
    x_mid, y_mid = int((x_min+x_max)/2), int((y_min+y_max)/2)
    img = ndimage.binary_closing(img, structure=np.ones((16,16)))
    # because img.shape[0]/16 = 16
    min_exist_ratio = 0.0
    min_nn_1 = 0
    for j in range(16, y_mid-y_min):
        for i in range(16, x_mid-x_min):
            temp_img_bbox = img[y_mid-j:y_mid+j, x_mid-i:x_mid+i]
            nn_0 = temp_img_bbox.size - np.count_nonzero(temp_img_bbox)
            nn_1 = np.count_nonzero(temp_img_bbox)
            exist_ratio = (nn_1 - nn_0)/nn_1
            if exist_ratio>=min_exist_ratio and nn_1>min_nn_1:
                min_exist_ratio = exist_ratio
                min_nn_1 = nn_1
                best_i = i
                best_j = j
    return y_mid-best_j, y_mid+best_j, x_mid-best_i, x_mid+best_i

def find_mask_rectangle(img, mask):
    y_min, y_max, x_min, x_max = find_outer_rectangle(img)
    y_mask_min = int((y_min+y_max)/2) - int(mask/2)
    y_mask_max = int((y_min+y_max)/2) + int(mask/2)
    x_mask_min = int((x_min+x_max)/2) - int(mask/2)
    x_mask_max = int((x_min+x_max)/2) + int(mask/2)
    return y_mask_min, y_mask_max, x_mask_min, x_mask_max


def rmse_bbox(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_outer_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    rmsebb = rmse(img_a_bbox, img_b_bbox)
    return rmsebb

def rmse_bbox_inner(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_inner_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    rmsebb = rmse(img_a_bbox, img_b_bbox)
    return rmsebb

def rmse_mid(img_ref, img_a, img_b, mask):
    y_mask_min, y_mask_max, x_mask_min, x_mask_max = find_mask_rectangle(img_ref, mask)
    ROI_a = img_a[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    ROI_b = img_b[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    rmsemid = rmse(ROI_a, ROI_b)
    return rmsemid


def rmspe_bbox(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_outer_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    rmspebb = rmspe(img_a_bbox, img_b_bbox)
    return rmspebb

def rmspe_bbox_inner(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_inner_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    rmspebb = rmspe(img_a_bbox, img_b_bbox)
    return rmspebb

def rmspe_mid(img_ref, img_a, img_b, mask):
    y_mask_min, y_mask_max, x_mask_min, x_mask_max = find_mask_rectangle(img_ref, mask)
    ROI_a = img_a[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    ROI_b = img_b[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    rmspemid = rmspe(ROI_a, ROI_b)
    return rmspemid


def mae_bbox(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_outer_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    maebb = mae(img_a_bbox, img_b_bbox)
    return maebb

def mae_bbox_inner(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_inner_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    maebb = mae(img_a_bbox, img_b_bbox)
    return maebb

def mae_mid(img_ref, img_a, img_b, mask):
    y_mask_min, y_mask_max, x_mask_min, x_mask_max = find_mask_rectangle(img_ref, mask)
    ROI_a = img_a[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    ROI_b = img_b[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    maemid = mae(ROI_a, ROI_b)
    return maemid


def psnr(img_a, img_b):
    err = np.sum( (img_a - img_b)*(img_a - img_b) )/(img_a.shape[0]*img_a.shape[1])
    return 10 * np.log10((1.0 ** 2) / err)

def psnr_bbox(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_outer_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    psnrbb = psnr(img_a_bbox, img_b_bbox)
    return psnrbb

def psnr_bbox_inner(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_inner_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    psnrbb = psnr(img_a_bbox, img_b_bbox)
    return psnrbb

def psnr_mid(img_ref, img_a, img_b, mask):
    y_mask_min, y_mask_max, x_mask_min, x_mask_max = find_mask_rectangle(img_ref, mask)
    ROI_a = img_a[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    ROI_b = img_b[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    psnrmid = psnr(ROI_a, ROI_b)
    return psnrmid


def ssim_bbox(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_outer_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    ssimbb = ssim(img_a_bbox, img_b_bbox)
    return ssimbb

def ssim_bbox_inner(img_ref, img_a, img_b):
    y_min, y_max, x_min, x_max = find_inner_rectangle(img_ref)
    img_a_bbox = img_a[y_min:y_max, x_min:x_max]
    img_b_bbox = img_b[y_min:y_max, x_min:x_max]
    ssimbb = ssim(img_a_bbox, img_b_bbox)
    return ssimbb

def ssim_mid(img_ref, img_a, img_b, mask):
    y_mask_min, y_mask_max, x_mask_min, x_mask_max = find_mask_rectangle(img_ref, mask)
    ROI_a = img_a[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    ROI_b = img_b[y_mask_min:y_mask_max, x_mask_min:x_mask_max]
    ssimmid = ssim(ROI_a, ROI_b)
    return ssimmid


def gaussiannoise0(img_array, sd):
    noise = np.random.normal(0, sd, img_array.shape)
    outimg = img_array + noise
    return outimg

def gaussiannoise0_4d(img, sd):
    outimg = np.zeros(img.shape)
    for c in range(img.shape[3]):
        for k in range(img.shape[0]):
            for i in range(img.shape[1]):
                for j in range(img.shape[2]):
                    outimg[k,i,j,c] = np.random.normal(img[k,i,j,c], sd, 1)
    return outimg

if __name__ == '__main__':
    np.random.seed(10)
    
    IMG_HEIGHT, IMG_WIDTH = 400, 400#detector size
    INPUT_HEIGHT, INPUT_WIDTH = 256, 256#load_image.py->def load_Raw-> cv2.resize  (this is resize for down)
    #INPUT_HEIGHT, INPUT_WIDTH = 400, 400#load_image.py->def load_Raw-> cv2.resize
    VAL_HEIGHT, VAL_WIDTH = 512, 512# not important now(on 2022/1/6)
    spectral = ['100kV']
    
    element = ['H']###anything is ok.
    
   
    #### U-net scattering 
    INPUT_CHANNELS = len(spectral)
    INPUT_ELEMENTS = len(element)
    FILTER_COUNT = 16 #first convolution filter shimo
    #FILTER_COUNT = 64 #first convolution filter shimo
    CONCATENATE_AXIS = -1
    CONV_FILTER_SIZE = 3####shimo
    CONV_STRIDE = 2
    UPSAMPLE_SIZE = (2, 2)
    DECONV_FILTER_SIZE = 3
    USE_WholeBody = False ########## switching whole body/30 slices for validation ##########
    
    ### switch the following only in the case of multiple inputs ###
    USE_indiv = False ########## change!!! switching mix/individual 2D network ##########
    USE_3D = False ########## change!!! switching 2D/3D network ##########
    #USE_noise = True ########## change!!! added noise in image-domain ##########
    USE_noise = False ########## change!!! added noise in image-domain ##########
    CONV3D_FILTER_SIZE = (4, 4, 1)
    CONV3D_STRIDE = (2, 2, 1)
    CONV3D_PADDING = (1, 1, 0)
    UPSAMPLE3D_SIZE = (2, 2, 1)
    DECONV3D_FILTER_SIZE = (3, 3, 1)

    

    
    ###In train
    model_name = 'Unet_RICAP_scatter' ########## change!!! ##########
    folder_name = '/home/dl-box/U-net/scatter_2D_unet/' + model_name
    file_model = folder_name + '/test_' + model_name + '_model_architecture.json'
    file_weight = folder_name + '/test_' + model_name + '_model_weights.hdf5'
    
    file_csv = folder_name + '/id_SL_' + model_name + '_history.csv'
    plot_loss = folder_name + '/id_SL_' + model_name + '_plot_loss.png'
    plot_acc = folder_name + '/id_SL_' + model_name + '_plot_accuracy.png'



    
    
    #####In TRAINING  Phase#########  
    X = load_image.load_Raw('*******', 'float32', 1, IMG_HEIGHT, IMG_WIDTH, INPUT_HEIGHT, 'human')#direct projection


    print("1st data input finished")
    if USE_noise==True: X = gaussiannoise0_4d(X, 0.015)
    print("1st data input gaussian noise on")
    if INPUT_CHANNELS>1:
        X2 = load_image.load_Raw('FBP_SL_%s_RICAP_BHC_1_200/train' % spectral[1], 'float32', 1, IMG_HEIGHT, IMG_WIDTH, INPUT_HEIGHT, 'human')
        if USE_noise==True: X2 = gaussiannoise0_4d(X2, 0.005)
        X = np.concatenate([X, X2], axis=CONCATENATE_AXIS)
    if INPUT_CHANNELS>2:
        X3 = load_image.load_Raw('FBP_SL_%s_RICAP_BHC_1_200/train' % spectral[2], 'float32', 1, IMG_HEIGHT, IMG_WIDTH, INPUT_HEIGHT, 'human')
        if USE_noise==True: X3 = gaussiannoise0_4d(X3, 0.015)
        X = np.concatenate([X, X3], axis=CONCATENATE_AXIS)
    if INPUT_CHANNELS>3:
        X4 = load_image.load_Raw('FBP_SL_%s_RICAP_BHC_1_370/train' % spectral[3], 'float32', 1, IMG_HEIGHT, IMG_WIDTH, INPUT_HEIGHT, 'human')
        if USE_noise==True: X4 = gaussiannoise0_4d(X4, 0.005)
        X = np.concatenate([X, X4], axis=CONCATENATE_AXIS)

   
    Y = load_image.load_Raw('********', 'float32', 1, IMG_HEIGHT, IMG_WIDTH, INPUT_HEIGHT, 'human')#scattered images
    
   
    batch_size = 4#
    n_epoch = 30 
    LR = 0.0004 
  
    
    
    if USE_3D==False: model = Unet2D(INPUT_CHANNELS, INPUT_ELEMENTS, USE_indiv)
    if USE_3D==True: model = Unet3D(INPUT_CHANNELS, INPUT_ELEMENTS)
    model.summary()
    #plot_model(model, to_file='%s/id_SL_Unet_model.png' % folder_name,show_shapes=True)
    #model.compile(loss='mse', optimizer=Adam(lr=LR), metrics=['accuracy'])
    model.compile(loss='mae', optimizer=Adam(lr=LR), metrics=['accuracy'])#shimo for U-net scattering
    json_string = model.to_json()
    open(file_model,'w').write(json_string)



    #######paper!!!!#####
    train_step_cv(X, Y, batch_size, LR, file_weight, file_csv, 1, 0, 0, 0, USE_indiv, USE_3D)
    plot_history(file_csv, plot_loss, plot_acc, 0.0025)


   
    
    
    

