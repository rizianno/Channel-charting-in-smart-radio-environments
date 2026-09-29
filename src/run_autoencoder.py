# Step 1: Import Libraries
import tensorflow as tf
from tensorflow.keras import layers, models
import numpy as np
# import matplotlib.pyplot as plt
from scipy import io
from keras.layers import Dense, Conv1DTranspose
from tensorflow.keras.optimizers import Adam
from sklearn.utils.class_weight import compute_class_weight
import numpy as np
from tensorflow.keras.layers import BatchNormalization
from scipy.linalg import logm

print("AE")

# step 1.5
############## 2800 #########################
# file_path = "/home/reza/jupyter/"
# file_path_1 = "/home/reza/jupyter/scenario_position/"
# dataset = io.loadmat(file_path+'R_Cov_BS12NovLOSNLOS_blocker_fullarea_2RIS4phase1bounceuser_3200.mat')
# R_Cov_BS = dataset['R_Cov_BS']

# aa = io.loadmat(file_path_1+'R_xy_LOSNLOS12Nov_3200_blocker_fullarea.mat')
# R_xy = aa['R_xy']

# bb = io.loadmat(file_path_1+'index_uncoveredUEByRIS_block_fullarea.mat')
# count = bb['count']
# count -=1 # convert to python form
###########################################
# Load MATLAB file
file_path = "/home/reza/jupyter/"
file_path_1 = "/home/reza/jupyter/scenario_position/"
file_path_2 = "/home/reza/jupyter/ChannelRISCode/"
# dataset = io.loadmat(file_path+
# 'R_Cov_BS126NovLOSNLOS_blocker_fullarea_2RIS_4phase_uniform32_8noise_RISy4z5751bounceuser_320019tx_power.mat')

# dataset = io.loadmat(file_path+
# 'R_Cov_BS126NovLOSNLOS_blocker_fullarea_2RIS_mirror_uniform32_8noise_RISy4z5751bounceuser_320019tx_power.mat')

# dataset = io.loadmat(file_path+
# 'R_Cov_BS126NovLOSNLOS_blocker_fullarea_NoRIS_uniform32_8noise_RISy4z5751bounceuser_320019tx_power.mat')

dataset = io.loadmat(
'input_data.mat')


# dataset = io.loadmat(file_path+
# 'R_Cov_BS16NovLOSNLOS_blocker_fullarea_2RIS_4phase1bounceuser_9900.mat')

# dataset_1 = io.loadmat(file_path_1+
#                       'R_xy_LOSNLOS15Nov_9900_blocker.mat')
# dataset_2 = io.loadmat(file_path_1+
# 'index_uncoveredUEByRIS_block_fullarea9900.mat')

# dataset_1 = io.loadmat(file_path_2+
#     'DismilarityCLogEuc_LEO_NLOS_blocker_BS1_26Nov_filledpoints_1pbuilding_1bounce_NLOS_fullarea_user3200uniform_2RIS_4phase_org_8upa_noise_middleRISy4z575.mat')
R_Cov_BS = dataset['R_Cov_BS']
# aa = io.loadmat(file_path_1+'R_xy_LOSNLOS12Nov_3200_blocker_fullarea.mat')

# aa = io.loadmat(file_path_1+'R_xy_LOSNLOS15Nov_9900_blocker.mat')
R_xy = dataset['R_xy']
ll = dataset['Supervision_percentage']
# D_Diss = dataset_1['D_LogEuc']
# bb = io.loadmat(file_path_1+'index_uncoveredUEByRIS_block_fullarea.mat')
# # bb = io.loadmat(file_path_1+'index_uncoveredUEByRIS_block_fullarea9900.mat')
# count = dataset_2['count']
# count -=1 # convert to python form



np.random.seed(42)
##### data sepratation to train and test for semi-supervised laearning
# Assuming X_train and X_test contain the input images

# train_data, test_data = train_test_split(features, test_size=0.15, random_state=42) # samples
# list_ = [i for i in range(3000)]
# R_xy = np.delete(R_xy,count,axis=0)
# R_xy = np.delete(R_xy,list_,axis=0)

# R_Cov_BS = np.delete(R_Cov_BS,count,axis=0)
# R_Cov_BS = np.delete(R_Cov_BS,list_,axis=0)

R_logm = np.array([logm(R_Cov_BS[i]) for i in range(R_Cov_BS.shape[0])])

############# cov matrix
# main_diagonals = np.array([np.diag(R_Cov_BS[i]) for i in range(R_Cov_BS.shape[0])]).real
upper_triangles_maindiag = np.array(
    [R_logm[i][np.triu_indices(R_logm.shape[1], k=0)] for i in range(R_logm.shape[0])])
# (np.log(main_diagonals+1e-20),
# upper_triangles = np.abs  , (np.log1p(np.abs(upper_triangles)) np.log1p(main_diagonals)
# normalized_angle = np.angle(upper_triangles)   # Map [-pi, pi] to [0, 1]

R_Cov_BS = np.concatenate((upper_triangles_maindiag.real,
                           upper_triangles_maindiag.imag), axis=1)

# R_Cov_BS = np.concatenate((np.log(np.abs(upper_triangles)+1e-40),
#                            normalized_angle),axis=1)
# R_Cov_BS = np.log(np.abs(upper_triangles)+1e-40)
# R_Cov_BS = np.angle(upper_triangles)



length_feature = R_Cov_BS.shape[1]
# ll = int(.15*2800)


# reshapping
# features = D_Diss+1e-40
features = R_Cov_BS
# features = np.log1p(R_Cov_BS.reshape(len(R_Cov_BS),length_feature))
# assert features.shape == R_Cov_BS.shape
locations = R_xy



l_mean = locations.mean(axis=0)
l_std = locations.std(axis=0)
f_mean = features.mean(axis=0)
f_std = features.std(axis=0)
# # print(f_std)
# # Normalization step
locations = (locations - l_mean) / l_std
features = (features - f_mean) / f_std

# features = np.concatenate((features,normalized_angle),axis=1)

l_min = locations.min(axis=0)
l_max = locations.max(axis=0)
f_min = features.min(axis=0)
f_max = features.max(axis=0)
# Normalization step
# locations = (locations - l_min) / (l_max - l_min)
# features = (features - f_min) / (f_max - f_min)



# Step 2: Generate Synthetic Data
np.random.seed(42)
X = features#[:100,:]  # Input data (500 samples, 10 features)
Y = locations#[:100,:]   # Target labels (500 samples, 2D space for dimensionality reduction)

# Replace with your real data
# X = np.random.rand(1000, num_features).astype(np.float32)  # Example unlabeled data
# y = np.random.rand(100, num_features).astype(np.float32)   # Example anchor data (labeled)

features = R_Cov_BS
# features = np.log1p(R_Cov_BS.reshape(len(R_Cov_BS),length_feature))
# assert features.shape == R_Cov_BS.shape
locations = R_xy

X = features#[:100,:]  # Input data (500 samples, 10 features)
Y = locations#[:100,:]   # Target labels (500 samples, 2D space for dimensionality reduction)


# Combine and split data
# X_train, X_val = train_test_split(X, test_size=0.2, random_state=42)
# ll = .15

# X_train = X[int(ll*3200):,:]
# X_labeled = X[:int(ll*3200),:]
# Y_labeled = locations[:int(ll*3200),:]
# Y_eval = locations[int(ll*3200):,:]

# Y_labeled_mean = Y_labeled.mean(axis=0)
# Y_labeled_std = Y_labeled.std(axis=0)
# l_mean = Y_eval.mean(axis=0)
# l_std = Y_eval.std(axis=0)
# x_labeled_mean = X_labeled.mean(axis=0)
# x_labeled_std = X_labeled.std(axis=0)
# x_mean = X_train.mean(axis=0)
# x_std = X_train.std(axis=0)

X_mean = X.mean(axis=0)
X_std = X.std(axis=0)
Y_mean = Y.mean(axis=0)
Y_std = Y.std(axis=0)

# # print(f_std)
# # Normalization step
# Y_eval = (Y_eval - l_mean) / l_std
# Y_labeled = (Y_labeled - Y_labeled_mean) / Y_labeled_std
# X_train = (X_train - x_mean) / x_std
# X_labeled = (X_labeled - x_labeled_mean) / x_labeled_std

Y = (Y - Y_mean) / Y_std
X = (X - X_mean) / X_std

X_train = X[int(ll*len(X)):,:]
X_labeled = X[:int(ll*len(X)),:]
Y_labeled = Y[:int(ll*len(X)),:]
Y_eval = Y[int(ll*len(X)):,:]

# a = np.array([[[3,4,5],[1,5,3]],[[6,4,5],[8,6,7]]])
# logm(a)
# np.nonzero(features==0)
# normalized_angle.max()
# R_Cov_BS.min()
# print(f_std)
# R_Cov_BS[0,:]
# upper_triangles[:10,0].imag
# features.max()
# plt.scatter(X[:,0],X[:,1])
# plt.show()
# features.max()
# print(X[0,:])
# print(X[1,:])
# X.max()
# print(features.std())
# np.count_nonzero(np.isnan(features))
# np.count_nonzero(features==0)
# np.where(features == 0)[0]
# R_Cov_BS[3,:]
# print((X[0,:]))

# print(np.log(features[0,:]))
# Y.shape
# print(X)


# class_weights = compute_class_weight('balanced', classes=np.unique(X), y=X)
# class_weight_dict = dict(enumerate(class_weights))


# Step 3: Build the Supervised Autoencoder
input_dim = X.shape[1]  # Input feature dimension
latent_dim = 2          # Dimension of latent space (same as Y)
seed = 42

# Encoder
inputs = tf.keras.Input(shape=(input_dim,))
# inputs = tf.keras.Input(shape=(input_dim,))
# noisy_inputs = layers.GaussianNoise(0.05, seed=seed)(inputs)  # Add noise
# encoder part
# Bs 1
Xencoded = Dense(1024, activation='relu')(inputs)
Xencoded = Dense(512, activation='relu')(Xencoded)
Xencoded = Dense(256, activation='tanh')(Xencoded)
Xencoded = Dense(128, activation='tanh')(Xencoded)
Xencoded = Dense(64, activation='tanh')(Xencoded) # 'tanh'
latent_space = layers.Dense(latent_dim, activation='linear' # , kernel_regularizer=tf.keras.regularizers.l2(0.5)
                            )(Xencoded)
# latent_space = Dense(latent_dim, activation='linear', 
#                      kernel_regularizer=tf.keras.regularizers.l2(0.001))(Xencoded)

# Decoder
# Decoder part
# Bs 1
Xdecoded = Dense(64, activation='tanh')(latent_space) #32,64,128,256,1024
Xdecoded = Dense(128, activation='tanh')(Xdecoded)    
Xdecoded = Dense(256, activation='tanh')(Xdecoded)
Xdecoded = Dense(512, activation='relu')(Xdecoded)
Xdecoded = Dense(1024, activation='relu')(Xdecoded)
reconstructed = layers.Dense(input_dim, activation='linear', name='reconstructed')(Xdecoded)

##########
# # Supervised branch (predict target labels)
# def custom_loss(y_true, y_pred):
#     return tf.reduce_mean(tf.square(y_true - y_pred)) + tf.reduce_mean(tf.abs(y_true - y_pred))

# # Define the model
# autoencoder = models.Model(inputs, [reconstructed], name="supervised_autoencoder")
# autoencoder.compile(optimizer=Adam(learning_rate=0.0001), 
#                     loss={'reconstructed': 'custom_loss'},
#                     loss_weights={'reconstructed': 1})

# # encoder = models.Model(inputs, [latent_space], name="supervised_latent_space")
# decoder = models.Model(latent_space, reconstructed, name="supervised_latent_space")

# decoder.compile(optimizer=Adam(learning_rate=0.01), 
#                     loss='mse')

# # Extract the encoder part for dimensionality reduction
# # encoder = models.Model(inputs, latent_space, name="encoder")

# # Add early stopping
# # from tensorflow.keras.callbacks import EarlyStopping
# # early_stopping = EarlyStopping(monitor='val_loss', patience=20, restore_best_weights=True)

# history = decoder.fit(Y_labeled, {'reconstructed': X_labeled}, 
#                           batch_size=1024, epochs=300, verbose=1)
# # Step 4: Train the Supervised Autoencoder
# # history = encoder.fit(X_labeled, {'latent_space': Y_labeled}, 
# #                           batch_size=1024, epochs=300, verbose=1)
#                           #batch_size=1024, epochs=300, verbose=1, validation_split=0.4, callbacks=[early_stopping])
#####################
autoencoder = tf.keras.Model(inputs=inputs, outputs=reconstructed)

# # Compile the model
# autoencoder.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), 
#                     loss='mse')

# Optimizer
optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
encoder_optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
decoder_optimizer = tf.keras.optimizers.Adam(learning_rate=0.05)

# # Define custom loss to handle semi-supervised learning
# def semi_supervised_loss(X_true, X_pred, X_anchor, latent_output, alpha=0.5):
#     # MSE loss for reconstruction
#     reconstruction_loss = tf.reduce_mean(tf.square(X_true - X_pred))
    
#     # Supervised loss for anchors (minimize distance between labeled latent space)
#     latent_anchor_output = encoder(X_anchor)
#     supervised_loss = tf.reduce_mean(tf.square(latent_anchor_output - X_anchor[:,:latent_dim]))
    
#     # Combined loss
#     return alpha * reconstruction_loss + (1 - alpha) * supervised_loss

# Custom training loop

epochs = 200 #len(X)/160
batch_size = 1024 # int(len(X) / epochs) #1024
alpha = .5  # encoder weights 0.9
beta = .25    # autoencoder weights 0.5
gamma = .25   # decoder weights 0.5

encoder = tf.keras.Model(inputs=autoencoder.input, outputs=latent_space)
decoder = tf.keras.Model(inputs=latent_space,outputs=autoencoder.output)

for epoch in range(epochs):
    print(f"Epoch {epoch+1}/{epochs}")
    for step in range(len(X) // batch_size):
        # Get a batch of unlabeled data
        idx = np.random.choice(len(X), batch_size, replace=False)
        X_batch = X[idx]

        with tf.GradientTape(persistent=True) as tape:
            # Forward pass for reconstruction
            reconstruction = autoencoder(X_batch, training=True)
            reconstruction_loss = tf.reduce_mean(tf.square(X_batch - reconstruction))

            # Forward pass for anchor data-encoder part (latent representation)
            latent_labeled_output = encoder(X_labeled)
            anchor_labels = Y_labeled  # Replace with true anchor labels if available
            anchor_loss_enc = tf.reduce_mean(tf.square(latent_labeled_output - anchor_labels))
            
            # Forward pass for anchor data-decoder part (feature reconstruction)
            feature_out = decoder(Y_labeled)
            anchor_loss_dec = tf.reduce_mean(tf.square(feature_out - X_labeled))
            anchor_loss = alpha*anchor_loss_enc+gamma*anchor_loss_dec

            
            # Combine losses
            total_loss = beta*reconstruction_loss + anchor_loss

        # Backpropagation for whole AE
        gradients = tape.gradient(total_loss, autoencoder.trainable_variables)
        optimizer.apply_gradients(zip(gradients, autoencoder.trainable_variables))

        # Compute gradients specifically for the encoder part
        encoder_gradients = tape.gradient(anchor_loss_enc, encoder.trainable_variables)
        encoder_optimizer.apply_gradients(zip(encoder_gradients, encoder.trainable_variables))

        # Compute gradients specifically for the decoder part
        decoder_gradients = tape.gradient(anchor_loss_dec, decoder.trainable_variables)
        decoder_optimizer.apply_gradients(zip(decoder_gradients, decoder.trainable_variables))
#         del tape
        print(f"Rec Loss: {reconstruction_loss.numpy()}, "
          f"Loss (Encoder): {anchor_loss_enc.numpy()}, "
          f"Loss (Decoder): {anchor_loss_dec.numpy()}")

    from sklearn.metrics import mean_squared_error

# Step 1: Compute MSE for Training Data
X_reduced_train = encoder.predict(X)  # Predicted latent representations
X_reduced_train = X_reduced_train*(l_std)+l_mean   #*(l_max-l_min)+l_min

Y = Y *(l_std)+l_mean
# X_reduced_train = X_reduced_train*(l_max-l_min)+l_min
# X_reduced_train  = np.expm1(X_reduced_train*(l_max-l_min)+l_mean)
# Y= np.expm1(Y *(l_max-l_min)+l_mean)

# X_reduced_train = (X_reduced_train-l_min)/(l_max-l_min)
# Y = Y *(l_max-l_min)+l_min

# mse = mean_squared_error(Y, X_reduced_train)  # Compute MSE between predicted and true labels
z1 = np.sqrt((X_reduced_train[:,0]-Y[:,0])**2
                    +(X_reduced_train[:,1]-Y[:,1])**2)
MLE = np.sum(z1)/len(z1)
# Step 2: Display the MSE
print(f"Mean localization Error (MLE) between predicted and ground truth labels: {MLE}")

# Optional: Visualize Predicted vs Ground Truth
# plt.figure(figsize=(8, 6))
# plt.scatter(Y[:, 0], Y[:, 1], label='Ground Truth', alpha=0.7)
# plt.scatter(X_reduced_train[:, 0], X_reduced_train[:, 1], label='Predicted (Encoded)', alpha=0.7)
# plt.title("Comparison of Ground Truth and Predicted Latent Space")
# plt.xlabel("Latent Dimension 1")
# plt.ylabel("Latent Dimension 2")
# plt.legend()
# plt.show()

io.savemat('ae_output.mat',dict(R_xy=Y,xy_FCAE=X_reduced_train))
