import os
import zipfile
from io import BytesIO
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from Bio import SeqIO

# Define paths
ZIP_FILE = "data.zip"
EXTRACT_FOLDER = "data"
os.makedirs(EXTRACT_FOLDER, exist_ok=True)

def unzip_file(file_path, output_folder):
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(output_folder)
    return os.listdir(output_folder)

def encode_sequence(sequence, max_length):
    """Encodes a DNA sequence into numeric format."""
    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    encoded = [mapping.get(base, 0) for base in sequence.upper()]  # Map bases to integers
    # Pad sequences to the same length
    return tf.keras.preprocessing.sequence.pad_sequences([encoded], maxlen=max_length, padding='post')[0]

def load_data(data_folder, max_length):
    sequences = []
    labels = []

    for file in os.listdir(data_folder):
        label = 1 if "diseased" in file.lower() else 0  # Label: 1 = Diseased, 0 = Healthy
        file_path = os.path.join(data_folder, file)
        sequence = str(SeqIO.read(file_path, "fasta").seq)
        encoded_sequence = encode_sequence(sequence, max_length)
        sequences.append(encoded_sequence)
        labels.append(label)
    
    return np.array(sequences), np.array(labels)

unzipped_files = unzip_file(ZIP_FILE, EXTRACT_FOLDER)
MAX_LENGTH = 500  # Set maximum sequence length
X, y = load_data(EXTRACT_FOLDER, MAX_LENGTH)
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
l
model = tf.keras.Sequential([
    tf.keras.layers.Embedding(input_dim=5, output_dim=16, input_length=MAX_LENGTH),  # Embedding for DNA bases
    tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences=True)),
    tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(32)),
    tf.keras.layers.Dense(32, activation='relu'),
    tf.keras.layers.Dropout(0.5),
    tf.keras.layers.Dense(1, activation='sigmoid')  # Output: probability of disease
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

history = model.fit(
    X_train, y_train,
    validation_data=(X_val, y_val),
    epochs=10,  # Adjust the number of epochs as needed
    batch_size=32
)

MODEL_PATH = "my_model.h5"
model.save(MODEL_PATH)
print(f"Model saved to {MODEL_PATH}")
