import argparse
from api_tensorflow import audioProcess
from input_lpc_output_weight import WeightsAnimation
from toarkit import main
import numpy as np
"""
python main.py wav/xuanyi.wav res.csv
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "wav_path"
    )
    parser.add_argument(
        "output_path"
    )
    args = parser.parse_args()

    pb_path = "./best_model/Audio2Face"
    bs_name = np.load("./shirley_1119_bs_name.npy")
    audio_data = audioProcess(args.wav_path)
    weight_data = WeightsAnimation('./best_model/Audio2Face.tflite', pb_path).get_weight(audio_data)
    arkit_df = main(weight_data, bs_name)
    arkit_df.to_csv(args.output_path, index = True) 