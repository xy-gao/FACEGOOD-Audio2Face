import numpy as np
import pandas as pd
import argparse

mapping_one = {
    'jawForward': 'jaw_thrust_c',
    'jawLeft': 'jaw_sideways_l',
    'jawRight': 'jaw_sideways_r',
    'jawOpen': 'mouth_stretch_c',
    'mouthClose': 'mouth_chew_c',
    'mouthLeft': 'mouth_sideways_l',
    'mouthRight': 'mouth_sideways_r',
    'mouthSmileLeft': 'mouth_lipCornerPull_l',
    'mouthSmileRight': 'mouth_lipCornerPull_r',
    'mouthDimpleLeft': 'mouth_dimple_l',
    'mouthDimpleRight': 'mouth_dimple_r',
    'mouthStretchLeft': 'mouth_lipStretch_l',
    'mouthStretchRight': 'mouth_lipStretch_r',
    'mouthShrugLower': 'mouth_chinRaise_d',
    'mouthShrugUpper': 'mouth_chinRaise_u',
    'mouthPressLeft': 'mouth_press_l',
    'mouthPressRight': 'mouth_press_r',
    'mouthLowerDownLeft': 'mouth_lowerLipDepress_l',
    'mouthLowerDownRight': 'mouth_lowerLipDepress_r',
    'mouthUpperUpLeft': 'mouth_upperLipRaise_l',
    'mouthUpperUpRight': 'mouth_upperLipRaise_r'
}

mapping_max = {
    "mouthFunnel": ["mouth_funnel_dl", "mouth_funnel_dr", "mouth_funnel_ul", "mouth_funnel_ur"],
    "mouthPucker": ["mouth_pucker_l", "mouth_pucker_r"],
    "mouthFrownLeft": ["mouth_lipCornerDepress_l", "mouth_lipCornerDepressFix_l"],
    "mouthFrownRight": ["mouth_lipCornerDepress_r", "mouth_lipCornerDepressFix_r"],
    "mouthRollLower": ["mouth_suck_dl", "mouth_suck_dr"],
    "mouthRollUpper": ["mouth_suck_ul", "mouth_suck_ur"],
}

ARKIT_names = [
    "eyeBlinkLeft",
    "eyeLookDownLeft",
    "eyeLookInLeft",
    "eyeLookOutLeft",
    "eyeLookUpLeft",
    "eyeSquintLeft",
    "eyeWideLeft",
    "eyeBlinkRight",
    "eyeLookDownRight",
    "eyeLookInRight",
    "eyeLookOutRight",
    "eyeLookUpRight",
    "eyeSquintRight",
    "eyeWideRight",
    "jawForward",
    "jawLeft",
    "jawRight",
    "jawOpen",
    "mouthClose",
    "mouthFunnel",
    "mouthPucker",
    "mouthLeft",
    "mouthRight",
    "mouthSmileLeft",
    "mouthSmileRight",
    "mouthFrownLeft",
    "mouthFrownRight",
    "mouthDimpleLeft",
    "mouthDimpleRight",
    "mouthStretchLeft",
    "mouthStretchRight",
    "mouthRollLower",
    "mouthRollUpper",
    "mouthShrugLower",
    "mouthShrugUpper",
    "mouthPressLeft",
    "mouthPressRight",
    "mouthLowerDownLeft",
    "mouthLowerDownRight",
    "mouthUpperUpLeft",
    "mouthUpperUpRight",
    "browDownleft",
    "browDownRight",
    "browInnerUp",
    "browOuterUpLeft",
    "browOuterUpRight",
    "cheekPuff",
    "cheekSquintLeft",
    "cheekSquintRight",
    "noseSneerLeft",
    "noseSneerRight"
]


def load_bs_to_df(bs_value, bs_name):
    weights = np.zeros((bs_value.shape[0],len(bs_name)))
    bs_name_index = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 94, 93, 95, 96, 97, 98, 99, 100, 101, 102, 103, 105, 104, 106, 107, 108, 109, 110, 111, 112, 113, 114, 1, 115]
    for i in range(len(bs_name_index)):
        weights[:,i] = bs_value[:,bs_name_index[i]]
    df = pd.DataFrame(weights, columns=bs_name)
    return df

def convert_to_arkit_df(df, ARKIT_names, mapping_one, mapping_max):
    weights_zeros = np.zeros((df.shape[0], len(ARKIT_names)))
    df_arkit = pd.DataFrame(weights_zeros, columns=ARKIT_names)
    
    for k, v in mapping_one.items():
        df_arkit[k] = df[v]
    
    for k, v in mapping_max.items():
        df_arkit[k] = df[v].max(axis=1)
    return df_arkit


def main(bs_value, bs_name):
    BS_CONUNT = 116
    const_bs_value = [0.,0.,-0.,0.,-0.,0.,-0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.,-0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.,0.,-0.,0.,-0.,0.,-0.,0.,0.,-0.,0.,-0.,0.,-0.,0.,-0.,0.,-0.,0.]
    var_bs_index = [10, 13, 14, 15, 18, 33, 38, 40, 41, 42, 43, 44, 45, 51, 52, 53, 54, 57, 58, 59, 60, 63, 64, 66, 67, 68, 69, 74, 75, 76, 77, 78, 79, 80, 81, 82, 84]
    const_bs_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 39, 46, 47, 48, 49, 50, 55, 56, 61, 62, 65, 70, 71, 72, 73, 83, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115]

    weights = np.zeros((bs_value.shape[0],BS_CONUNT))
    # print(weights.shape)

    weights[:,var_bs_index] = bs_value
    weights[:,const_bs_index] = const_bs_value

    

    df_src = load_bs_to_df(weights, bs_name)
    df_arkit = convert_to_arkit_df(df_src, ARKIT_names, mapping_one, mapping_max)

    return df_arkit
    # plot converted bs values
    #df_arkit.iloc[:400, :].plot(figsize=(20, 4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="inference")

    parser.add_argument(
        "--bs_value",
        default="./model_output/xuanyi_output.npy",
        type=str,
        metavar="PATH",
        help="path to FACEGOOD bs value in npy",
    )

    parser.add_argument(
        "--bs_name",
        default="./shirley_1119_bs_name.npy",
        type=str,
        metavar="PATH",
        help="path to FACEGOOD bs_name in npy",
    )

    args = parser.parse_args()
    bs_value = np.load(args.bs_value)
    bs_name = np.load(args.bs_name)
    res = main(bs_value, bs_name)
    res.to_csv('xuanyi_res.csv', index = True) 