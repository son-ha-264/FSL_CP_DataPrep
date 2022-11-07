import pandas as pd
import os
import numpy as np
from PIL import Image 

def list_tif_to_multichannel_npz(list_tif_path, npz_save_path):
    """Open 5 tif images, perform image processing, convert to npz, and save to 'npz_save_path'
    """
    list_npz = []
    for tif_path in list_tif_path:
        img = Image.open(tif_path)
        img_arr = np.array(img)
        # Conv img to 8bit and remove 0.0028% outlier bits
        list_npz.append((np.minimum(img_arr / np.percentile(img_arr, q=99.9972), 1.0) * 255).astype(np.uint8))

    # Save .npz

    return None



# Init paths
base_image_path = "/home/son.ha/gigascience_30k/images_toy"

# Other inits
df = pd.read_csv("temp/cp_inchi_smiles_sample_key_map.csv")
plate_list = list(set(df["Metadata_Plate"]))
sample_key_plate = list(set(df["SAMPLE_KEY"]))

# Loop through all plates
for plate in plate_list:
    list_dyes_in_plate = [
        str(plate)+"-ERSyto",
        str(plate)+"-ERSytoBleed",
        str(plate)+"-Hoechst",
        str(plate)+"-Mito",
        str(plate)+"-Ph_golgi",
    ]
    for plate_dye in list_dyes_in_plate:
        current_folder = os.path.join(base_image_path, plate_dye)
        if not os.path.isdir(current_folder):
            continue
        for filename in os.listdir(current_folder):
            if filename.endswith(".tif"):
                well = filename.split("_")[1].upper()
                view = filename.split("_")[2][1]
                sample_key = str(plate)+"-"+well
                if sample_key in sample_key_plate:
                    #print(sample_key)
                    #remember to sort the tif_list of 5 dyes so that the dyes are in same order
                    break