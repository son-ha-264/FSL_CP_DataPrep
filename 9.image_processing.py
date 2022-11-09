import pandas as pd
import os
import numpy as np
from PIL import Image
from tqdm import tqdm

# Helper Functions
def list_tif_to_multichannel_npz(dict_tif_path, npz_save_path):
    """Open 5 tif images, perform image processing, convert to npz, and save to 'npz_save_path'
    """
    list_npz = []
    dict_filenames = {}
    list_dyes = ['ERSyto', 'ERSytoBleed', 'Hoechst', 'Mito', 'Ph_golgi']

    for dye in list_dyes:
        tif_path = dict_tif_path[dye]

        # Save original .tiff name 
        tif_name = tif_path.split("/")[-1]
        dict_filenames[dye] = tif_name
        
        # Convert to npz and save
        img = Image.open(tif_path)
        img_arr = np.array(img)
        # Conv img to 8bit and remove 0.0028% outlier bits
        list_npz.append((np.minimum(img_arr / np.percentile(img_arr, q=99.9972), 1.0) * 255).astype(np.uint8))

    final_array = np.stack(list_npz, axis=2)
    np.savez(npz_save_path, sample=final_array, filenames=dict_filenames)
    return None


# Init paths
base_image_path = "/home/son.ha/gigascience_30k/images_unzipped"
final_image_path = "/mnt/scratch/Son_cellpainting/my_cp_images"

# Other inits
df = pd.read_csv("temp/cp_inchi_smiles_sample_key_map.csv")
plate_list = list(set(df["Metadata_Plate"]))
sample_key_plate = list(set(df["SAMPLE_KEY"]))
dict_sample_key_and_paths = {}


### Program start
# Loop through all plates, group all images of the same plate-view but different dyes 
for plate in tqdm(plate_list, desc= 'Iterate through all plates'):
    """
    list_dyes_in_plate = [
        str(plate)+"-ERSyto",
        str(plate)+"-ERSytoBleed",
        str(plate)+"-Hoechst",
        str(plate)+"-Mito",
        str(plate)+"-Ph_golgi", # This exact order for dyes
    ]
    """
    list_dyes = ['ERSyto', 'ERSytoBleed', 'Hoechst', 'Mito', 'Ph_golgi'] # This exact order for dyes
    for dye in list_dyes:
        plate_dye = str(plate)+'-'+dye
        current_folder = os.path.join(base_image_path, plate_dye)
        if not os.path.isdir(current_folder):
            continue
        for filename in os.listdir(current_folder):
            if filename.endswith(".tif"):
                well = filename.split("_")[1].upper()
                view = filename.split("_")[2][1]
                sample_key = str(plate)+"-"+well
                if sample_key in sample_key_plate:
                    image_dir = os.path.join(current_folder, filename)
                    sample_key_with_view = sample_key+'-'+view
                    if sample_key_with_view in dict_sample_key_and_paths.keys():
                        dict_sample_key_and_paths[sample_key_with_view][dye] = image_dir
                    else:
                        dict_sample_key_and_paths[sample_key_with_view] = {dye:image_dir}

# Save each view into a 5-channel .npz image (because there are 5 dyes for each view)
for sample_key_with_view in tqdm(dict_sample_key_and_paths.keys(), desc= 'Saving each view as .npz'):
    list_tif_to_multichannel_npz(
        dict_tif_path=dict_sample_key_and_paths[sample_key_with_view],
        npz_save_path=os.path.join(final_image_path, sample_key_with_view+'.npz')
    )
