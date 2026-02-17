import os
import numpy as np
import pickle
import uuid
import shutil
import gzip
import tarfile
import trimUtils
import sys

if len(sys.argv) < 2:
  print("Error! Pass in config!\n")
  sys.exit()

inpFilename = sys.argv[1]
mass_dict = trimUtils.createMassDict("mass_1.mas20.txt")
materials_dict = trimUtils.createMaterialsDict()
input_data = trimUtils.parseConfig(inpFilename)
trimUtils.checkTrimArgs(input_data,mass_dict,materials_dict)

SRIM_EXE_PATH= "/Users/patrickstengel/Documents/SRIM-2013"
SRIM_TMP_PATH = "/Users/patrickstengel/Documents/TRIM_tracks_CCs/trimUtils/trimRunner/SRIM_tmp"

target_name = input_data["material"]
ion_name = input_data["ionSymbol"]
outputFolder = input_data["outputPath"]
runMode = input_data["runMode"]
energy = input_data["energy_keV"]
nps = input_data["nps"]
seed = input_data["seed"]

output_base_name_coll = f"COLLISON_{target_name}_{ion_name}_seed_{seed}"
output_base_name_exyz = f"EXYZ_{target_name}_{ion_name}_seed_{seed}"

# Create new SRIM run folder, go there
runFolder = trimUtils.makeTempSRIMFolder(f"{target_name}_{ion_name}",SRIM_TMP_PATH,SRIM_EXE_PATH)
try:
  os.chdir(runFolder)

  if runMode=="damage":
    txt_path_coll = os.path.join(outputFolder, f"{output_base_name_coll}.txt")
    tar_path_coll = os.path.join(outputFolder, f"{output_base_name_coll}.tar.gz")
    if os.path.exists(txt_path_coll):
      os.remove(txt_path_coll)
    if os.path.exists(tar_path_coll):
      os.remove(tar_path_coll)

    txt_path_exyz = os.path.join(outputFolder, f"{output_base_name_exyz}.txt")
    tar_path_exyz = os.path.join(outputFolder, f"{output_base_name_exyz}.tar.gz")
    if os.path.exists(txt_path_exyz):
      os.remove(txt_path_exyz)
    if os.path.exists(tar_path_exyz):
      os.remove(tar_path_exyz)

    # Run SRIM return location of collison.txt [sic] 
    collison_path, exyz_path = trimUtils.runSRIM(runFolder,input_data,mass_dict,materials_dict)

    # Move and compress output
    os.replace(collison_path, txt_path_coll)
    
    os.replace(exyz_path, txt_path_exyz)

    # Create .tar.gz files
    with tarfile.open(tar_path_coll, "w:gz") as tar:
      tar.add(txt_path_coll, arcname=os.path.basename(txt_path_coll))
    os.remove(txt_path_coll)

    with tarfile.open(tar_path_exyz, "w:gz") as tar:
      tar.add(txt_path_exyz, arcname=os.path.basename(txt_path_exyz))
    os.remove(txt_path_exyz)

finally:
  os.chdir(SRIM_TMP_PATH)
  shutil.rmtree(runFolder, ignore_errors=True)
