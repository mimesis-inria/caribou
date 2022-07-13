import subprocess
import os

material_path = "./FEniCS_code_generation"
material_directories = []
for root, dirnames, filenames in os.walk(material_path):
    if not len(dirnames) == 0:
        material_directories.append(dirnames)
    else:
        assert "No material model in this directory"
for material in material_directories[0]:
    path = os.path.join(os.path.dirname(__file__), material_path + "/" + material)
    os.chdir(path)
    subprocess.Popen(["ffcx", material + "_Tetra.py", "-o", "../../FEniCS_Generated_code"])
    subprocess.Popen(["ffcx", material + "_Tetra_Order2.py", "-o", "../../FEniCS_Generated_code"])
    subprocess.Popen(["ffcx", material + "_Hexa.py", "-o", "../../FEniCS_Generated_code"])
    subprocess.Popen(["ffcx", material + "_Hexa_Order2.py", "-o", "../../FEniCS_Generated_code"])
