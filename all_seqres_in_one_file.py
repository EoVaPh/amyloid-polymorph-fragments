from pathlib import Path


input_folder = Path('extracted_chains')
output_folder = Path(r"C:\Users\User\Documents\bioinf\smtb")

output_folder.mkdir(exist_ok=True)

output_file = output_folder / "all_seqres_from_cifs.txt"


with open(output_file, "w", encoding="utf8") as output:

    for file in input_folder.glob("*_seq.txt"):
        with open(file, "r", encoding="utf8") as input_file:
            lines = input_file.readlines()

        sequence = lines[0].strip()
        file_name = file.stem[:4]
