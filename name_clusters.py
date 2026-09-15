from pathlib import Path


clusters_file = Path("clusters_final.txt")
names_file = Path("clusters_names.txt")
output_file = Path("clusters_renamed.txt")


cluster_names = {}

with open(names_file, "r", encoding="utf-8") as file:
    for line in file:
        line = line.strip()

        if not line:
            continue

        parts = line.split()

        if len(parts) < 2:
            continue

        cluster_id = parts[0]
        cluster_name = parts[1]
        cluster_names[cluster_id] = cluster_name

output_lines = []

with open(clusters_file, "r", encoding="utf-8") as file:
    for line in file:
        stripped = line.strip()

        if stripped.startswith(">Cluster"):
            cluster_number = stripped.replace(">Cluster", "").strip()
            cluster_id = f"Cluster_{cluster_number}"

            if cluster_id in cluster_names:
                new_name = cluster_names[cluster_id]
                output_lines.append(f">{new_name}\n")
            else:
                print(f"WARNING: {cluster_id} not found in name mapping")
                output_lines.append(line)

        elif stripped:
            pdb_id = stripped.split()[0]
            output_lines.append(f"{pdb_id}\n")

with open(output_file, "w", encoding="utf-8") as file:
    file.writelines(output_lines)

print(f"Done. Renamed clusters saved to: {output_file}")