with open("data/srs_nlos_los.csv", "r") as f:
    data = f.readlines()
    
data = data[1:]  # Skip header

srs_nlos_los = []
for dat in data:
    dat.strip()
    s_data = dat.split(",")
    srs_nlos_los.append([float(s_data[0]), float(s_data[1]), float(s_data[2])])
    
with open("data/prs_nlos_los.csv", "r") as f:
    data = f.readlines()
    
data = data[1:]  # Skip header

prs_nlos_los = []
for dat in data:
    dat.strip()
    s_data = dat.split(",")
    prs_nlos_los.append([float(s_data[0]), float(s_data[1]), float(s_data[2])])
    
base_counts = []
with open("data/seq_09_gps_with_site_count_voronoi.csv", "r") as f:
    data = f.readlines()


data = data[1:]  # Skip header
for dat in data:
    dat.strip()
    s_data = dat.split(",")
    base_counts.append(int(s_data[5]))

#print how many unique items there are and how many of each unique item there are
unique_base_counts = set(base_counts)
print(f"Unique base station counts: {len(unique_base_counts)}")
# Print the counts of each unique item
from collections import Counter
base_count_counter = Counter(base_counts)
print("Base station counts:")
for count, freq in base_count_counter.items():
    print(f"{count}: {freq} times")
    
with open("data/srs_nlos_los_new.csv", "w") as f:
    f.write("x,y,z,base_stations\n")
    for idx,row in enumerate(srs_nlos_los):
        idx = idx + 1
        f.write(f"{row[0]},{row[1]},{row[2]},{base_counts[len(base_counts)-idx]}\n")
        
with open("data/prs_nlos_los_new.csv", "w") as f:
    f.write("x,y,z,base_stations\n")
    for idx,row in enumerate(prs_nlos_los):
        idx = idx + 1
        f.write(f"{row[0]},{row[1]},{row[2]},{base_counts[len(base_counts)-idx]}\n")