#%%
import json
import numpy as np
import matplotlib.pyplot as plt

def plot_scatter(ax, data, label, color):
    data_np = np.array(list(data.values()))
    x = np.array(list(data.keys()), dtype=float)
    for i in range(data_np.shape[1]):
        a = ax.scatter(x, data_np[:, i], color=color)
    a.set_label(label)

def plot_line(ax, data, label, color):
    data_np = np.array(list(data.values()))
    x = np.array(list(data.keys()), dtype=float)
    b = ax.plot(x, data_np, color=color)
    b[0].set_label(label)

fig, ax = plt.subplots()
plt.xlabel('Bond Length')
plt.ylabel('Energy')

with open('data/H2_sfuccsd_sto3g_PES.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    plot_scatter(ax, data, 'SSVQE (UCCSD)', 'red')

with open('data/H2_exact_PES.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    plot_line(ax, data, 'Exact', 'blue')

# with open('data/H2_uccsd_hf_sto3g_PES.json', 'r', encoding='utf-8') as f:
#     data = json.load(f)
#     plot_line(ax, data, 'SSVQE (HF)', 'black')

# with open('4-31g_test/H2_sfuccsd_sto3g_PES.json', 'r', encoding='utf-8') as f:
#     data = json.load(f)
#     plot_scatter(ax, data, 'SSVQE (SFUCCSD)', 'green')

plt.legend()
plt.savefig('figures/H2_total_modified_PES.png')

#%%

import json
import numpy as np
import matplotlib.pyplot as plt

def plot_scatter(ax, data, label, color):
    data_np = np.array(list(data.values()))
    x = np.array(list(data.keys()), dtype=float)
    for i in range(data_np.shape[1]):
        a = ax.scatter(x, data_np[:, i], color=color)
    a.set_label(label)

def plot_line(ax, data, label, color):
    data_np = np.array(list(data.values()))
    x = np.array(list(data.keys()), dtype=float)
    b = ax.plot(x, data_np, color=color)
    b[0].set_label(label)

fig, ax = plt.subplots()
plt.xlabel('Bond Length')
plt.ylabel('Energy')

# with open('data/H2_sfuccsd_sto3g_PES.json', 'r', encoding='utf-8') as f:
#     data = json.load(f)
#     plot_scatter(ax, data, 'SSVQE (UCCSD)', 'red')

with open('4-31g_test/H2_exact_PES.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    plot_line(ax, data, 'Exact', 'blue')

with open('4-31g_test/H2_sfuccsd_sto3g_PES.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    plot_scatter(ax, data, 'SSVQE (SFUCCSD)', 'green')

plt.legend()
plt.savefig('4-31g_test/H2_total_PES.png')
