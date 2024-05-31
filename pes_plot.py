#%%
import json
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

with open('data/H2_vqe_sto3g_PES.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    data_np = np.array(list(data.values()))
    x = np.array(list(data.keys()), dtype=float)
    for i in range(data_np.shape[1]):
        a = ax.scatter(x, data_np[:, i], color='red')
    a.set_label('VQE (UCCSD)')

# with open('data/H2_vqe_sto3g_hf_PES.json', 'r', encoding='utf-8') as f:
#     data = json.load(f)
#     data_np = np.array(list(data.values()))
#     x = np.array(list(data.keys()), dtype=float)
#     b = ax.plot(x, data_np, color='black')
#     b[0].set_label('VQE (HF)')

with open('data/H2_exact_PES.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    data_np = np.array(list(data.values()))
    x = np.array(list(data.keys()), dtype=float)
    c = ax.plot(x, data_np, color='blue')
    c[0].set_label('Exact')

plt.xlabel('Bond Length')
plt.ylabel('Energy')
plt.legend()
plt.savefig('figures/H2_total_PES.png')