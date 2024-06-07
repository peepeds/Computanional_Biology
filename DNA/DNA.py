# %%
# %%

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import numpy

# %%
dna_1 = Seq('ATGATCTCGTAACAGGTAACAAACC')
dna_2 = Seq('ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA')

# %% [markdown]
# Nomor 1

# %%

len_dna_1 = len(dna_1)
len_dna_2 = len(dna_2)
print("Panjang DNA 1 : ", len_dna_1)
print("Panjang DNA 2 : ", len_dna_2)

# %% [markdown]
# Nomor 2

# %%
def count_freq(dna) : 
    freq_A = dna.count('A')
    freq_T = dna.count('T')
    freq_G = dna.count('G')
    freq_C = dna.count('C')
    freq = {
        'A' : freq_A,
        'T' : freq_T,
        'G' : freq_G,
        'C' : freq_C
    }
    return freq

# %%
freq_DNA_1 = count_freq(dna_1)
print("Frekuensi basa dna 1 : \n",freq_DNA_1)

# %%
freq_DNA_2 = count_freq(dna_2)
print("Frekuensi basa dna 2 : \n",freq_DNA_2)

# %%
# Mengatur lebar batang
bar_width = 0.15

# Posisi batang pada sumbu x
indices = numpy.arange(len(freq_DNA_1))

# Membuat grafik batang berdampingan
plt.bar(indices - bar_width, list(freq_DNA_1.values()), bar_width, label="DNA 1")
plt.bar(indices + bar_width, list(freq_DNA_2.values()), bar_width, label="DNA 2")

# Menambahkan label dan judul
plt.xlabel('basa')
plt.ylabel('freq')
plt.xticks(indices, list(freq_DNA_1.keys()))
plt.legend()

# Menampilkan grafik
plt.show()

# %% [markdown]
# Nomor 3

# %%
def content(dna):
    at_content = float(dna.count('A') + dna.count('T')) / len(dna) * 100
    gc_content = 100 - at_content
    return at_content , gc_content


# %%
at_content_dna1 , gc_content_dna1 = content(dna_1)
at_content_dna2 , gc_content_dna2 = content(dna_2)

print(f"AT Content DNA 1 :{at_content_dna1:.2f}%")
print(f"GC Content DNA 1 :{gc_content_dna1:.2f}%")
print(f"AT Content DNA 2 :{at_content_dna2:.2f}%")
print(f"GC Content DNA 2 :{gc_content_dna2:.2f}%")


# %%
def melting(dna):
    Wallace = mt.Tm_Wallace(dna)
    GC = mt.Tm_GC(dna)
    NN = mt.Tm_NN(dna)
    melting_point = {
        'Wallace' : str(Wallace)+'°C',
        'GC' : str(GC) + '°C' ,
        'NN' : str(NN) + '°C'
    }
    return melting_point


# %%
melting_dna_1 = melting(dna_1)
melting_dna_2 = melting(dna_2)
print("Melting point DNA 1 :\n",melting_dna_1)
print("Melting point DNA 2 :\n",melting_dna_2)

# %% [markdown]
# Nomor 4

# %%
def gc_skew(dna, window=2):
    skews = []
    for i in range(0, len(dna), window):
        subseq = dna[i:i+window]
        g_count = subseq.count('G')
        c_count = subseq.count('C')
        if g_count + c_count !=0:
            skew = (g_count - c_count)/(g_count + c_count) 
            skews.append(skew)
        else:
            skews.append(0)

    return skews

# %%
skew1 = gc_skew(dna_1)
skew2 = gc_skew(dna_2)

# %%
plt.plot(skew1, label='generic_dna1' , color = 'red')
plt.plot(skew2, label='generic_dna2' , color = 'blue')

plt.xlabel('Posisi')
plt.ylabel('GC skewness')
plt.title('GC Skewness dengan windows 2')

plt.legend()
plt.show()

# %% [markdown]
# Nomor 5

# %%
global_alignments = pairwise2.align.globalxx(dna_1, dna_2) 
local_alignments = pairwise2.align.localxx(dna_1, dna_2)
print('Global alignment:\n') 
print(format_alignment(*global_alignments[0]))
print('Local alignment:\n') 
print(format_alignment(*local_alignments[0]))


