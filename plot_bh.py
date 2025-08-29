import pandas as pd
import matplotlib.pyplot as plt

# Load data from file
data = pd.read_csv('bh_data.csv')

# Creating a graph
plt.figure(figsize=(10, 6))
plt.plot(data['H'], data['B'], label='B-H curve', color='blue')
plt.xlabel('Magnetic Field Strength (H), A/m')
plt.ylabel('Magnetic Induction (B), T')
plt.title('B-H Hysteresis Loop')
plt.grid(True)
plt.legend()
plt.show()

