import numpy as np
import matplotlib.pyplot as plt

def hill_equation(tf1_concentration, tf2_concentration, n, m, K):
    numerator = tf1_concentration**n * tf2_concentration**m
    denominator = K**n + tf1_concentration**n
    gene_expression = numerator / denominator
    return gene_expression

def activate_2(A1, A2, Kd1, n1, Kd2=0, n2=0):

    if not Kd2:
        Kd2 = Kd1
    if not n2:
        n2 = n1
     
    return pow(A1/Kd1, n1) * pow(A2/Kd2, n2)/(1 + pow(A1/Kd1, n1) + pow(A2/Kd2, n2) + pow(A1/Kd1, n1) * pow(A2/Kd2, n2))

def activate_1(A, Kd, n):
    return pow(A/Kd, n)/(1 + pow(A/Kd, n))

def simplify_and(A1, A2, Kd, n):
    return activate_1(A1, Kd, n) * activate_1(A2, Kd, n)

def xor_gate(R1, R2, Kd1, n1):
    return activate_1((pow(R1/Kd1, n1)+pow(R2/Kd1, n1))/(1 + pow(R1/Kd1, n1)*pow(R2/Kd1, n1)), Kd1, n1)

# Parameters
m = 2  # Hill coefficient for TF2
n = 2 # Hill coefficient for TF1
K = 500  # Dissociation constant

# Concentrations of TF1 and TF2
tf1_concentration = np.linspace(0, 1000, 100)
tf2_concentration = np.linspace(0, 1000, 100)

# Create a 2D grid of concentrations
tf1_grid, tf2_grid = np.meshgrid(tf1_concentration, tf2_concentration)

# Calculate gene expression for each combination of TF1 and TF2 concentrations
# gene_expression_grid = hill_equation(tf1_grid, tf2_grid, n, m, K)
# gene_expression_grid = activate_2(tf1_grid, tf2_grid, K, n) * 2000
gene_expression_grid = xor_gate(tf1_grid, tf2_grid, K, n) * 20000000




# # Plot the results
# plt.figure(figsize=(8, 6))
# plt.contourf(tf1_grid, tf2_grid, gene_expression_grid, cmap='viridis', levels=20)
# # plt.colorbar(label='Gene Expression')
# plt.xlabel('Koncentracija proteina x')
# plt.ylabel('Koncentracija proteina y')
# plt.show()




fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot in 3D
sc = ax.scatter(tf1_grid, tf2_grid, gene_expression_grid, c=gene_expression_grid, cmap='viridis')

# Set labels and title
ax.set_xlabel('Protein x [nM]') # nM = number of molecules
ax.set_ylabel('Protein y [nM]')
ax.set_zlabel('Protein z [nM]')
plt.show()
