import matplotlib.pyplot as plt
import numpy as np

max_volumes = np.array([
    0.0416667,
    0.00294767,
    0.000527096,
    3.36425e-05
])

average_volumes = np.array([
    0.0416667,
    0.00114025,
    0.000202675,
    1.36218e-05
])

l2_res = np.array([
    1.04427,
    0.101343,
    0.0284829,
    0.00427462
])

plt.figure(figsize=(7, 5))
plt.loglog(max_volumes, l2_res, 'bo-', label='Measured Error')

ref_line_max = [l2_res[0] * (v / max_volumes[0])**(2/3) for v in max_volumes]
plt.loglog(max_volumes, ref_line_max, 'r--', label='O(2/3)')

plt.xlabel('Max Volume of Elements')
plt.ylabel('L2 Residual')
plt.title('Convergence: L2 Residual vs Max Volume')
plt.grid(True, which="both", ls="-", alpha=0.4)
plt.legend()
plt.show()

plt.figure(figsize=(7, 5))
plt.loglog(average_volumes, l2_res, 'bo-', label='Measured Error')

ref_line_avg = [l2_res[0] * (v / average_volumes[0])**(2/3) for v in average_volumes]
plt.loglog(average_volumes, ref_line_avg, 'r--', label='O(2/3)')

plt.xlabel('Average Volume of Elements')
plt.ylabel('L2 Residual')
plt.title('Convergence: L2 Residual vs Average Volume')
plt.grid(True, which="both", ls="-", alpha=0.4)
plt.legend()
plt.show()
