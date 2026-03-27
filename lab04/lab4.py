import random
import numpy as np

class CustomRandomGenerator:
    def __init__(self, seed=42):
        self.M = 2**63
        self.beta = 2**32 + 3
        self.x = seed

    def next_double(self):
        self.x = (self.beta * self.x) % self.M
        return self.x / self.M


N = 100000

custom_rnd = CustomRandomGenerator(seed=42)
built_in_rnd = random.Random(42)


custom_sample = [custom_rnd.next_double() for _ in range(N)]
built_in_sample = [built_in_rnd.random() for _ in range(N)]

custom_mean = np.mean(custom_sample)
custom_var = np.var(custom_sample, ddof=1)

built_in_mean = np.mean(built_in_sample)
built_in_var = np.var(built_in_sample, ddof=1)


theory_mean = 0.5
theory_var = 1 / 12


print(f"Размер выборки: {N} значений\n")

print("Теоретические значения")
print(f"Среднее:   {theory_mean:.10f}")
print(f"Дисперсия: {theory_var:.10f}\n")

print("Реализованный датчик (MCG)")
print(f"Выборочное среднее:   {custom_mean:.10f} (Погрешность: {abs(custom_mean - theory_mean):.10f})")
print(f"Выборочная дисперсия: {custom_var:.10f} (Погрешность: {abs(custom_var - theory_var):.10f})\n")

print("Встроенный датчик (Python random)")
print(f"Выборочное среднее:   {built_in_mean:.10f} (Погрешность: {abs(built_in_mean - theory_mean):.10f})")
print(f"Выборочная дисперсия: {built_in_var:.10f} (Погрешность: {abs(built_in_var - theory_var):.10f})")