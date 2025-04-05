from sympy import symbols, Eq, sin, solve, Sum, pi, Matrix, factorial, sqrt

import time

start_time = time.time()

N = 1
K = 10
b = 0.5
L = 2 * N * b
M = 60.21
EI = 6.38 * 10**6
V = 138.89
P0 = 83385
nsmal = 11
# فرمول‌های جدید
p_bar = 2 * P0 / (M * L)
omega_n = lambda n: (n**2 * pi**2 / L**2) * sqrt(EI / M)  # مقدار \omega_n
omega_bar = lambda n: (n * pi * V) / L  # مقدار \overline{\omega_n}

# تعریف پارامترهای a_{i, 2p-1} برای i از 1 تا 2N-1 و p
a_ij = Matrix([[symbols(f'a_{i}_{2 * p - 1}') for p in range(1, K + 1)] for i in range(1, 2 * N)])

# تعریف متغیر i به صورت نمادین
i = symbols('i')
print(a_ij)
# تابعی برای محاسبه C_n^{2p-1}
def compute_C_n_2p_1(n, p, M, L, N, b):
    return Sum((2 / (M * L)) * a_ij[i - 1, p - 1] * sin(n * i * b * pi / L),
                  (i, 1, 2 * N - 1)).doit()

# متغیر برای ذخیره معادلات
equations = []

# حلقه برای تغییر مقدار k از 1 تا maxk
for k in range(1, 1 + K):
    for l in range(1, 2 * N):
        result_l = 0
        for n in range(1, nsmal+1):  # حلقه n
            for p in range(1, 1 + k):  # حلقه p
                term3 = 0
                # محاسبه C_n^{2p-1}
                C_n_2p_1 = compute_C_n_2p_1(n, p, M, L, N, b)
                # محاسبه قسمت اول معادله
                term1 = ((-1) ** (k + p) * C_n_2p_1 * (omega_n(n)/100) ** (2 * (k - p)) * factorial(2 * p - 1)) / factorial(2 * k + 1)
                term3 += term1

            term3 *= sin(n * pi * l / (2 * N))
            result_l += term3

            # محاسبه قسمت دوم معادله
            term2 = (-1) * (((-1) ** (k + 1)) / factorial(2 * k + 1)) * (
                        (p_bar * omega_bar(n) * ((omega_n(n)/100) ** (2 * k) - omega_bar(n) ** (2 * k))) / (
                        (omega_n(n)/100) ** 2 - omega_bar(n) ** 2))

            term2 *= sin(n * pi * l / (2 * N))
            result_l += term2
        # ساخت معادله برای هر l به فرم معادله = 0
        equations.append(Eq(result_l, 0))
# بررسی تعداد معادلات و متغیرها
print(f"Number of equations: {len(equations)}")
print(f"Number of variables: {a_ij.shape[0] * a_ij.shape[1]}")

# ذخیره معادلات در فایل
with open('../equations.txt', 'w') as file:
    for eq in equations:
        file.write(str(eq) + '\n')
solution = solve(equations, a_ij)
#def limit_precision(expr, digits=4):
   # """محدود کردن دقت عددی به تعداد مشخصی از ارقام"""
    #return expr.evalf(n=digits)

# حل دستگاه معادلات با دقت محدود به 10 رقم
#solution = solve([limit_precision(eq) for eq in equations], a_ij)
# چاپ حل دستگاه خط به خط
with open('../aij_solutions_with_N_K.txt', 'w') as file:
    # ذخیره مقادیر N، K و a_ij در فایل خروجی
    output_file = '../aij_solutions_with_N_K.txt'
    with open(output_file, 'w') as file:
        file.write(f"N = {N}\n")  # ذخیره مقدار N
        file.write(f"K = {K}\n")  # ذخیره مقدار K
        file.write(f"b = {b}\n")  # ذخیره مقدار N
        file.write(f"L = {L}\n")  # ذخیره مقدار K
        file.write(f"V = {V}\n")  # ذخیره مقدار K
        file.write("Solutions for a_ij:\n")
        for sol in solution:
            file.write(f"{sol} = {solution[sol]}\n")

end_time = time.time()
execution_time = end_time - start_time

print(execution_time)
print("مقادیر a_ij و زمان اجرا در فایل 'aij_solutions_of_8.txt' ذخیره شدند.")
