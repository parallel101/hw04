# Init Version: main.cpp
Initial energy: -13.414000
Final energy: -13.356842
Time elapsed: 2545 ms

# 打开-fast-math和-march=native: opt_main_vec.cpp
使用数值加速，依然使用vector。
Initial energy: -10.980484
Final energy: -10.608335
Time elapsed: 1594 ms

# 使用数组代替vector: opt_main_array.cpp
Initial energy: -13.787754
Final energy: -13.404387
Time elapsed: 333 ms

# 将内部循环体内的大计算转为多个循环体的简单计算: opt_main_seperate.cpp
step的大循环体一次对多个值进行了计算，其实它们是可以分开计算的，这样对cache更友好。
Initial energy: -13.787754
Final energy: -13.404387
Time elapsed: 278 ms

# 对外层循环体进行unroll: opt_main_unroll.cpp
对外层循环体进行unroll。
Initial energy: -13.787754
Final energy: -13.404387
Time elapsed: 215 ms

# 最终结果
2545 / 215 = 11.8倍加速。