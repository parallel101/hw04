origin time consume
---

```
Initial energy: -13.414000
Final energy: -13.356842
Time elapsed: 1519 ms
```
- [x]因爲只有48個star，因此可以考慮unroll
- [x] 表格中速度中最快的是soa_unroll,修改爲soa格式的速度应该会快一点。是因为SIMD吧
- [x]别名问题使用__restrict or progma omp simd  . 这里没有函数传参数。
- [x]除法变乘法 - 这个是不是编译器自己做的呀。
- [x]-ffast-math
- [x]sqrt前+std::
- [x]可不可以改成array --> 我直接没用array直接struct嵌套数组了。
- constexper?
- [x]随机访问？ none
- [x]对齐？
