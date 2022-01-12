设备： intel i7-6700HQ, ubuntu20.10(with GCC-9)

origin time consume
---

```
Before:

Initial energy: -13.414000
Final energy: -13.356842
Time elapsed: 1519 ms
-------------------------------------
After:

Initial energy: -13.414012
Final energy: -13.356912
Time elapsed: 174 ms
```

**有用的：**
---
- [x]因爲只有48個star，因此可以考慮unroll。
- [x]修改爲soa格式,使用SIMD。直接struct嵌套数组了，没使用array。
- [x]开启-ffast-math。
- [x]sqrt前+std:: 使用模板完成传入参数匹配。
- [x]constexper，提高了3ms。



**没用上的：**
---
- [ ]别名问题使用__restrict or progma omp simd  . 这里没有函数传参数。
- [ ]除法变乘法 - 这个是不是编译器自己做的呀。
- [ ]删除随机访问。
- [ ]尝试了一下对齐，就添加一个padding到512byte。没有提升。

