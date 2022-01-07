# hw4

最终的加速比=1030/104=9.9

使用的方法为 -fast-math + -march=native + SOA + 提取不变量

## 原始时间
```
Initial energy: -13.414000
Final energy: -13.356842
Time elapsed: 1030 ms
```

## 增加 -fast-math
```
Initial energy: -13.413998
Final energy: -13.392558
Time elapsed: 333 ms
```
部分汇编代码，可以看到汇编代码使用的指令是并行指令集，使用的寄存器是ymm
```x86asm
vinsertf128	$1, %xmm4, %ymm0, %ymm10
vbroadcastss	196(%rsi), %ymm1
vmovaps	%ymm21, %ymm12
vpermt2ps	%ymm30, %ymm5, %ymm12
vmovups	304(%rsi), %xmm24
vblendps	$24, %ymm3, %ymm0, %ymm0        # ymm0 = ymm0[0,1,2],ymm3[3,4],ymm0[5,6,7]
vunpcklpd	%ymm24, %ymm22, %ymm3   # ymm3 = ymm22[0],ymm24[0],ymm22[2],ymm24[2]
vmovups	384(%rsi), %xmm5

```

## 增加 SOA结构
```
# 原始代码，只增加了-fast-math和-march=native编译选项
Initial energy: -13.413998
Final energy: -13.392558
Time elapsed: 324 ms

# SOA代码
Initial energy: -13.413998
Final energy: -13.388204
Time elapsed: 113 ms
```

## 使用数组和std::array\<float>
```
Initial energy: -13.413998
Final energy: -13.392558
Time elapsed: 322 ms

Initial energy: -13.413998
Final energy: -13.388204
Time elapsed: 116 ms
```
运行时间基本没有区别

## 给不变量加括号
```
Initial energy: -13.413998
Final energy: -13.392558
Time elapsed: 329 ms

Initial energy: -13.413998
Final energy: -13.388204
Time elapsed: 104 ms
```
提升了大概10ms
